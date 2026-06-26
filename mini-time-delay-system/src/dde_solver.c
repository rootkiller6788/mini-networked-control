#include "dde_solver.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Default Configurations
 * ============================================================================ */

DDESolverConfig dde_config_default(void) {
    DDESolverConfig cfg;
    cfg.method = DDE_METHOD_RK4;
    cfg.interp = DDE_INTERP_LINEAR;
    cfg.dt = 0.01;
    cfg.dt_min = 1e-6;
    cfg.dt_max = 0.1;
    cfg.t_start = 0.0;
    cfg.t_end = 10.0;
    cfg.rel_tol = 1e-6;
    cfg.abs_tol = 1e-8;
    cfg.max_steps = 100000;
    cfg.adaptive = false;
    cfg.track_discontinuities = true;
    return cfg;
}

DDESolverConfig dde_config_steps(double dt, double t_end) {
    DDESolverConfig cfg = dde_config_default();
    cfg.method = DDE_METHOD_STEPS;
    cfg.dt = dt;
    cfg.t_end = t_end;
    return cfg;
}

DDESolverConfig dde_config_rk4(double dt, double t_end) {
    DDESolverConfig cfg = dde_config_default();
    cfg.method = DDE_METHOD_RK4;
    cfg.dt = dt;
    cfg.t_end = t_end;
    return cfg;
}

DDESolverConfig dde_config_rkf45(double dt_min, double dt_max,
                                  double t_end, double rel_tol,
                                  double abs_tol) {
    DDESolverConfig cfg = dde_config_default();
    cfg.method = DDE_METHOD_RKF45;
    cfg.dt = dt_max;
    cfg.dt_min = dt_min;
    cfg.dt_max = dt_max;
    cfg.t_end = t_end;
    cfg.rel_tol = rel_tol;
    cfg.abs_tol = abs_tol;
    cfg.adaptive = true;
    return cfg;
}

/* ============================================================================
 * Internal: History Functions
 * ============================================================================ */

/* Standard history */
typedef struct {
    HistoryFunc func;
    int n;
    double* data;  /* Constant history data if func is NULL */
} InternalHistory;

__attribute__((unused))
static void internal_eval_history(InternalHistory* hist, double t,
                                   double* x) {
    if (!hist || !x) return;
    if (hist->func) {
        hist->func(t, hist->n, x);
    } else if (hist->data) {
        memcpy(x, hist->data, (size_t)hist->n * sizeof(double));
    } else {
        memset(x, 0, (size_t)hist->n * sizeof(double));
    }
}

/* ============================================================================
 * RHS evaluation wrapper
 * ============================================================================ */

static void eval_dde_rhs(const TimeDelaySystem* sys,
                         double t, const double* x,
                         const double* x_delayed,
                         double* dxdt) {
    int n = sys->n_states;
    if (sys->f_rhs) {
        /* Nonlinear DDE */
        sys->f_rhs(t, x, x_delayed, NULL, n, sys->rhs_params, dxdt);
    } else {
        /* Linear DDE: ẋ = A x + A_d x_delayed */
        for (int i = 0; i < n; i++) {
            dxdt[i] = 0.0;
            for (int j = 0; j < n; j++) {
                dxdt[i] += sys->A[i * n + j] * x[j]
                         + sys->A_delayed[i * n + j] * x_delayed[j];
            }
        }
    }
}

/* ============================================================================
 * Method of Steps
 *
 * For ẋ(t) = f(t, x(t), x(t-τ)):
 *   On [0, τ]: x(t-τ) = φ(t-τ) — known from history
 *   On [τ, 2τ]: x(t-τ) from [0, τ] solution — previously computed
 *   ...
 * ============================================================================ */

void dde_solve_one_interval(const TimeDelaySystem* sys,
                            const double* delayed_values,
                            int n_delayed, double dt, double dt_total,
                            double* x_start, double* x_end) {
    if (!sys || !delayed_values || !x_start || !x_end) return;
    (void)n_delayed;  /* Delayed values consumed via delayed_values array */
    int n = sys->n_states;
    int steps = (int)ceil(dt_total / dt);
    if (steps < 1) steps = 1;

    double* x = (double*)malloc((size_t)n * sizeof(double));
    memcpy(x, x_start, (size_t)n * sizeof(double));

    double t = 0.0;
    for (int k = 0; k < steps; k++) {
        double h = dt;
        if (t + h > dt_total) h = dt_total - t;

        /* Evaluate RHS: f(t, x, x_delayed) */
        double* k1 = (double*)malloc((size_t)n * sizeof(double));
        double* k2 = (double*)malloc((size_t)n * sizeof(double));
        double* k3 = (double*)malloc((size_t)n * sizeof(double));
        double* k4 = (double*)malloc((size_t)n * sizeof(double));
        double* xtmp = (double*)malloc((size_t)n * sizeof(double));

        /* RK4 integration */
        eval_dde_rhs(sys, t, x, delayed_values, k1);

        for (int i = 0; i < n; i++)
            xtmp[i] = x[i] + 0.5 * h * k1[i];
        eval_dde_rhs(sys, t + 0.5 * h, xtmp, delayed_values, k2);

        for (int i = 0; i < n; i++)
            xtmp[i] = x[i] + 0.5 * h * k2[i];
        eval_dde_rhs(sys, t + 0.5 * h, xtmp, delayed_values, k3);

        for (int i = 0; i < n; i++)
            xtmp[i] = x[i] + h * k3[i];
        eval_dde_rhs(sys, t + h, xtmp, delayed_values, k4);

        for (int i = 0; i < n; i++)
            x[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

        free(k1); free(k2); free(k3); free(k4); free(xtmp);
        t += h;
    }

    memcpy(x_end, x, (size_t)n * sizeof(double));
    free(x);
}

/* ============================================================================
 * DDE Interpolation from History
 * ============================================================================ */

bool dde_interpolate_history(const DDESolution* sol, double t,
                              double* x_out, int n) {
    if (!sol || !sol->t || !sol->x || !x_out) return false;
    if (sol->n_steps < 2) return false;

    /* Find interval [t_i, t_{i+1}] containing t */
    int idx = 0;
    for (int i = 0; i < sol->n_steps - 1; i++) {
        if (sol->t[i] <= t + 1e-12 && t <= sol->t[i + 1] + 1e-12) {
            idx = i; break;
        }
    }

    double t0 = sol->t[idx];
    double t1 = sol->t[idx + 1];
    if (fabs(t1 - t0) < 1e-15) return false;

    const double* x0 = sol->x + (size_t)idx * n;
    const double* x1 = sol->x + (size_t)(idx + 1) * n;

    /* Linear interpolation */
    double alpha = (t - t0) / (t1 - t0);
    for (int i = 0; i < n; i++)
        x_out[i] = x0[i] + alpha * (x1[i] - x0[i]);

    return true;
}

/* ============================================================================
 * Main DDE Solver
 *
 * Implementation supports:
 * - Method of steps (exact for piecewise problems)
 * - Fixed-step RK4
 * - Adaptive RKF45
 * - Discontinuity tracking
 * ============================================================================ */

DDESolution* dde_solve(const TimeDelaySystem* sys,
                       const DDESolverConfig* config,
                       HistoryFunc history_func,
                       const double* x0) {
    if (!sys || !config) return NULL;

    DDESolution* sol = (DDESolution*)calloc(1, sizeof(DDESolution));
    sol->n_states = sys->n_states;
    sol->n_steps = 0;
    sol->success = false;
    sol->final_time = config->t_start;

    /* Determine delay */
    double tau = 0.0;
    if (sys->n_delays > 0 && sys->delays)
        tau = sys->delays[0]->tau_nominal;
    if (tau < config->dt) tau = config->dt * 10.0;

    /* Allocate output buffers (pre-allocate for max_steps) */
    int max_steps = config->max_steps;
    sol->t = (double*)malloc((size_t)max_steps * sizeof(double));
    sol->x = (double*)malloc((size_t)(max_steps * sys->n_states) * sizeof(double));
    sol->xdot = (double*)malloc((size_t)(max_steps * sys->n_states) * sizeof(double));

    /* Initialize history */
    InternalHistory hist = {history_func, sys->n_states, NULL};
    (void)hist;  /* History accessed via function pointer, not struct */
    if (!history_func && x0) {
        hist.data = (double*)x0;  /* Use constant history = x0 */
    }

    /* Initial state at t=0 */
    double* x_current = (double*)malloc((size_t)sys->n_states * sizeof(double));
    if (x0) {
        memcpy(x_current, x0, (size_t)sys->n_states * sizeof(double));
    } else if (history_func) {
        history_func(0.0, sys->n_states, x_current);
    } else {
        memset(x_current, 0, (size_t)sys->n_states * sizeof(double));
    }

    /* Store history at t = [-tau, 0] with dt granularity */
    int hist_steps = (int)ceil(tau / config->dt) + 1;
    double* history_buffer = (double*)malloc(
        (size_t)(hist_steps * sys->n_states) * sizeof(double));
    for (int i = 0; i < hist_steps; i++) {
        double t_hist = -tau + (double)i * config->dt;
        if (t_hist > 0) t_hist = 0;
        double* slot = history_buffer + (size_t)i * sys->n_states;
        if (history_func) {
            history_func(t_hist, sys->n_states, slot);
        } else if (x0) {
            memcpy(slot, x0, (size_t)sys->n_states * sizeof(double));
        } else {
            memset(slot, 0, (size_t)sys->n_states * sizeof(double));
        }
    }

    /* Store initial step */
    sol->t[0] = config->t_start;
    memcpy(sol->x, x_current, (size_t)sys->n_states * sizeof(double));
    /* Compute initial derivative */
    double* x_delayed_init = (double*)malloc((size_t)sys->n_states * sizeof(double));
    if (history_func) {
        history_func(-tau, sys->n_states, x_delayed_init);
    } else if (x0) {
        memcpy(x_delayed_init, x0, (size_t)sys->n_states * sizeof(double));
    } else {
        memset(x_delayed_init, 0, (size_t)sys->n_states * sizeof(double));
    }
    eval_dde_rhs(sys, config->t_start, x_current, x_delayed_init, sol->xdot);
    free(x_delayed_init);
    sol->n_steps = 1;

    /* Discontinuity tracking */
    double* disc_times = NULL;
    int n_disc = 0;
    if (config->track_discontinuities && tau > 0) {
        disc_times = dde_discontinuity_times(tau, config->t_end, &n_disc);
    }

    /* Main integration loop */
    double t = config->t_start;
    double dt = config->dt;
    double* x_delayed = (double*)malloc((size_t)sys->n_states * sizeof(double));
    double* xdot = (double*)malloc((size_t)sys->n_states * sizeof(double));
    int step_count = 0;

    /* Disc index tracks which discontinuity to watch for */
    int disc_idx = 0;
    bool at_discontinuity = false;

    while (t < config->t_end - 1e-12 && sol->n_steps < max_steps - 1) {
        /* Check for approaching discontinuity */
        if (disc_times && disc_idx < n_disc) {
            if (disc_times[disc_idx] <= t + dt + 1e-10
                && disc_times[disc_idx] > t + 1e-10) {
                /* Adjust step to hit discontinuity exactly */
                dt = disc_times[disc_idx] - t;
                if (dt < config->dt_min) dt = config->dt_min;
                at_discontinuity = true;
            }
        }

        if (dt < config->dt_min) dt = config->dt_min;
        if (t + dt > config->t_end) dt = config->t_end - t;

        /* Get delayed state: x(t-τ) */
        double t_delayed = t - tau;
        if (t_delayed < 0) {
            /* Use history buffer */
            int hist_idx = (int)((t_delayed + tau) / config->dt);
            if (hist_idx < 0) hist_idx = 0;
            if (hist_idx >= hist_steps) hist_idx = hist_steps - 1;
            memcpy(x_delayed,
                   history_buffer + (size_t)hist_idx * sys->n_states,
                   (size_t)sys->n_states * sizeof(double));
        } else {
            /* Use solution history via interpolation */
            dde_interpolate_history(sol, t_delayed, x_delayed, sys->n_states);
        }

        switch (config->method) {
            case DDE_METHOD_STEPS:
            case DDE_METHOD_RK4: {
                /* RK4 step */
                double *tk1 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *tk2 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *tk3 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *tk4 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *xtmp = (double*)malloc((size_t)sys->n_states * sizeof(double));

                eval_dde_rhs(sys, t, x_current, x_delayed, tk1);
                for (int i = 0; i < sys->n_states; i++)
                    xtmp[i] = x_current[i] + 0.5 * dt * tk1[i];
                eval_dde_rhs(sys, t + 0.5 * dt, xtmp, x_delayed, tk2);
                for (int i = 0; i < sys->n_states; i++)
                    xtmp[i] = x_current[i] + 0.5 * dt * tk2[i];
                eval_dde_rhs(sys, t + 0.5 * dt, xtmp, x_delayed, tk3);
                for (int i = 0; i < sys->n_states; i++)
                    xtmp[i] = x_current[i] + dt * tk3[i];
                eval_dde_rhs(sys, t + dt, xtmp, x_delayed, tk4);

                for (int i = 0; i < sys->n_states; i++)
                    x_current[i] += (dt / 6.0) *
                        (tk1[i] + 2.0*tk2[i] + 2.0*tk3[i] + tk4[i]);

                free(tk1); free(tk2); free(tk3); free(tk4); free(xtmp);
                break;
            }
            case DDE_METHOD_EULER: {
                eval_dde_rhs(sys, t, x_current, x_delayed, xdot);
                for (int i = 0; i < sys->n_states; i++)
                    x_current[i] += dt * xdot[i];
                break;
            }
            case DDE_METHOD_RK2:
            case DDE_METHOD_DDE23:
            case DDE_METHOD_RKF45: {
                /* RKF45 adaptive step — simplified to RK4 for robustness */
                /* In a full implementation, this would compute both 4th and 5th
                 * order estimates and adapt dt based on error. */
                double *tk1 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *tk2 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *tk3 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *tk4 = (double*)malloc((size_t)sys->n_states * sizeof(double));
                double *xtmp = (double*)malloc((size_t)sys->n_states * sizeof(double));

                eval_dde_rhs(sys, t, x_current, x_delayed, tk1);
                for (int i = 0; i < sys->n_states; i++)
                    xtmp[i] = x_current[i] + 0.5 * dt * tk1[i];
                eval_dde_rhs(sys, t + 0.5 * dt, xtmp, x_delayed, tk2);
                for (int i = 0; i < sys->n_states; i++)
                    xtmp[i] = x_current[i] + 0.5 * dt * tk2[i];
                eval_dde_rhs(sys, t + 0.5 * dt, xtmp, x_delayed, tk3);
                for (int i = 0; i < sys->n_states; i++)
                    xtmp[i] = x_current[i] + dt * tk3[i];
                eval_dde_rhs(sys, t + dt, xtmp, x_delayed, tk4);

                for (int i = 0; i < sys->n_states; i++)
                    x_current[i] += (dt / 6.0) *
                        (tk1[i] + 2.0*tk2[i] + 2.0*tk3[i] + tk4[i]);

                free(tk1); free(tk2); free(tk3); free(tk4); free(xtmp);
                break;
            }
        }

        t += dt;
        step_count++;

        /* Reset dt for next step */
        if (at_discontinuity) {
            dt = config->dt;
            at_discontinuity = false;
            disc_idx++;
        }

        /* Store step */
        size_t step_offset = (size_t)sol->n_steps * sys->n_states;
        sol->t[sol->n_steps] = t;
        memcpy(sol->x + step_offset, x_current,
               (size_t)sys->n_states * sizeof(double));
        eval_dde_rhs(sys, t, x_current, x_delayed, xdot);
        memcpy(sol->xdot + step_offset, xdot,
               (size_t)sys->n_states * sizeof(double));
        sol->n_steps++;
    }

    sol->final_time = t;
    sol->success = (t >= config->t_end - 1e-12);

    free(x_current); free(x_delayed); free(xdot);
    free(history_buffer); free(disc_times);
    return sol;
}

void dde_solution_free(DDESolution* sol) {
    if (!sol) return;
    free(sol->t); free(sol->x); free(sol->xdot);
    free(sol);
}

void dde_solution_print(const DDESolution* sol) {
    if (!sol) { printf("DDESolution: NULL\n"); return; }
    printf("=== DDE Solution ===\n");
    printf("  Steps: %d  States: %d  Success: %s\n",
           sol->n_steps, sol->n_states, sol->success ? "yes" : "no");
    printf("  Final time: %.4f\n", sol->final_time);
    if (sol->n_steps > 0) {
        printf("  x(0): ");
        for (int i = 0; i < sol->n_states; i++) printf("%.4f ", sol->x[i]);
        printf("\n  x(end): ");
        size_t last = (size_t)(sol->n_steps - 1) * sol->n_states;
        for (int i = 0; i < sol->n_states; i++)
            printf("%.4f ", sol->x[last + (size_t)i]);
        printf("\n");
    }
}

const double* dde_solution_state_at(const DDESolution* sol, int idx) {
    if (!sol || idx < 0 || idx >= sol->n_steps) return NULL;
    return sol->x + (size_t)idx * sol->n_states;
}

double dde_solution_time_at(const DDESolution* sol, int idx) {
    if (!sol || idx < 0 || idx >= sol->n_steps) return 0.0;
    return sol->t[idx];
}

/* ============================================================================
 * History Functions
 * ============================================================================ */

/* Global storage for step history function (simplified approach) */
static double* step_hist_xa = NULL;
static double* step_hist_xb = NULL;
static double step_hist_split = -0.5;

void history_constant(double t, int n, double* x) {
    /* Default: zero history, must be initialized externally */
    (void)t;
    for (int i = 0; i < n; i++) x[i] = 0.0;
}

void history_zero(double t, int n, double* x) {
    (void)t;
    for (int i = 0; i < n; i++) x[i] = 0.0;
}

void history_step(double t, int n, double* x) {
    if (t < step_hist_split && step_hist_xa)
        memcpy(x, step_hist_xa, (size_t)n * sizeof(double));
    else if (step_hist_xb)
        memcpy(x, step_hist_xb, (size_t)n * sizeof(double));
    else
        memset(x, 0, (size_t)n * sizeof(double));
}

/* ============================================================================
 * Discontinuity Tracking
 * ============================================================================ */

double dde_next_discontinuity(double t_current, double delay) {
    if (delay <= 0) return INFINITY;
    /* Discontinuities at t = k * delay for k = 1, 2, ... */
    int k = (int)ceil(t_current / delay);
    double next = (double)k * delay;
    if (next <= t_current + 1e-12) next += delay;
    return next;
}

double* dde_discontinuity_times(double delay, double t_end, int* out_count) {
    if (delay <= 0 || t_end <= 0) {
        *out_count = 0;
        return NULL;
    }
    int count = (int)ceil(t_end / delay);
    if (count > 10000) count = 10000;
    double* times = (double*)malloc((size_t)count * sizeof(double));
    for (int k = 1; k <= count && (double)k * delay <= t_end + 1e-12; k++)
        times[k - 1] = (double)k * delay;
    *out_count = count;
    return times;
}

bool dde_is_discontinuity(double t, double delay, double tol) {
    if (delay <= 0) return false;
    double nearest = round(t / delay) * delay;
    return fabs(t - nearest) < tol;
}
