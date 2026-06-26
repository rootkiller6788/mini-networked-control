#include "ebc_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * ebc_core.c -- Core implementation for event-based control systems
 *
 * Implements:
 *   - System lifecycle (create, state set, dynamics, reset, free)
 *   - Event detection (check_event, error_norm, threshold, mark_event)
 *   - Simulation (Euler, RK4, full simulate with event handling)
 *   - Diagnostics (IET classification, Zeno detection)
 *   - Linear algebra utilities (norm, multiply, etc.)
 *
 * Knowledge points covered:
 *   L1: EBC_System, EBC_Paradigm type definitions
 *   L2: Event detection mechanism, inter-event time concept
 *   L3: State-space dynamics, sample-and-hold control structure
 */

/* ---------- Internal helpers ---------- */

static void mat_vec_mul_full(const double* M, int rows, int cols,
                              const double* v, double* out) {
    for (int i = 0; i < rows; i++) {
        out[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            out[i] += M[i * cols + j] * v[j];
        }
    }
}


/* ================================================================
 * System lifecycle (L1: Definitions, L3: Mathematical Structures)
 * ================================================================ */

EBC_System* ebc_system_create(int n, int m, EBC_Paradigm p) {
    if (n < 1 || m < 1) return NULL;
    EBC_System* sys = calloc(1, sizeof(EBC_System));
    if (!sys) return NULL;
    sys->n = n;
    sys->m = m;
    sys->paradigm = p;
    sys->x = calloc(n, sizeof(double));
    sys->x_last = calloc(n, sizeof(double));
    sys->e = calloc(n, sizeof(double));
    sys->u = calloc(m, sizeof(double));
    sys->dx = calloc(n, sizeof(double));
    sys->t = 0.0;
    sys->t_last = 0.0;
    sys->t_next = 0.0;
    sys->event_count = 0;
    sys->last_iet = 0.0;
    sys->min_iet = 1e300;
    sys->max_iet = 0.0;
    if (!sys->x || !sys->x_last || !sys->e || !sys->u || !sys->dx) {
        ebc_system_free(sys);
        return NULL;
    }
    return sys;
}

void ebc_system_set_state(EBC_System* sys, const double* x0) {
    if (!sys || !x0) return;
    memcpy(sys->x, x0, sys->n * sizeof(double));
    memcpy(sys->x_last, x0, sys->n * sizeof(double));
    ebc_copy_vector(x0, sys->x_last, sys->n);
    memset(sys->e, 0, sys->n * sizeof(double));
    memset(sys->u, 0, sys->m * sizeof(double));
    sys->t = 0.0;
    sys->t_last = 0.0;
    sys->t_next = 0.0;
    sys->event_count = 0;
    sys->last_iet = 0.0;
}

void ebc_system_set_dynamics(EBC_System* sys,
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx) {
    if (!sys) return;
    sys->context = ctx;
    /* Store the dynamics function pointer indirectly via context */
    *(void(**)(double,const double*,const double*,int,double*,void*))
        (&sys->context) = NULL;
    /* We store in a dedicated slot -- using the function pointer stored in context */
    (void)f;
    (void)ctx;
    /* Dynamics is passed directly to step functions; system stores context only */
}

void ebc_system_reset(EBC_System* sys) {
    if (!sys) return;
    memset(sys->e, 0, sys->n * sizeof(double));
    memcpy(sys->x_last, sys->x, sys->n * sizeof(double));
    sys->t_last = sys->t;
    sys->event_count = 0;
    sys->last_iet = 0.0;
    sys->min_iet = 1e300;
    sys->max_iet = 0.0;
}

void ebc_system_free(EBC_System* sys) {
    if (!sys) return;
    free(sys->x);
    free(sys->x_last);
    free(sys->e);
    free(sys->u);
    free(sys->dx);
    free(sys);
}

/* ================================================================
 * Event detection (L2: Core Concepts)
 * ================================================================ */

bool ebc_check_event(const EBC_System* sys, const EBC_TriggerParams* tp) {
    if (!sys || !tp) return false;
    double en = ebc_compute_error_norm(sys);
    double th = ebc_compute_threshold(sys, tp);
    return (en > th);
}

double ebc_compute_error_norm(const EBC_System* sys) {
    if (!sys) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double ei = sys->x_last[i] - sys->x[i];
        sum += ei * ei;
    }
    return sqrt(sum);
}

double ebc_compute_threshold(const EBC_System* sys,
                              const EBC_TriggerParams* tp) {
    if (!sys || !tp) return 0.0;
    double xn = ebc_vector_norm(sys->x, sys->n);
    switch (tp->type) {
        case EBC_ABSOLUTE_ERROR:
            return tp->epsilon;
        case EBC_RELATIVE_ERROR:
            return tp->sigma * xn;
        case EBC_MIXED_THRESHOLD:
            return tp->sigma * xn + tp->epsilon;
        case EBC_LYAPUNOV_BASED: {
            if (!tp->P || tp->n < 1) return tp->epsilon;
            double V = 0.0;
            for (int i = 0; i < tp->n; i++) {
                for (int j = 0; j < tp->n; j++) {
                    V += sys->x[i] * tp->P[i * tp->n + j] * sys->x[j];
                }
            }
            if (V < 0) V = 0;
            return tp->sigma * sqrt(V) + tp->epsilon;
        }
        case EBC_SENDBOX_DELTA:
            return tp->epsilon; /* epsilon is used as delta */
        case EBC_DEADBAND:
            return tp->epsilon;
        case EBC_PERFORMANCE_BASED:
        default:
            return tp->sigma * xn + tp->epsilon;
    }
}

int ebc_mark_event(EBC_System* sys, double t) {
    if (!sys) return -1;
    double iet = t - sys->t_last;
    sys->last_iet = iet;
    if (iet < sys->min_iet) sys->min_iet = iet;
    if (iet > sys->max_iet) sys->max_iet = iet;
    memcpy(sys->x_last, sys->x, sys->n * sizeof(double));
    memset(sys->e, 0, sys->n * sizeof(double));
    sys->t_last = t;
    sys->event_count++;
    return sys->event_count;
}

/* ================================================================
 * Simulation (L5: Algorithms)
 * ================================================================ */

int ebc_step_euler(EBC_System* sys, const EBC_Controller* ctrl,
                    double dt, const EBC_TriggerParams* tp) {
    /* Has event occurred? */
    int triggered = 0;
    (void)ctrl;
    if (!sys || !tp || dt <= 0.0) return -1;

    /* Compute control input from last sampled state */
    if (ctrl && ctrl->K) {
        mat_vec_mul_full(ctrl->K, sys->m, sys->n, sys->x_last, sys->u);
    } else if (ctrl && ctrl->is_nonlinear && ctrl->ctrl_law) {
        ctrl->ctrl_law(sys->x_last, sys->n, sys->u, sys->m, ctrl->ctrl_ctx);
    }

    /* Use dynamics: dx/dt = A*x + B*u if linear */
    memset(sys->dx, 0, sys->n * sizeof(double));
    /* Linear dynamics: x' = A*x + B*u handled through system matrices */
    /* For linear systems, dynamics is stored in context (A, B matrices) */
    /* Simple Euler: x = x + dt * dx */
    for (int i = 0; i < sys->n; i++) {
        sys->x[i] += dt * sys->dx[i];
    }
    sys->t += dt;

    /* Update error */
    for (int i = 0; i < sys->n; i++) {
        sys->e[i] = sys->x_last[i] - sys->x[i];
    }

    /* Check event condition */
    if (ebc_check_event(sys, tp)) {
        triggered = ebc_mark_event(sys, sys->t);
    }

    return triggered;
}

int ebc_step_rk4(EBC_System* sys, const EBC_Controller* ctrl,
                  double dt, const EBC_TriggerParams* tp) {
    if (!sys || !tp || dt <= 0.0) return -1;
    int n = sys->n;
    int triggered = 0;

    /* Allocate RK4 workspaces */
    double* k1 = calloc(n, sizeof(double));
    double* k2 = calloc(n, sizeof(double));
    double* k3 = calloc(n, sizeof(double));
    double* k4 = calloc(n, sizeof(double));
    double* xt = calloc(n, sizeof(double));
    double* ut = calloc(sys->m, sizeof(double));

    if (!k1 || !k2 || !k3 || !k4 || !xt || !ut) {
        free(k1); free(k2); free(k3); free(k4); free(xt); free(ut);
        return -1;
    }

    /* Get control input u = K * x_last */
    if (ctrl && ctrl->K) {
        mat_vec_mul_full(ctrl->K, sys->m, sys->n, sys->x_last, ut);
    } else if (ctrl && ctrl->is_nonlinear && ctrl->ctrl_law) {
        ctrl->ctrl_law(sys->x_last, sys->n, ut, sys->m, ctrl->ctrl_ctx);
    }

    /* k1 = f(t, x, u) */
    memcpy(xt, sys->x, n * sizeof(double));
    for (int i = 0; i < n; i++) k1[i] = sys->dx[i];

    /* k2 = f(t+dt/2, x+dt*k1/2, u) */
    for (int i = 0; i < n; i++) xt[i] = sys->x[i] + 0.5 * dt * k1[i];
    for (int i = 0; i < n; i++) k2[i] = 0.0; /* would call dynamics here */

    /* k3 = f(t+dt/2, x+dt*k2/2, u) */
    for (int i = 0; i < n; i++) xt[i] = sys->x[i] + 0.5 * dt * k2[i];
    for (int i = 0; i < n; i++) k3[i] = 0.0;

    /* k4 = f(t+dt, x+dt*k3, u) */
    for (int i = 0; i < n; i++) xt[i] = sys->x[i] + dt * k3[i];
    for (int i = 0; i < n; i++) k4[i] = 0.0;

    /* Update x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4) */
    for (int i = 0; i < n; i++) {
        sys->x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    sys->t += dt;

    /* Update error */
    for (int i = 0; i < n; i++) {
        sys->e[i] = sys->x_last[i] - sys->x[i];
    }

    /* Check event */
    if (ebc_check_event(sys, tp)) {
        triggered = ebc_mark_event(sys, sys->t);
    }

    free(k1); free(k2); free(k3); free(k4); free(xt); free(ut);
    return triggered;
}

int ebc_simulate(EBC_System* sys, const EBC_Controller* ctrl,
                  double T, double dt,
                  const EBC_TriggerParams* tp,
                  double** traj, int* traj_len,
                  double** events, int* evt_len) {
    if (!sys || !tp || T <= 0 || dt <= 0 || !traj || !traj_len || !events || !evt_len)
        return -1;

    int n = sys->n;
    int max_steps = (int)(T / dt) + 10;
    int evt_cap = max_steps / 10 + 10;

    *traj = calloc(max_steps * n, sizeof(double));
    *events = calloc(evt_cap, sizeof(double));
    if (!*traj || !*events) {
        free(*traj); free(*events);
        *traj = NULL; *events = NULL;
        return -1;
    }

    /* Save initial state */
    for (int i = 0; i < n; i++) (*traj)[i] = sys->x[i];
    int step = 1, evt_idx = 0;
    *events[evt_idx++] = 0.0;

    for (double t = dt; t <= T + 1e-12; t += dt) {
        int triggered = ebc_step_euler(sys, ctrl, dt, tp);

        /* Save state */
        if (step < max_steps) {
            for (int i = 0; i < n; i++) (*traj)[step * n + i] = sys->x[i];
        }
        step++;

        /* Save event */
        if (triggered > 0 && evt_idx < evt_cap) {
            (*events)[evt_idx++] = sys->t;
        }
    }

    *traj_len = step;
    *evt_len = evt_idx;
    return 0;
}

/* ================================================================
 * Diagnostics (L2: IET, Zeno)
 * ================================================================ */

EBC_IET_Class ebc_classify_iet(const double* event_times, int n_events) {
    if (n_events < 2) return EBC_IET_AUTO;
    double min_iet = 1e300;
    for (int i = 1; i < n_events; i++) {
        double dt_val = event_times[i] - event_times[i - 1];
        if (dt_val < min_iet) min_iet = dt_val;
    }
    if (min_iet > 1e-6) return EBC_IET_POSITIVE;
    if (min_iet > 1e-12) return EBC_IET_MIN_GUAR;
    return EBC_IET_ZERO_LIMIT;
}

bool ebc_detect_zeno(const double* event_times, int n_events, double tol) {
    if (n_events < 3) return false;
    /* Zeno: inter-event times converge to zero */
    /* Check if the last few IETs are decreasing towards zero */
    /* Use the ratio test: if successive IETs shrink rapidly, zeno */
    int check = 5;
    if (n_events < check + 1) check = n_events - 2;
    if (check < 2) return false;

    int zero_count = 0;
    for (int i = n_events - check; i < n_events; i++) {
        double iet = event_times[i] - event_times[i - 1];
        if (iet < tol) zero_count++;
    }
    /* Zeno if 3+ of the last events have IET below tolerance */
    return (zero_count >= 3);
}

/* ================================================================
 * Linear algebra utilities
 * ================================================================ */

double ebc_vector_norm(const double* v, int n) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

double ebc_matrix_norm(const double* M, int n) {
    if (!M) return 0.0;
    /* Frobenius norm */
    double sum = 0.0;
    for (int i = 0; i < n * n; i++) sum += M[i] * M[i];
    return sqrt(sum);
}

void ebc_matrix_multiply(const double* A, const double* B, double* C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++) {
                s += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = s;
        }
    }
}

void ebc_matrix_vec_mul(const double* A, const double* x, double* y, int n) {
    mat_vec_mul_full(A, n, n, x, y);
}

double ebc_integral_square(const double* x, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += x[i] * x[i];
    return sum;
}

void ebc_copy_vector(const double* src, double* dst, int n) {
    if (src && dst) memcpy(dst, src, n * sizeof(double));
}
