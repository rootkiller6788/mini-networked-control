#include "ebc_trigger.h"
#include "ebc_core.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * ebc_trigger.c -- Event-triggering condition implementations (L5)
 *
 * Implements various event-triggering conditions:
 *   - Absolute, relative, mixed threshold
 *   - Lyapunov-based triggering
 *   - Send-on-delta
 *   - Hysteresis triggering
 *   - Dynamic event-triggering (Girard 2015)
 *   - Trigger margin analysis and sigma optimization
 *
 * Knowledge points:
 *   L2: Event condition, threshold, inter-event time
 *   L3: Quadratic forms, Lyapunov functions
 *   L5: Trigger evaluation algorithms
 *   L8: Dynamic triggering (advanced)
 */

/* ---------- Trigger parameter constructors ---------- */

EBC_TriggerParams ebc_trigger_default(void) {
    EBC_TriggerParams tp;
    memset(&tp, 0, sizeof(tp));
    tp.type = EBC_MIXED_THRESHOLD;
    tp.sigma = 0.1;
    tp.epsilon = 0.01;
    tp.eta = 0.0;
    tp.P = NULL;
    tp.n = 0;
    return tp;
}

EBC_TriggerParams ebc_trigger_make_relative(double sigma) {
    EBC_TriggerParams tp;
    memset(&tp, 0, sizeof(tp));
    tp.type = EBC_RELATIVE_ERROR;
    tp.sigma = (sigma > 0 && sigma < 1) ? sigma : 0.1;
    tp.epsilon = 0.0;
    return tp;
}

EBC_TriggerParams ebc_trigger_mixed(double sigma, double epsilon) {
    EBC_TriggerParams tp;
    memset(&tp, 0, sizeof(tp));
    tp.type = EBC_MIXED_THRESHOLD;
    tp.sigma = (sigma > 0 && sigma < 1) ? sigma : 0.1;
    tp.epsilon = (epsilon > 0) ? epsilon : 0.01;
    return tp;
}

EBC_TriggerParams ebc_trigger_lyapunov(double sigma, const double* P, int n) {
    EBC_TriggerParams tp;
    memset(&tp, 0, sizeof(tp));
    tp.type = EBC_LYAPUNOV_BASED;
    tp.sigma = (sigma > 0) ? sigma : 0.1;
    tp.n = n;
    if (P && n > 0) {
        tp.P = malloc(n * n * sizeof(double));
        if (tp.P) memcpy(tp.P, P, n * n * sizeof(double));
    }
    tp.epsilon = 0.001;
    return tp;
}

EBC_TriggerParams ebc_trigger_make_send_on_delta(double delta, int n) {
    EBC_TriggerParams tp;
    memset(&tp, 0, sizeof(tp));
    tp.type = EBC_SENDBOX_DELTA;
    tp.epsilon = (delta > 0) ? delta : 0.1;
    tp.n = n;
    return tp;
}

void ebc_trigger_free(EBC_TriggerParams* tp) {
    if (tp && tp->P) {
        free(tp->P);
        tp->P = NULL;
    }
}

/* ---------- Trigger evaluation dispatcher ---------- */

bool ebc_trigger_evaluate(const EBC_System* sys, const EBC_TriggerParams* tp) {
    if (!sys || !tp) return false;
    if (tp->ext_gamma) {
        double g = tp->ext_gamma(sys->t, sys->x, sys->e, sys->n, tp->ctx);
        return (g >= 0.0);
    }
    switch (tp->type) {
        case EBC_ABSOLUTE_ERROR:  return ebc_trigger_absolute(sys, tp);
        case EBC_RELATIVE_ERROR:  return ebc_trigger_relative(sys, tp);
        case EBC_MIXED_THRESHOLD: return ebc_trigger_mixed_threshold(sys, tp);
        case EBC_LYAPUNOV_BASED:  return ebc_trigger_lyapunov_based(sys, tp);
        case EBC_SENDBOX_DELTA:   return ebc_trigger_send_on_delta(sys, tp);
        case EBC_DEADBAND:        return ebc_trigger_absolute(sys, tp);
        case EBC_PERFORMANCE_BASED:
        default:                  return ebc_trigger_mixed_threshold(sys, tp);
    }
}

/* ---------- Specific trigger conditions ---------- */

bool ebc_trigger_absolute(const EBC_System* sys, const EBC_TriggerParams* tp) {
    if (!sys || !tp) return false;
    double err_norm = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double ei = sys->x_last[i] - sys->x[i];
        err_norm += ei * ei;
    }
    err_norm = sqrt(err_norm);
    return (err_norm > tp->epsilon);
}

bool ebc_trigger_relative(const EBC_System* sys, const EBC_TriggerParams* tp) {
    if (!sys || !tp) return false;
    double err_norm = 0.0, state_norm = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double ei = sys->x_last[i] - sys->x[i];
        err_norm += ei * ei;
        state_norm += sys->x[i] * sys->x[i];
    }
    err_norm = sqrt(err_norm);
    state_norm = sqrt(state_norm);
    if (state_norm < 1e-12) state_norm = 1e-12;
    return (err_norm > tp->sigma * state_norm);
}

bool ebc_trigger_mixed_threshold(const EBC_System* sys,
                                  const EBC_TriggerParams* tp) {
    if (!sys || !tp) return false;
    double err_norm = 0.0, state_norm = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double ei = sys->x_last[i] - sys->x[i];
        err_norm += ei * ei;
        state_norm += sys->x[i] * sys->x[i];
    }
    err_norm = sqrt(err_norm);
    state_norm = sqrt(state_norm);
    double threshold = tp->sigma * state_norm + tp->epsilon;
    return (err_norm > threshold);
}

bool ebc_trigger_lyapunov_based(const EBC_System* sys,
                                 const EBC_TriggerParams* tp) {
    if (!sys || !tp || !tp->P || tp->n < 1) return false;
    int n = tp->n;
    /* Compute V(x) = x' P x */
    double V = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V += sys->x[i] * tp->P[i * n + j] * sys->x[j];
        }
    }
    if (V < 0.0) V = 0.0;
    double err_norm = ebc_compute_error_norm(sys);
    double threshold = tp->sigma * sqrt(V) + tp->epsilon;
    return (err_norm > threshold);
}

bool ebc_trigger_send_on_delta(const EBC_System* sys,
                                const EBC_TriggerParams* tp) {
    if (!sys || !tp) return false;
    double delta_sq = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double d = sys->x[i] - sys->x_last[i];
        delta_sq += d * d;
    }
    return (sqrt(delta_sq) > tp->epsilon);
}

bool ebc_trigger_hysteresis(const EBC_System* sys,
                             double sigma_upper, double sigma_lower,
                             double epsilon) {
    if (!sys) return false;
    double err_norm = 0.0, state_norm = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double ei = sys->x_last[i] - sys->x[i];
        err_norm += ei * ei;
        state_norm += sys->x[i] * sys->x[i];
    }
    err_norm = sqrt(err_norm);
    state_norm = sqrt(state_norm);
    if (state_norm < 1e-12) state_norm = 1e-12;
    double th_upper = sigma_upper * state_norm + epsilon;
    double th_lower = sigma_lower * state_norm + epsilon;
    /* If error exceeds upper, fire event; if below lower, reset */
    /* Use the error norm vs upper threshold as the decision */
    /* The state machine would track whether we are in "armed" state */
    return (err_norm > th_upper);
}

/* ---------- Dynamic event-triggering (Girard 2015) ---------- */

EBC_DynamicTrigger ebc_dynamic_trigger_create(double eta0, double beta,
                                               double theta, double sigma) {
    EBC_DynamicTrigger dt;
    dt.eta = (eta0 > 0) ? eta0 : 1.0;
    dt.eta0 = dt.eta;
    dt.beta = (beta > 0) ? beta : 1.0;
    dt.theta = (theta > 0) ? theta : 1.0;
    dt.sigma = (sigma > 0 && sigma < 1) ? sigma : 0.1;
    return dt;
}

void ebc_dynamic_trigger_update(EBC_DynamicTrigger* dt,
                                 const double* x, const double* e,
                                 int n, double h) {
    if (!dt || !x || !e || n < 1) return;
    /* Compute |x| and |e| */
    double xn = 0.0, en = 0.0;
    for (int i = 0; i < n; i++) {
        xn += x[i] * x[i];
        en += e[i] * e[i];
    }
    xn = sqrt(xn);
    en = sqrt(en);
    /* d(eta)/dt = -beta*eta + (sigma*|x| - |e|) */
    double deta = -dt->beta * dt->eta + (dt->sigma * xn - en);
    /* Euler update */
    dt->eta += h * deta;
    if (dt->eta < 0.0) dt->eta = 0.0;
}

bool ebc_dynamic_trigger_evaluate(const EBC_DynamicTrigger* dt,
                                   const double* x, const double* e, int n) {
    if (!dt || !x || !e) return false;
    double xn = 0.0, en = 0.0;
    for (int i = 0; i < n; i++) {
        xn += x[i] * x[i];
        en += e[i] * e[i];
    }
    xn = sqrt(xn);
    en = sqrt(en);
    double margin = dt->sigma * xn - en;
    /* Event fires when eta + theta * margin <= 0 */
    return (dt->eta + dt->theta * margin <= 0.0);
}

/* ---------- Trigger margin analysis ---------- */

void ebc_trigger_margin_trace(const double* x_traj,
                               const double* e_traj,
                               int n, int len, double sigma, double epsilon,
                               double* margins) {
    if (!x_traj || !e_traj || !margins || n < 1 || len < 1) return;
    for (int k = 0; k < len; k++) {
        double xn = 0.0, en = 0.0;
        for (int i = 0; i < n; i++) {
            xn += x_traj[k * n + i] * x_traj[k * n + i];
            en += e_traj[k * n + i] * e_traj[k * n + i];
        }
        xn = sqrt(xn);
        en = sqrt(en);
        margins[k] = sigma * xn + epsilon - en;
    }
}

/* ---------- Sigma optimization ---------- */

double ebc_trigger_optimize_sigma(
    void (*system_dynamics)(double, const double*, const double*, int,
                             double*, void*),
    void* ctx, int n, int m,
    const double* K, const double* x0,
    double T, double dt,
    double target_iet, double epsilon) {
    if (!system_dynamics || n < 1 || m < 1 || !K || !x0 || target_iet <= 0)
        return 0.1;

    /* Binary search for sigma that achieves target IET */
    double lo = 0.001, hi = 0.999;
    int max_iter = 30;
    double best_sigma = 0.1, best_diff = 1e300;

    for (int iter = 0; iter < max_iter; iter++) {
        double mid = (lo + hi) / 2.0;

        /* Simulate with this sigma and measure avg IET */
        EBC_System* sys = ebc_system_create(n, m, EBC_CONTINUOUS_ETC);
        if (!sys) break;

        ebc_system_set_state(sys, x0);
        EBC_TriggerParams tp = ebc_trigger_mixed(mid, epsilon);

        /* Simulate with simple Euler */
        int n_steps = (int)(T / dt);
        int events_count = 0;
        double last_event_t = 0.0;

        /* Local state buffer */
        double* x = malloc(n * sizeof(double));
        double* u = malloc(m * sizeof(double));
        double* dx = malloc(n * sizeof(double));
        double* x_last = malloc(n * sizeof(double));
        if (x && u && dx && x_last) {
            memcpy(x, x0, n * sizeof(double));
            memcpy(x_last, x0, n * sizeof(double));
            for (int step = 0; step < n_steps; step++) {
                double t = step * dt;
                /* Compute control */
                for (int i = 0; i < m; i++) {
                    u[i] = 0.0;
                    for (int j = 0; j < n; j++)
                        u[i] += K[i * n + j] * x_last[j];
                }
                /* Compute dynamics */
                system_dynamics(t, x, u, n, dx, ctx);
                /* Euler step */
                for (int i = 0; i < n; i++) x[i] += dt * dx[i];
                /* Check event */
                double en = 0.0, xn = 0.0;
                for (int i = 0; i < n; i++) {
                    double ei = x_last[i] - x[i];
                    en += ei * ei;
                    xn += x[i] * x[i];
                }
                en = sqrt(en);
                xn = sqrt(xn);
                if (en > mid * xn + epsilon) {
                    memcpy(x_last, x, n * sizeof(double));
                    events_count++;
                    last_event_t = t;
                }
            }
        }
        free(x); free(u); free(dx); free(x_last);
        ebc_system_free(sys);

        double avg_iet = (events_count > 0) ? T / events_count : T;
        double diff = fabs(avg_iet - target_iet);

        if (diff < best_diff) {
            best_diff = diff;
            best_sigma = mid;
        }

        if (avg_iet > target_iet) {
            hi = mid;  /* Reduce sigma to decrease IET */
        } else {
            lo = mid;  /* Increase sigma to increase IET */
        }

        if (hi - lo < 0.001) break;
    }

    return best_sigma;
}
