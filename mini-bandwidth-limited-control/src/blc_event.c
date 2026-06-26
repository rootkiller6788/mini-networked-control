/**
 * blc_event.c — Event-Triggered and Aperiodic Control
 *
 * Implementation of event-triggered communication strategies
 * for bandwidth-limited control:
 *  - Send-on-Delta (transmit on significant change)
 *  - Lebesgue sampling (transmit on level crossing)
 *  - Event-triggered state feedback (Tabuada, 2007)
 *  - Self-triggered control (pre-compute next transmission)
 *
 * Event-triggered control is complementary to quantization:
 * quantization reduces bits per sample, event-triggering
 * reduces samples per second. Together they minimize total bit rate.
 *
 * Knowledge coverage: L5 (Event-Triggered Algorithms), L2 (Concepts)
 */

#include "blc_core.h"
#include "blc_datarate.h"
#include "blc_event.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * Send-on-Delta Implementation
 *
 * Trigger: ||x - x_last||_norm > δ
 *
 * The threshold δ should be chosen relative to the quantization
 * step Δ to avoid transmitting quantization noise:
 *   δ = k · Δ · sqrt(n)   where k ≥ 1
 *
 * For Chebyshev norm: δ = k · Δ (each component separately)
 * For L2 norm: δ = k · Δ · sqrt(n)
 * ================================================================ */

int blc_sod_init(BLCSendOnDelta* sod, double delta, int norm_type) {
    if (!sod || delta <= 0.0) return -1;

    sod->delta         = delta;
    sod->norm_type     = norm_type;
    sod->error_norm    = 0.0;
    sod->transmissions = 0;
    sod->suppressed    = 0;
    sod->bandwidth_saved = 0.0;
    sod->min_interval   = 0.001;  /** 1ms minimum */
    sod->last_tx_time  = -1.0;
    sod->forced_tx     = false;

    memset(sod->x_last, 0, sizeof(sod->x_last));
    return 0;
}

bool blc_sod_should_transmit(BLCSendOnDelta* sod, const double* x,
                              int n, double current_time, bool force) {
    if (!sod || !x) return true;
    if (force) {
        sod->forced_tx = true;
        return true;
    }

    /** Check minimum interval */
    if (sod->last_tx_time >= 0.0 &&
        current_time - sod->last_tx_time < sod->min_interval) {
        sod->suppressed++;
        return false;
    }

    /** Compute error norm */
    double en = 0.0;
    switch (sod->norm_type) {
        case 1:  /** L1 */
            for (int i = 0; i < n; i++)
                en += fabs(x[i] - sod->x_last[i]);
            break;
        case 2:  /** L2 */
            for (int i = 0; i < n; i++) {
                double e = x[i] - sod->x_last[i];
                en += e * e;
            }
            en = sqrt(en);
            break;
        case 0:  /** Chebyshev (L∞) */
        default:
            for (int i = 0; i < n; i++) {
                double e = fabs(x[i] - sod->x_last[i]);
                if (e > en) en = e;
            }
            break;
    }

    sod->error_norm = en;

    if (en > sod->delta) {
        return true;
    } else {
        sod->suppressed++;
        return false;
    }
}

void blc_sod_transmitted(BLCSendOnDelta* sod, const double* x,
                          double current_time) {
    if (!sod || !x) return;

    for (int i = 0; i < BLC_MAX_STATES; i++) {
        sod->x_last[i] = x[i];
    }
    sod->transmissions++;
    sod->last_tx_time = current_time;
    sod->forced_tx = false;

    /** Update bandwidth saved percentage */
    int total = sod->transmissions + sod->suppressed;
    if (total > 0) {
        sod->bandwidth_saved = 100.0 * (double)sod->suppressed / (double)total;
    }
}

void blc_sod_stats(const BLCSendOnDelta* sod, double* pct_saved,
                    int* tx_count, int* suppressed) {
    if (!sod) return;
    if (pct_saved)  *pct_saved  = sod->bandwidth_saved;
    if (tx_count)   *tx_count   = sod->transmissions;
    if (suppressed) *suppressed = sod->suppressed;
}

double blc_sod_optimal_delta(double quant_step, double tradeoff_factor) {
    /** δ_opt = k · Δ
     *  k ∈ [1, 3]: k=1 gives high fidelity, k=3 saves more bandwidth */
    if (tradeoff_factor < 1.0) tradeoff_factor = 1.0;
    if (tradeoff_factor > 5.0) tradeoff_factor = 5.0;
    return tradeoff_factor * quant_step;
}

/* ================================================================
 * Lebesgue Sampler Implementation
 *
 * Level set: L_i = {x : μ^{i-1} < |x| ≤ μ^i}
 *
 * μ = e^{|λ|·T_desired} gives average sampling period T_desired
 * for a stable first-order system with pole λ.
 * ================================================================ */

int blc_lebesgue_init(BLCLebesgueSampler* ls, double mu,
                       int n_levels, double max_val) {
    if (!ls || mu <= 1.0 || n_levels < 2 || n_levels > 32
        || max_val <= 0.0) return -1;

    ls->mu         = mu;
    ls->n_levels   = n_levels;
    ls->current_level = 0;
    ls->last_level    = 0;
    ls->transmissions = 0;
    ls->avg_interval  = 0.0;
    ls->last_tx_time  = -1.0;

    /** Build level thresholds: L_i = μ^i */
    for (int i = 0; i < n_levels; i++) {
        ls->levels[i] = max_val * pow(mu, (double)(i - n_levels/2));
    }

    return 0;
}

bool blc_lebesgue_check(BLCLebesgueSampler* ls, double x, double time) {
    if (!ls) return false;

    double abs_x = fabs(x);

    /** Find current level */
    int new_level = 0;
    for (int i = 1; i < ls->n_levels; i++) {
        if (abs_x > ls->levels[i]) {
            new_level = i;
        }
    }
    if (abs_x <= ls->levels[0]) new_level = 0;

    ls->current_level = new_level;

    /** Transmit if level changed */
    if (new_level != ls->last_level) {
        double interval = time - ls->last_tx_time;
        if (interval > 0) {
            /** Exponential moving average of interval */
            double alpha = 0.2;
            ls->avg_interval = (1.0 - alpha) * ls->avg_interval + alpha * interval;
        }
        return true;
    }
    return false;
}

int blc_lebesgue_get_level(const BLCLebesgueSampler* ls, double x) {
    if (!ls) return 0;
    double abs_x = fabs(x);
    int level = 0;
    for (int i = 1; i < ls->n_levels; i++) {
        if (abs_x > ls->levels[i]) level = i;
    }
    return level;
}

void blc_lebesgue_stats(const BLCLebesgueSampler* ls,
                         double* avg_interval, int* tx_count) {
    if (!ls) return;
    if (avg_interval) *avg_interval = ls->avg_interval;
    if (tx_count)     *tx_count     = ls->transmissions;
}

/* ================================================================
 * Event-Triggered State Feedback (Tabuada, 2007)
 *
 * Theorem: For an ISS system with Lyapunov function V(x),
 * there exists σ > 0 such that the event-triggered implementation
 * with trigger ||e|| ≤ σ||x|| is asymptotically stable and has
 * a strictly positive minimum inter-event time.
 *
 * σ = (√(λ_min(Q) · (1-β)) - ||PB||) / (||PB||)
 * where β ∈ (0, 1) is a design parameter, P is the Lyapunov matrix,
 * and Q = -(A_c'P + PA_c) is positive definite.
 *
 * Minimum inter-event time:
 *   τ_min = (1/λ) · ln(1 + σ·β)
 * ================================================================ */

int blc_etf_init(BLCEventTriggeredFeedback* etf, double sigma,
                  const double* Ac, int n) {
    if (!etf || sigma <= 0.0 || sigma >= 1.0 || !Ac || n < 1) return -1;

    etf->sigma           = sigma;
    etf->beta            = 1.0;
    etf->lambda_cl       = 1.0;
    etf->min_inter_event  = 0.0;
    etf->last_event_time  = -1.0;
    etf->event_count      = 0;
    etf->V_last           = 0.0;
    etf->is_iss           = true;

    memset(etf->x_last, 0, sizeof(etf->x_last));
    memset(etf->e, 0, sizeof(etf->e));

    /** Estimate λ_cl = -max(Re(eig(A_c))) */
    double ev_re[BLC_MAX_STATES], ev_im[BLC_MAX_STATES];
    blc_eigenvalues(Ac, n, ev_re, ev_im);
    double max_re = -1e10;
    for (int i = 0; i < n; i++) {
        if (ev_re[i] > max_re) max_re = ev_re[i];
    }
    /** λ_cl is the convergence rate = -max(Re(eig)) */
    etf->lambda_cl = -max_re;
    if (etf->lambda_cl < 0.01) etf->lambda_cl = 0.01;

    /** Compute min inter-event time bound */
    etf->min_inter_event = (1.0 / etf->lambda_cl)
                            * log(1.0 + sigma * etf->beta);

    return 0;
}

bool blc_etf_check(BLCEventTriggeredFeedback* etf, const double* x,
                    int n, double time) {
    (void)time;
    if (!etf || !x) return true;

    /** Compute error: e = x_last - x */
    double e_norm = 0.0, x_norm = 0.0;
    for (int i = 0; i < n; i++) {
        double ei = etf->x_last[i] - x[i];
        etf->e[i] = ei;
        e_norm += ei * ei;
        x_norm += x[i] * x[i];
    }
    e_norm = sqrt(e_norm);
    x_norm = sqrt(x_norm);

    /** Trigger: ||e|| > σ ||x|| */
    return (e_norm > etf->sigma * x_norm);
}

void blc_etf_triggered(BLCEventTriggeredFeedback* etf,
                        const double* x, double time) {
    if (!etf || !x) return;

    for (int i = 0; i < BLC_MAX_STATES; i++) {
        etf->x_last[i] = x[i];
        etf->e[i] = 0.0;
    }
    etf->last_event_time = time;
    etf->event_count++;
}

double blc_etf_min_inter_event(const BLCEventTriggeredFeedback* etf) {
    return etf ? etf->min_inter_event : 0.0;
}

void blc_etf_stats(const BLCEventTriggeredFeedback* etf,
                    int* event_count, double* avg_interval) {
    if (!etf) return;
    if (event_count) *event_count = etf->event_count;
    /** Average interval = total time / event_count */
    if (avg_interval && etf->event_count > 0) {
        *avg_interval = etf->last_event_time / (double)etf->event_count;
    }
}

/* ================================================================
 * Self-Triggered Control
 *
 * At time t_k, the next transmission time is:
 *   t_{k+1} = t_k + τ(x(t_k))
 *
 * For quadratic Lyapunov function V(x) = x'Px:
 *   τ(x) = min{ t > 0 : V(e^{A_c t} x) ≤ ρ·V(x) }
 *
 * We approximate τ using a Taylor expansion of V(x(t)):
 *   V(x(t)) ≈ V(x₀) + V̇(x₀)·t + V̈(x₀)·t²/2
 * where V̇ ≤ -λ_min(Q)·||x₀||² (for stable A_c).
 *
 * Self-triggering eliminates continuous monitoring but is more
 * conservative (larger min inter-event time than event-triggering).
 * ================================================================ */

int blc_self_trig_init(BLCSelfTriggered* st, double rho,
                        double tau_min, double tau_max,
                        const double* Ac, const double* P, int n) {
    if (!st || rho <= 0.0 || rho >= 1.0 || !Ac || !P || n < 1) return -1;

    st->rho        = rho;
    st->tau_min    = tau_min;
    st->tau_max    = tau_max;
    st->next_event_time = 0.0;
    st->event_count = 0;
    st->avg_interval = 0.0;

    /** Store A_c and P */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            st->A_c[i][j] = Ac[i * n + j];
            st->P[i][j]   = P[i * n + j];
        }
    }

    return 0;
}

double blc_self_trig_next_event(BLCSelfTriggered* st, const double* x,
                                 int n, double current_time) {
    if (!st || !x) return current_time + st->tau_max;

    /** Compute V(x) = x'Px */
    double V = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V += x[i] * st->P[i][j] * x[j];
        }
    }

    /** Compute V̇ = x'(A_c'P + P A_c)x = -x'Qx
     *  For stable A_c, V̇ < 0, giving τ = (1-ρ)·V / |V̇| */
    double V_dot = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            /** A_c'P + P A_c */
            double Q_ij = 0.0;
            for (int k = 0; k < n; k++) {
                Q_ij += st->A_c[k][i] * st->P[k][j]  /** A_c'(i,k)P(k,j) */
                      + st->P[i][k] * st->A_c[k][j];  /** P(i,k)A_c(k,j) */
            }
            V_dot += x[i] * Q_ij * x[j];
        }
    }

    /** τ = (1-ρ) V / |V_dot| for linear bound */
    double tau;
    if (fabs(V_dot) > 1e-12) {
        tau = (1.0 - st->rho) * V / fabs(V_dot);
    } else {
        tau = st->tau_max;
    }

    /** Clamp to [tau_min, tau_max] */
    if (tau < st->tau_min) tau = st->tau_min;
    if (tau > st->tau_max) tau = st->tau_max;

    st->next_event_time = current_time + tau;
    st->event_count++;

    /** Update running average */
    double alpha = 0.1;
    st->avg_interval = (1.0 - alpha) * st->avg_interval + alpha * tau;

    return st->next_event_time;
}

void blc_self_trig_stats(const BLCSelfTriggered* st,
                          double* avg_interval, int* event_count) {
    if (!st) return;
    if (avg_interval) *avg_interval = st->avg_interval;
    if (event_count)  *event_count  = st->event_count;
}