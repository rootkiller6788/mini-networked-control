/* cloud_control_delay.c - Network Delay Models & Compensation for Cloud Control */
#include "cloud_control_delay.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <assert.h>

/* LCG for reproducible random delay generation */
static unsigned int dlcg_state = 42;
static void dlcg_seed(unsigned int s) { dlcg_state = s; }
static double dlcg_rand(void) {
    dlcg_state = (1103515245 * dlcg_state + 12345) & 0x7fffffff;
    return (double)dlcg_state / 0x7fffffff;
}
/* Box-Muller for normal distribution */
static double rand_normal(double mean, double std) {
    double u1 = dlcg_rand(), u2 = dlcg_rand();
    if (u1 < 1e-10) u1 = 1e-10;
    return mean + std * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* -------- Delay Statistics -------- */

int delay_stats_compute(const double *delays, int count,
                         DelayStatistics *stats_out) {
    if (!delays || count <= 0 || !stats_out) return -1;
    double sum = 0.0, sum2 = 0.0, min_v = 1e300, max_v = 0.0;
    for (int i = 0; i < count; i++) {
        double v = delays[i];
        sum += v; sum2 += v*v;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    double mean = sum / count;
    double var = sum2/count - mean*mean;
    if (var < 0) var = 0;
    stats_out->mean_us = mean;
    stats_out->variance_us = var;
    stats_out->min_us = min_v;
    stats_out->max_us = max_v;
    stats_out->sample_count = count;
    /* Sort a copy for percentiles */
    double *sorted = (double*)malloc((size_t)count * sizeof(double));
    if (!sorted) return -1;
    for (int i = 0; i < count; i++) sorted[i] = delays[i];
    for (int i = 0; i < count-1; i++)
        for (int j = i+1; j < count; j++)
            if (sorted[i] > sorted[j]) { double t = sorted[i]; sorted[i] = sorted[j]; sorted[j] = t; }
    stats_out->median_us = (count % 2) ? sorted[count/2] : (sorted[count/2-1] + sorted[count/2])/2.0;
    int p95_idx = (int)(count * 0.95); if (p95_idx >= count) p95_idx = count-1;
    int p99_idx = (int)(count * 0.99); if (p99_idx >= count) p99_idx = count-1;
    stats_out->p95_us = sorted[p95_idx];
    stats_out->p99_us = sorted[p99_idx];
    /* Jitter: mean absolute deviation */
    double jit_sum = 0.0;
    for (int i = 0; i < count; i++) jit_sum += fabs(delays[i] - mean);
    stats_out->jitter_us = jit_sum / count;
    /* Lag-1 autocorrelation */
    double num = 0.0, den = 0.0;
    for (int i = 0; i < count-1; i++) {
        num += (delays[i] - mean) * (delays[i+1] - mean);
        den += (delays[i] - mean) * (delays[i] - mean);
    }
    stats_out->autocorrelation = (den > 1e-12) ? num/den : 0.0;
    /* Outlier count */
    double std = sqrt(var);
    stats_out->outliers_count = 0;
    for (int i = 0; i < count; i++)
        if (fabs(delays[i] - mean) > 3.0 * std) stats_out->outliers_count++;
    free(sorted);
    return 0;
}

void delay_stats_print(const DelayStatistics *stats) {
    if (!stats) return;
    printf("Delay Statistics (n=%d):\n", stats->sample_count);
    printf("  Mean=%.1fus Median=%.1fus\n", stats->mean_us, stats->median_us);
    printf("  Min=%.1fus Max=%.1fus P95=%.1fus P99=%.1fus\n",
           stats->min_us, stats->max_us, stats->p95_us, stats->p99_us);
    printf("  Jitter=%.1fus Autocorr=%.3f Outliers=%d\n",
           stats->jitter_us, stats->autocorrelation, stats->outliers_count);
}

/* -------- Delay Stability Analysis -------- */

int delay_is_stable(const CloudControlSystem *ccs, double delay_us) {
    if (!ccs) return -1;
    int n = ccs->n, m = ccs->m;
    double *A_cl = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!A_cl) return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            A_cl[i*n+j] = ccs->A[i][j];
            for (int k = 0; k < m; k++)
                A_cl[i*n+j] -= ccs->B[i][k] * ccs->K[k][j];
        }
    /* Power iteration to estimate spectral radius */
    double *v = (double*)calloc((size_t)n, sizeof(double));
    double *Av = (double*)calloc((size_t)n, sizeof(double));
    if (!v || !Av) { free(v); free(Av); free(A_cl); return -1; }
    for (int i = 0; i < n; i++) v[i] = 1.0/sqrt((double)n);
    double lambda = 0.0;
    for (int iter = 0; iter < 100; iter++) {
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++) Av[i] += A_cl[i*n+j] * v[j];
        }
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += Av[i]*Av[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        double new_lambda = 0.0;
        for (int i = 0; i < n; i++) { v[i] = Av[i]/norm; new_lambda += v[i]*Av[i]; }
        if (fabs(new_lambda - lambda) < 1e-10) { lambda = new_lambda; break; }
        lambda = new_lambda;
    }
    free(v); free(Av);
    /* Delay-dependent condition: tau * |lambda_max| < 1 (crude bound) */
    double tau_s = delay_us / 1e6;
    int stable = (fabs(lambda) * tau_s < 0.5) ? 1 : 0;
    free(A_cl);
    return stable;
}

double delay_mati_compute(const CloudControlSystem *ccs) {
    if (!ccs || ccs->n <= 0) return -1.0;
    int n = ccs->n, m = ccs->m;
    double *A_cl = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!A_cl) return -1.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            A_cl[i*n+j] = ccs->A[i][j];
            for (int k = 0; k < m; k++)
                A_cl[i*n+j] -= ccs->B[i][k] * ccs->K[k][j];
        }
    /* 2-norm of A_cl */
    double norm = 0.0;
    for (int i = 0; i < n*n; i++) norm += A_cl[i]*A_cl[i];
    norm = sqrt(norm);
    free(A_cl);
    if (norm < 1e-12) return CCS_MAX_DELAY_US;
    double tau_max = 1.0 / (2.0 * norm);
    return tau_max * 1e6;
}

int delay_lyapunov_krasovskii(const CloudControlSystem *ccs, double delay_us) {
    if (!ccs) return -1;
    /* Simplified LK check: verify P = I yields negative derivative for small tau */
    int n = ccs->n, m = ccs->m;
    double tau = delay_us / 1e6;
    double *A_cl = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!A_cl) return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            A_cl[i*n+j] = ccs->A[i][j];
            for (int k = 0; k < m; k++)
                A_cl[i*n+j] -= ccs->B[i][k] * ccs->K[k][j];
        }
    /* Check if A_cl + A_cl^T is negative definite */
    double max_eig = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++)
            row_sum += fabs(A_cl[i*n+j] + A_cl[j*n+i]);
        if (row_sum > max_eig) max_eig = row_sum;
    }
    free(A_cl);
    return (max_eig * tau < 1.0) ? 1 : 0;
}
/* -------- Smith Predictor Implementation -------- */

SmithPredictor* smith_create(CloudControlSystem *ccs) {
    if (!ccs) return NULL;
    SmithPredictor *sp = (SmithPredictor*)calloc(1, sizeof(SmithPredictor));
    if (!sp) return NULL;
    sp->ccs = ccs;
    sp->n_m = ccs->n; sp->m_m = ccs->m; sp->p_m = ccs->p;
    /* Copy plant model as internal model */
    for (int i = 0; i < sp->n_m; i++) {
        for (int j = 0; j < sp->n_m; j++) sp->A_m[i][j] = ccs->A[i][j];
        sp->x_model[i] = ccs->plant_state.x[i];
        sp->x_delayed[i] = ccs->plant_state.x[i];
    }
    for (int i = 0; i < sp->n_m; i++)
        for (int j = 0; j < sp->m_m; j++) sp->B_m[i][j] = ccs->B[i][j];
    for (int i = 0; i < sp->p_m; i++)
        for (int j = 0; j < sp->n_m; j++) sp->C_m[i][j] = ccs->C[i][j];
    sp->estimated_delay_us = CCS_DEFAULT_DELAY_US;
    sp->delay_adaptation_rate = 0.01;
    return sp;
}

void smith_free(SmithPredictor *sp) { free(sp); }

void smith_set_delay(SmithPredictor *sp, double delay_us) {
    if (!sp) return;
    sp->estimated_delay_us = delay_us;
}

int smith_predict(SmithPredictor *sp, const double *y_meas,
                   const double *u, double dt, double *y_comp) {
    if (!sp || !y_meas || !u || !y_comp) return -1;
    int n = sp->n_m, m = sp->m_m, p = sp->p_m;
    /* Step 1: advance delay-free model x_m_dot = A_m*x_m + B_m*u */
    double dx_m[CCS_MAX_STATES] = {0};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) dx_m[i] += sp->A_m[i][j] * sp->x_model[j];
        for (int j = 0; j < m; j++) dx_m[i] += sp->B_m[i][j] * u[j];
        sp->x_model[i] += dt * dx_m[i];
    }
    /* Step 2: history-based delayed model (use stored old input) */
    double u_old[CCS_MAX_INPUTS] = {0};
    if (sp->u_history_count > 0) {
        int idx = sp->u_history_count - 1;
        for (int j = 0; j < m && j < CCS_MAX_INPUTS; j++)
            u_old[j] = sp->u_history[idx][j];
    }
    double dx_d[CCS_MAX_STATES] = {0};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) dx_d[i] += sp->A_m[i][j] * sp->x_delayed[j];
        for (int j = 0; j < m; j++) dx_d[i] += sp->B_m[i][j] * u_old[j];
        sp->x_delayed[i] += dt * dx_d[i];
    }
    /* Step 3: compute delay-free model output */
    double y_model[CCS_MAX_OUTPUTS] = {0};
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            y_model[i] += sp->C_m[i][j] * sp->x_model[j];
    /* Step 4: compute delayed model output */
    double y_delayed[CCS_MAX_OUTPUTS] = {0};
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            y_delayed[i] += sp->C_m[i][j] * sp->x_delayed[j];
    /* Step 5: Smith compensation y_comp = y_meas[i] - y_delayed[i] + y_model[i] */
    for (int i = 0; i < p; i++)
        y_comp[i] = y_meas[i] - y_delayed[i] + y_model[i];
    /* Step 6: prediction error tracking */
    double err_sum = 0.0;
    for (int i = 0; i < p; i++) err_sum += (y_meas[i] - y_delayed[i]) * (y_meas[i] - y_delayed[i]);
    sp->prediction_error = sqrt(err_sum);
    /* Store u in history */
    if (sp->u_history_count < CCS_MAX_DELAY_HISTORY) {
        for (int j = 0; j < m && j < CCS_MAX_INPUTS; j++)
            sp->u_history[sp->u_history_count][j] = u[j];
        sp->u_history_count++;
    }
    return 0;
}

double smith_adapt_delay(SmithPredictor *sp) {
    if (!sp) return -1.0;
    double grad = sp->prediction_error * sp->delay_adaptation_rate;
    sp->estimated_delay_us += grad * 100.0;
    if (sp->estimated_delay_us < 0) sp->estimated_delay_us = 0;
    if (sp->estimated_delay_us > CCS_MAX_DELAY_US) sp->estimated_delay_us = CCS_MAX_DELAY_US;
    return sp->estimated_delay_us;
}

double smith_get_prediction_error(const SmithPredictor *sp) {
    return sp ? sp->prediction_error : -1.0;
}

int smith_handle_oos_measurement(SmithPredictor *sp,
                                  const double *y_meas, double timestamp) {
    if (!sp || !y_meas) return -1;
    /* Insert OOS measurement at correct position in measurement buffer */
    if (sp->meas_count >= CCS_MAX_HISTORY) {
        for (int i = 0; i < CCS_MAX_HISTORY - 1; i++) {
            sp->meas_ts[i] = sp->meas_ts[i+1];
            for (int j = 0; j < sp->p_m; j++)
                sp->meas_buffer[i][j] = sp->meas_buffer[i+1][j];
        }
        sp->meas_count = CCS_MAX_HISTORY - 1;
    }
    int pos = sp->meas_count;
    for (int i = 0; i < sp->meas_count; i++) {
        if (timestamp < sp->meas_ts[i]) { pos = i; break; }
    }
    /* Shift to make room */
    for (int i = sp->meas_count; i > pos; i--) {
        sp->meas_ts[i] = sp->meas_ts[i-1];
        for (int j = 0; j < sp->p_m; j++)
            sp->meas_buffer[i][j] = sp->meas_buffer[i-1][j];
    }
    sp->meas_ts[pos] = timestamp;
    for (int j = 0; j < sp->p_m; j++)
        sp->meas_buffer[pos][j] = y_meas[j];
    sp->meas_count++;
    return 0;
}
/* -------- Delay Generation -------- */

double delay_generate(DelayModelType model, double param1, double param2) {
    double d = 0.0;
    switch (model) {
    case DELAY_MODEL_CONSTANT:
        d = param1; break;
    case DELAY_MODEL_UNIFORM:
        d = param1 + dlcg_rand() * (param2 - param1); break;
    case DELAY_MODEL_NORMAL: {
        d = rand_normal(param1, param2);
        if (d < 0) d = 0;
        break;
    }
    case DELAY_MODEL_EXPONENTIAL:
        d = param1 - param1 * log(dlcg_rand() + 1e-10); break;
    case DELAY_MODEL_GAMMA: {
        int k = (int)param1; if (k < 1) k = 1;
        double theta = param2;
        d = 0.0;
        for (int i = 0; i < k; i++)
            d += -theta * log(dlcg_rand() + 1e-10);
        break;
    }
    case DELAY_MODEL_PERIODIC:
        d = param1 + param2 * sin(2.0 * M_PI * dlcg_rand()); break;
    case DELAY_MODEL_MARKOV: {
        static double state = 0.0;
        if (dlcg_rand() < 0.1) state = param1 + dlcg_rand() * (param2 - param1);
        d = state;
        break;
    }
    case DELAY_MODEL_MEASURED:
    default:
        d = param1; break;
    }
    return d > 0 ? d : 0;
}

void delay_trace_generate(DelayModelType model, double param1, double param2,
                           int count, double *out) {
    if (!out || count <= 0) return;
    dlcg_seed((unsigned int)time(NULL));
    for (int i = 0; i < count; i++)
        out[i] = delay_generate(model, param1, param2);
}

double delay_estimate_online(double *current_ema, double new_measurement,
                              double alpha) {
    if (!current_ema) return new_measurement;
    *current_ema = alpha * new_measurement + (1.0 - alpha) * (*current_ema);
    return *current_ema;
}

double delay_jitter_compute(const double *delays, int count) {
    if (!delays || count < 2) return 0.0;
    double sum = 0.0;
    for (int i = 1; i < count; i++)
        sum += fabs(delays[i] - delays[i-1]);
    return sum / (count - 1);
}

double delay_packet_loss_correlate(const double *delays, const int *lost,
                                     int count) {
    if (!delays || !lost || count < 2) return 0.0;
    double mean_d = 0.0, mean_l = 0.0;
    for (int i = 0; i < count; i++) { mean_d += delays[i]; mean_l += lost[i]; }
    mean_d /= count; mean_l /= count;
    double num = 0.0, den_d = 0.0, den_l = 0.0;
    for (int i = 0; i < count; i++) {
        double dd = delays[i] - mean_d, dl = lost[i] - mean_l;
        num += dd * dl; den_d += dd*dd; den_l += dl*dl;
    }
    double den = sqrt(den_d * den_l);
    return (den > 1e-12) ? num/den : 0.0;
}

/* -------- Combined Delay-Compensated Step -------- */

int delay_compensated_step(CloudControlSystem *ccs, SmithPredictor *sp,
                            const double *measurement, double ts,
                            double dt, CompensationStrategy strategy) {
    if (!ccs) return -1;
    (void)ccs->p;
    double y_comp[CCS_MAX_OUTPUTS] = {0};
    double u[CCS_MAX_INPUTS] = {0};
    /* Generate simulated network delay for this cycle */
    double sim_delay = delay_generate(DELAY_MODEL_EXPONENTIAL,
                                       CCS_DEFAULT_DELAY_US, 0);
    if (measurement) {
        if (sp && strategy == COMPENSATE_SMITH) {
            /* Use Smith predictor to compensate delay */
            ccs_compute_control(ccs, u);
            smith_predict(sp, measurement, u, dt, y_comp);
            ccs_update_observer(ccs, y_comp, ts);
        } else if (strategy == COMPENSATE_NONE) {
            ccs_update_observer(ccs, measurement, ts);
        }
    }
    ccs_compute_control(ccs, u);
    ccs_apply_control(ccs, u, dt);
    ccs_step(ccs, measurement, ts, dt, sim_delay);
    return 0;
}
/* ============================================================================
 * Extended Delay Analysis Functions
 * ============================================================================ */

/* -------- Delay Distribution Fitting -------- */

int delay_fit_exponential(const double *delays, int count,
                           double *rate_out) {
    if (!delays || count <= 0 || !rate_out) return -1;
    double sum = 0.0;
    for (int i = 0; i < count; i++) sum += delays[i];
    if (sum < 1e-12) return -1;
    *rate_out = count / sum;
    return 0;
}

int delay_fit_normal(const double *delays, int count,
                      double *mean_out, double *std_out) {
    if (!delays || count < 2 || !mean_out || !std_out) return -1;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < count; i++) { sum += delays[i]; sum2 += delays[i]*delays[i]; }
    *mean_out = sum / count;
    double var = sum2/count - (*mean_out)*(*mean_out);
    *std_out = var > 0 ? sqrt(var) : 0.0;
    return 0;
}

/* -------- Delay Bound Prediction -------- */

double delay_predict_bound(const DelayStatistics *stats, double confidence) {
    if (!stats || confidence <= 0 || confidence >= 1.0) return -1.0;
    /* Chebyshev bound: P(|X-mu| >= k*sigma) <= 1/k^2 */
    double k = sqrt(1.0 / (1.0 - confidence));
    double std = sqrt(stats->variance_us);
    return stats->mean_us + k * std;
}

double delay_quantile_estimate(const double *sorted_delays, int count,
                                double quantile) {
    if (!sorted_delays || count <= 0 || quantile < 0 || quantile > 1) return -1.0;
    int idx = (int)(count * quantile);
    if (idx >= count) idx = count - 1;
    return sorted_delays[idx];
}

/* -------- Control-aware Delay Analysis -------- */

double delay_control_cost(const CloudControlSystem *ccs,
                           const double *delays, int count) {
    if (!ccs || !delays || count <= 0) return -1.0;
    /* Estimate control cost increase due to delay */
    /* Based on linear quadratic Gaussian (LQG) sensitivity to delay */
    double mean_delay = 0.0;
    for (int i = 0; i < count; i++) mean_delay += delays[i];
    mean_delay /= count;
    double tau_s = mean_delay / 1e6;
    /* Simplified cost model: J(tau) = J0 * exp(alpha * tau) */
    double J0 = ccs->ise > 0 ? ccs->ise : 1.0;
    double alpha = 0.1; /* Cost sensitivity parameter */
    return J0 * exp(alpha * tau_s);
}

/* -------- Adaptive Threshold Computation -------- */

double delay_adaptive_threshold(const DelayStatistics *stats,
                                 double safety_factor) {
    if (!stats || safety_factor <= 0) return -1.0;
    /* Adaptive threshold = mean + safety_factor * stddev */
    double std = sqrt(stats->variance_us);
    return stats->mean_us + safety_factor * std;
}

int delay_is_anomalous(const DelayStatistics *baseline,
                        double new_delay, double sigma_threshold) {
    if (!baseline) return 0;
    double std = sqrt(baseline->variance_us);
    if (std < 1e-9) return 0;
    return (fabs(new_delay - baseline->mean_us) > sigma_threshold * std) ? 1 : 0;
}

/* -------- Network Quality Scoring -------- */

double delay_network_quality_score(const DelayStatistics *stats,
                                    double target_latency_us) {
    if (!stats || target_latency_us <= 0) return 0.0;
    /* Score from 0 (bad) to 1 (excellent) based on how much delay margin exists */
    double margin = target_latency_us - stats->p95_us;
    if (margin <= 0) return 0.0;
    double score = margin / target_latency_us;
    /* Penalize jitter and outliers */
    double jitter_penalty = stats->jitter_us / (stats->mean_us + 1e-12);
    double outlier_penalty = (double)stats->outliers_count / (stats->sample_count + 1);
    score *= (1.0 - 0.3 * jitter_penalty - 0.2 * outlier_penalty);
    return score > 0 ? score : 0.0;
}

/* -------- Batch Delay Processing -------- */

int delay_process_batch(const double *delays, int count,
                         DelayStatistics *stats, double *fitted_rate) {
    if (!delays || count <= 0) return -1;
    if (stats) delay_stats_compute(delays, count, stats);
    if (fitted_rate) delay_fit_exponential(delays, count, fitted_rate);
    return 0;
}
