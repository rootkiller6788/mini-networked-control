/**
 * @file    qc_data_rate.c
 * @brief   Data Rate Theorem implementation and channel capacity analysis
 *
 * Implements the Data Rate Theorem (Nair & Evans, 2004):
 * For discrete-time LTI system x_{k+1} = A x_k + B u_k with quantized
 * state feedback at rate R bits/sample, asymptotic stabilizability
 * requires R > sum_i max(0, log2|lambda_i(A)|).
 *
 * Also implements channel capacity computations (AWGN, BSC, erasure),
 * rate-distortion theory (Gaussian source), reverse water-filling
 * for optimal rate allocation, and anytime capacity analysis.
 *
 * Key references:
 *   - Nair & Evans (2004). SIAM J. Control Optim. 42(6): 1953-1969.
 *   - Tatikonda & Mitter (2004). IEEE TAC 49(7): 1056-1068.
 *   - Cover & Thomas (2006). Elements of Information Theory. 2nd ed.
 *   - Sahai & Mitter (2006). IEEE TIT 52(8): 3664-3686.
 *
 * Applications (L7):
 *   - Smart grid: phasor measurement units (PMU) with GPS-synchronized sampling
 *   - Aerospace: Boeing 787 / Airbus flight control with ARINC 429 data buses
 *   - Nuclear: Fukushima Daiichi lessons — control under severe data loss
 *   - Climate monitoring: satellite telemetry with constrained downlink rates
 */

#include "quantized_control.h"
#include "qc_data_rate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

/* ================================================================
 * Data Rate Initialization and Configuration
 * ================================================================ */

void qc_data_rate_init(QCDataRate *dr, double total_rate, int num_ch) {
    if (!dr) return;
    memset(dr, 0, sizeof(QCDataRate));
    dr->total_rate_bps = total_rate;
    dr->num_channels = (num_ch > 0) ? num_ch : 1;
    dr->channel_bits = calloc(dr->num_channels, sizeof(int));
    dr->channel_rate = calloc(dr->num_channels, sizeof(double));
    dr->eig_magnitudes = calloc(dr->num_channels, sizeof(double));
    dr->min_rate = 0.0;
    dr->rate_margin = total_rate;
    dr->channel_capacity = total_rate;
}

int qc_data_rate_compute_min(QCDataRate *dr, const double *A, int n) {
    if (!dr || !A || n <= 0) return -1;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1; }
    qc_matrix_eigenvalues(A, n, er, ei);

    dr->num_unstable = 0;
    dr->min_rate = 0.0;
    for (int i = 0; i < n; i++) {
        double mag = sqrt(er[i]*er[i] + ei[i]*ei[i]);
        if (mag > 1.0) {
            dr->min_rate += log2(mag);
            dr->num_unstable++;
            if (dr->num_unstable <= dr->num_channels) {
                dr->eig_magnitudes[dr->num_unstable - 1] = mag;
            }
        }
    }

    dr->rate_margin = dr->total_rate_bps - dr->min_rate;
    free(er); free(ei);
    return 0;
}

int qc_data_rate_allocate(QCDataRate *dr) {
    if (!dr || dr->num_unstable <= 0) return -1;
    /* Allocate bits proportionally to log2|lambda| */
    double total_weight = 0.0;
    for (int i = 0; i < dr->num_unstable && i < dr->num_channels; i++) {
        total_weight += log2(dr->eig_magnitudes[i]);
    }
    if (total_weight <= 0.0) return -1;
    for (int i = 0; i < dr->num_unstable && i < dr->num_channels; i++) {
        double weight = log2(dr->eig_magnitudes[i]);
        dr->channel_rate[i] = dr->total_rate_bps * weight / total_weight;
        dr->channel_bits[i] = (int)ceil(dr->channel_rate[i]);
    }
    return 0;
}

int qc_data_rate_check_stabilizability(const QCDataRate *dr) {
    if (!dr) return 0;
    return (dr->total_rate_bps > dr->min_rate) ? 1 : 0;
}

/* ================================================================
 * Data Rate Theorem Analysis
 * ================================================================ */

int qc_data_rate_analyze(const double *A, int n, double dt,
                          QCDataRateResult *result) {
    if (!A || !result || n <= 0) return -1;
    memset(result, 0, sizeof(QCDataRateResult));

    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1; }
    qc_matrix_eigenvalues(A, n, er, ei);

    result->num_eigenvalues = n;
    result->eigenvalues_mag = calloc(n, sizeof(double));
    if (!result->eigenvalues_mag) {
        free(er); free(ei); return -1;
    }

    double min_rate = 0.0;
    int num_unstable = 0;
    result->instability_sum = 0.0;

    for (int i = 0; i < n; i++) {
        double mag = sqrt(er[i]*er[i] + ei[i]*ei[i]);
        result->eigenvalues_mag[i] = mag;
        if (mag > 1.0) {
            min_rate += log2(mag);
            num_unstable++;
            result->instability_sum += mag - 1.0;
        }
    }

    result->num_unstable = num_unstable;
    result->min_rate = min_rate;
    result->discretized = (dt > 0) ? 1 : 0;

    if (dt > 0) {
        /* Continuous-time case: eigenvalues of exp(A*dt)
         * log|lambda_dt| = log|exp(lambda_ct * dt)| = Re(lambda_ct) * dt * log2(e) */
        result->min_rate_ct = min_rate / dt;
    } else {
        result->min_rate_ct = min_rate;
    }

    free(er); free(ei);
    return 0;
}

double qc_data_rate_ct_minimum(const double *A, int n) {
    if (!A || n <= 0) return -1.0;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1.0; }
    qc_matrix_eigenvalues(A, n, er, ei);

    double rate = 0.0;
    for (int i = 0; i < n; i++) {
        /* For continuous-time, rate = sum max(0, Re(lambda)/ln(2)) */
        if (er[i] > 0) rate += er[i] / log(2.0);
    }
    free(er); free(ei);
    return rate;
}

int qc_data_rate_is_stabilizable(const double *A, int n, double rate,
                                  double dt) {
    if (!A || n <= 0 || rate <= 0.0) return 0;
    QCDataRateResult result;
    if (qc_data_rate_analyze(A, n, dt, &result) != 0) return 0;
    int stab = (rate > result.min_rate) ? 1 : 0;
    result.actual_rate = rate;
    result.stabilizable = stab;
    result.rate_surplus = rate - result.min_rate;
    qc_data_rate_result_free(&result);
    return stab;
}

int qc_data_rate_per_mode(const double *A, int n, double *mode_rates) {
    if (!A || !mode_rates || n <= 0) return -1;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1; }
    qc_matrix_eigenvalues(A, n, er, ei);

    for (int i = 0; i < n; i++) {
        double mag = sqrt(er[i]*er[i] + ei[i]*ei[i]);
        mode_rates[i] = (mag > 1.0) ? log2(mag) : 0.0;
    }
    free(er); free(ei);
    return 0;
}

double qc_data_rate_theoretical_min(const double *A, int n) {
    if (!A || n <= 0) return 0.0;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return 0.0; }
    qc_matrix_eigenvalues(A, n, er, ei);

    double rate = 0.0;
    for (int i = 0; i < n; i++) {
        double mag = sqrt(er[i]*er[i] + ei[i]*ei[i]);
        if (mag > 1.0) rate += log2(mag);
    }
    free(er); free(ei);
    return rate;
}

/* ================================================================
 * Channel Capacity Computations
 * ================================================================ */

void qc_channel_init(QCChannel *ch, QCChannelType type, double capacity,
                      double ber, double snr_db) {
    if (!ch) return;
    memset(ch, 0, sizeof(QCChannel));
    ch->type = type;
    ch->capacity = capacity;
    ch->bit_error_prob = ber;
    ch->snr = snr_db;
}

double qc_channel_capacity_awgn(double bandwidth, double snr_linear) {
    if (bandwidth <= 0.0 || snr_linear < 0.0) return 0.0;
    return bandwidth * log2(1.0 + snr_linear);
}

double qc_channel_capacity_bsc(double crossover_prob) {
    if (crossover_prob < 0.0 || crossover_prob > 1.0) return 0.0;
    if (crossover_prob == 0.0 || crossover_prob == 1.0) return 1.0;
    return 1.0 - qc_binary_entropy(crossover_prob);
}

double qc_channel_capacity_erasure(double erasure_prob) {
    if (erasure_prob < 0.0 || erasure_prob > 1.0) return 0.0;
    return 1.0 - erasure_prob;
}

double qc_binary_entropy(double p) {
    if (p <= 0.0 || p >= 1.0) return 0.0;
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p);
}

/* ================================================================
 * Rate-Distortion Theory
 * ================================================================ */

double qc_rate_distortion_gaussian_dr(double variance, double distortion) {
    /* R(D) = 0.5 * log2(variance / D) for 0 < D <= variance */
    if (variance <= 0.0 || distortion <= 0.0 || distortion >= variance) return 0.0;
    return 0.5 * log2(variance / distortion);
}

double qc_rate_distortion_gaussian(double variance, double distortion) {
    return qc_rate_distortion_gaussian_dr(variance, distortion);
}

int qc_bits_for_snr(double target_snr_db, double signal_std,
                     double range_half) {
    if (target_snr_db <= 0.0 || signal_std <= 0.0 || range_half <= 0.0) return 1;
    double target_snr_linear = pow(10.0, target_snr_db / 10.0);
    double noise_var = signal_std * signal_std / target_snr_linear;
    int bits = 1;
    for (; bits <= QC_MAX_BITS; bits++) {
        double step = 2.0 * range_half / (double)((1 << bits) - 1);
        double q_noise = step * step / 12.0;
        if (q_noise <= noise_var) break;
    }
    return bits;
}

int qc_reverse_waterfill(const double *noise_var, int num_ch,
                          double total_rate, double *rates) {
    if (!noise_var || !rates || num_ch <= 0 || total_rate <= 0.0) return -1;

    double *sorted_var = calloc(num_ch, sizeof(double));
    int *idx = calloc(num_ch, sizeof(int));
    if (!sorted_var || !idx) { free(sorted_var); free(idx); return -1; }
    for (int i = 0; i < num_ch; i++) { sorted_var[i] = noise_var[i]; idx[i] = i; }

    /* Sort by noise variance (descending) */
    for (int i = 0; i < num_ch; i++) {
        for (int j = i + 1; j < num_ch; j++) {
            if (sorted_var[j] > sorted_var[i]) {
                double tv = sorted_var[i]; sorted_var[i] = sorted_var[j]; sorted_var[j] = tv;
                int ti = idx[i]; idx[i] = idx[j]; idx[j] = ti;
            }
        }
    }

    /* Find water level w via bisection */
    double w_lo = 0.0, w_hi = sorted_var[0];
    for (int iter = 0; iter < 100; iter++) {
        double w = (w_lo + w_hi) / 2.0;
        double sum_rate = 0.0;
        for (int i = 0; i < num_ch; i++) {
            if (sorted_var[i] > w) {
                sum_rate += 0.5 * log2(sorted_var[i] / w);
            }
        }
        if (sum_rate > total_rate) w_lo = w;
        else w_hi = w;
    }

    /* Allocate rates */
    double w = (w_lo + w_hi) / 2.0;
    memset(rates, 0, num_ch * sizeof(double));
    for (int i = 0; i < num_ch; i++) {
        if (sorted_var[i] > w) {
            rates[idx[i]] = 0.5 * log2(sorted_var[i] / w);
        }
    }

    free(sorted_var); free(idx);
    return 0;
}

double qc_data_rate_entropy_power(const double *signal, int len, double var) {
    if (!signal || len <= 0) return 0.0;
    /* Entropy power of Gaussian: N(X) = (1 / (2*pi*e)) * exp(2*h(X))
     * For Gaussian: N(X) = var */
    double sum_sq = 0.0;
    for (int i = 0; i < len; i++) sum_sq += signal[i] * signal[i];
    double sq = sum_sq / len;
    /* If signal were optimally coded, its entropy power bounds distortion */
    (void)var;
    return sq * exp(2.0 * (-0.5 * log2(2.0 * M_PI * M_E * sq)));
}

int qc_data_rate_optimal_allocation(QCDataRate *dr, const double *ch_vars,
                                     int num_ch) {
    if (!dr || !ch_vars || num_ch <= 0) return -1;
    return qc_reverse_waterfill(ch_vars, num_ch, dr->total_rate_bps,
                                 dr->channel_rate);
}

double qc_anytime_capacity(double shannon_capacity, double delay,
                            double reliability_exponent) {
    /* Anytime capacity relates to error exponent:
     * For AWGN: C_anytime = C - (alpha / delay)
     * where alpha is the reliability exponent */
    if (delay <= 0.0) return 0.0;
    return shannon_capacity - reliability_exponent / delay;
}

int qc_data_rate_channel_feasible(double capacity, double required_rate) {
    return (capacity > required_rate + QC_EPSILON) ? 1 : 0;
}

void qc_data_rate_result_free(QCDataRateResult *r) {
    if (!r) return;
    free(r->eigenvalues_mag);
    r->eigenvalues_mag = NULL;
}

void qc_time_varying_result_free(QCTimeVaryingRateResult *r) {
    if (!r) return;
    free(r->window_min_rates);
    free(r->window_actual_rates);
    free(r->window_time);
    free(r->window_stable);
    r->window_min_rates = NULL;
    r->window_actual_rates = NULL;
    r->window_time = NULL;
    r->window_stable = NULL;
}

/* ================================================================
 * Sector Bound Analysis
 * ================================================================ */

void qc_sector_bound_init(QCSectorBoundResult *sr) {
    if (!sr) return;
    memset(sr, 0, sizeof(QCSectorBoundResult));
    sr->nominally_stable = 0;
    sr->robustly_stable = 0;
}

int qc_sector_bound_analyze(const QCSystem *sys, QCSectorBoundResult *sr) {
    if (!sys || !sr) return -1;
    sr->delta = sys->sector_delta;
    sr->is_sector_bounded = (sys->sector_delta > 0.0) ? 1 : 0;
    sr->nominally_stable = qc_is_hurwitz_stable(sys->A, sys->state_dim);

    /* Circle criterion: if nominal system is stable and
     * ||G(jw)||_inf * delta < 1, then robustly stable under
     * sector [delta, 1] quantization nonlinearity. */
    if (sr->nominally_stable && sr->delta > 0.0) {
        double rho = qc_spectral_radius(sys->A, sys->state_dim);
        sr->L2_gain_bound = rho / (1.0 - rho + QC_EPSILON);
        sr->circle_criterion_margin = 1.0 / sr->L2_gain_bound - sr->delta;
        sr->robustly_stable = (sr->circle_criterion_margin > 0.0) ? 1 : 0;
    }

    return 0;
}

double qc_sector_bound_delta_for_rho(double rho) {
    return qc_log_sector_delta(rho);
}

int qc_sector_bound_feasibility(const QCSystem *sys, double *feasibility_margin) {
    if (!sys) return 0;
    double delta = sys->sector_delta;
    double rho = qc_spectral_radius(sys->A, sys->state_dim);
    if (feasibility_margin) {
        *feasibility_margin = (1.0 - rho) - delta * rho;
    }
    return ((1.0 - rho) > delta * rho) ? 1 : 0;
}

int qc_circle_criterion_check(const QCSystem *sys, double delta,
                               double *margin) {
    if (!sys || delta < 0.0) return 0;
    /* Small-gain theorem via circle criterion:
     * Stability if |G(e^{jw})| * delta < 1 for all w.
     * Simplified: check sup norm bound via spectral radius. */
    double rho = qc_spectral_radius(sys->A, sys->state_dim);
    double l2_gain = rho / (1.0 - rho + QC_EPSILON);
    if (margin) *margin = 1.0 / l2_gain - delta;
    return (l2_gain * delta < 1.0) ? 1 : 0;
}
