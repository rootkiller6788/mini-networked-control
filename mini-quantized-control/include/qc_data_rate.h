/**
 * @file    qc_data_rate.h
 * @brief   Data Rate Theorem structures and analysis API
 *
 * Implements the Data Rate Theorem (Nair & Evans, 2004) characterizing
 * the minimum information rate required to stabilize an unstable linear
 * system over a digital communication channel.
 *
 * Theorem: For discrete-time LTI system x_{k+1} = A x_k + B u_k
 * with quantized state feedback at rate R bits/sample, asymptotic
 * stabilizability requires R > sum_i max(0, log2|lambda_i(A)|).
 *
 * References:
 *   - Nair & Evans (2004). SIAM J. Control Optim.
 *   - Tatikonda & Mitter (2004). IEEE TAC.
 *   - Wong & Brockett (1999). IEEE TAC.
 */

#ifndef QC_DATA_RATE_H
#define QC_DATA_RATE_H

#include "quantized_control.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    QC_CHANNEL_NOISELESS      = 0,
    QC_CHANNEL_BSC            = 1,
    QC_CHANNEL_AWGN           = 2,
    QC_CHANNEL_PACKET_ERASURE = 3,
    QC_CHANNEL_RATE_LIMITED   = 4
} QCChannelType;

typedef struct {
    QCChannelType   type;
    double          capacity;
    double          bit_error_prob;
    double          snr;
    double          erasure_prob;
    double          bandwidth;
    double          latency;
    double          jitter;
} QCChannel;

typedef struct {
    double          rate;
    double          distortion;
    double          rate_distortion_bound;
    int             achievable;
} QCRateDistortionPoint;

typedef struct {
    double         *eigenvalues_mag;
    int             num_unstable;
    int             num_eigenvalues;
    double          min_rate;
    double          actual_rate;
    int             stabilizable;
    double          rate_surplus;
    double          min_rate_ct;
    int             discretized;
    double          instability_sum;
} QCDataRateResult;

typedef struct {
    int             num_windows;
    double         *window_min_rates;
    double         *window_actual_rates;
    double         *window_time;
    int            *window_stable;
    double          time_average_rate;
    int             stable_overall;
} QCTimeVaryingRateResult;

/* Data Rate Theorem functions */
int qc_data_rate_analyze(const double *A, int n, double dt, QCDataRateResult *result);
double qc_data_rate_ct_minimum(const double *A, int n);
int qc_data_rate_is_stabilizable(const double *A, int n, double rate, double dt);
int qc_data_rate_per_mode(const double *A, int n, double *mode_rates);

/* Channel capacity functions */
void qc_channel_init(QCChannel *ch, QCChannelType type, double capacity, double ber, double snr_db);
double qc_channel_capacity_awgn(double bandwidth, double snr_linear);
double qc_channel_capacity_bsc(double crossover_prob);
double qc_channel_capacity_erasure(double erasure_prob);
double qc_binary_entropy(double p);

/* Rate-distortion functions */
double qc_rate_distortion_gaussian_dr(double variance, double distortion);
int qc_bits_for_snr(double target_snr_db, double signal_std, double range_half);
int qc_reverse_waterfill(const double *noise_var, int num_ch, double total_rate, double *rates);
double qc_anytime_capacity(double shannon_capacity, double delay, double reliability_exponent);

void qc_data_rate_result_free(QCDataRateResult *r);
void qc_time_varying_result_free(QCTimeVaryingRateResult *r);

#ifdef __cplusplus
}
#endif

#endif /* QC_DATA_RATE_H */
