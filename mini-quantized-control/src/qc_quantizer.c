/**
 * @file    qc_quantizer.c
 * @brief   Quantizer implementations: uniform, logarithmic, vector, dither
 *
 * Implements the fundamental quantizer operations including:
 *   - Uniform quantizer (mid-tread and mid-rise)
 *   - Logarithmic quantizer (Elia-Mitter formulation)
 *   - Vector quantizer with LBG training
 *   - Dither quantization (TPDF subtractive and non-subtractive)
 *   - Lloyd-Max optimal quantizer for Gaussian sources
 *   - Quantizer performance metrics
 *
 * Each quantizer maps continuous values to discrete codewords,
 * implementing the partition of R^n into quantization cells.
 *
 * Key references:
 *   - Widrow & Kollar (2008). Quantization Noise. Cambridge.
 *   - Gersho & Gray (1992). Vector Quantization and Signal Compression.
 *   - Elia & Mitter (2001). IEEE TAC.
 *   - Lloyd (1957). Least squares quantization in PCM. BSTJ.
 */

#include "quantized_control.h"
#include "qc_quantizer.h"
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
 * Core Quantizer Operations
 * ================================================================ */

void qc_quantizer_init(QCQuantizer *q, QCQuantizerType type, int bits) {
    if (!q) return;
    memset(q, 0, sizeof(QCQuantizer));
    q->type = type;
    q->bits = (bits > 0 && bits <= QC_MAX_BITS) ? bits : 8;
    q->num_levels = 1 << q->bits;
    q->range_min = -1.0;
    q->range_max = 1.0;
    q->step = (q->range_max - q->range_min) / (double)q->num_levels;
    q->rho = 0.5;
    q->delta = 0.01;
    q->overload = QC_OVERLOAD_SATURATE;
    q->zoom_factor = 2.0;
    q->is_logarithmic = (type == QC_QTYPE_LOGARITHMIC) ? 1 : 0;
    q->saturation_value = q->range_max;
}

int qc_quantizer_configure_range(QCQuantizer *q, double min_val, double max_val) {
    if (!q || min_val >= max_val) return -1;
    q->range_min = min_val;
    q->range_max = max_val;
    q->step = (max_val - min_val) / (double)q->num_levels;
    q->saturation_value = max_val;
    return 0;
}

int qc_quantizer_set_log_params(QCQuantizer *q, double rho, double delta) {
    if (!q || rho <= 0.0 || rho >= 1.0 || delta < 0.0) return -1;
    q->rho = rho;
    q->delta = delta;
    q->is_logarithmic = 1;
    return 0;
}

double qc_quantize_scalar(const QCQuantizer *q, double x) {
    if (!q) return 0.0;
    /* Clamp to quantizer range to handle saturation */
    double x_clamped = x;
    if (x_clamped > q->range_max) x_clamped = q->range_max;
    if (x_clamped < q->range_min) x_clamped = q->range_min;

    switch (q->type) {
        case QC_QTYPE_UNIFORM:
            return qc_uniform_quantize_midtread(x_clamped, q->step);
        case QC_QTYPE_LOGARITHMIC: {
            /* Logarithmic quantization:
             * q(x) = sign(x) * u_min * rho^{-round(log(x/u_min) / log(1/rho))}
             */
            if (fabs(x) < q->delta) return 0.0;
            double sign_x = (x > 0) ? 1.0 : -1.0;
            double abs_x = fabs(x);
            double log_val = log(abs_x / q->delta) / log(1.0 / q->rho);
            double idx = round(log_val);
            double level = q->delta * pow(1.0 / q->rho, idx);
            return sign_x * level;
        }
        case QC_QTYPE_RANDOM_DITHER:
            return qc_dither_quantize_nonsubtractive(x, q->step);
        default:
            return qc_uniform_quantize_midtread(x, q->step);
    }
}

void qc_quantize_vector(const QCQuantizer *q, const double *x, int dim, double *xq) {
    if (!q || !x || !xq || dim <= 0) return;
    for (int i = 0; i < dim; i++) {
        xq[i] = qc_quantize_scalar(q, x[i]);
    }
}

double qc_quantization_error(const QCQuantizer *q, double x) {
    if (!q) return 0.0;
    double xq = qc_quantize_scalar(q, x);
    return x - xq;
}

int qc_quantizer_num_levels(const QCQuantizer *q) {
    if (!q) return 0;
    return q->num_levels;
}

double qc_quantizer_snr(const QCQuantizer *q, double signal_power) {
    if (!q || signal_power <= 0.0) return 0.0;
    double noise_var = qc_quantization_noise_variance(q->step);
    if (noise_var <= 0.0) return INFINITY;
    return 10.0 * log10(signal_power / noise_var);
}

void qc_quantizer_free(QCQuantizer *q) {
    if (!q) return;
    free(q->levels);
    q->levels = NULL;
    q->levels_len = 0;
}

int qc_quantizer_find_cell(const QCQuantizer *q, double x,
                            QCQuantizationCell *cell) {
    if (!q || !cell) return -1;
    double xq = qc_quantize_scalar(q, x);
    cell->index = (int)((xq - q->range_min) / q->step);
    cell->dimension = 1;
    cell->center[0] = xq;
    cell->half_width[0] = q->step / 2.0;
    cell->code = cell->index;
    cell->volume = q->step;
    return 0;
}

double qc_quantizer_max_error(const QCQuantizer *q) {
    if (!q) return INFINITY;
    switch (q->type) {
        case QC_QTYPE_UNIFORM: return q->step / 2.0;
        case QC_QTYPE_LOGARITHMIC: return q->delta * q->range_max;
        default: return q->step / 2.0;
    }
}

/* ================================================================
 * Uniform Quantizer
 * ================================================================ */

double qc_uniform_quantize_midtread(double x, double step) {
    if (step <= 0.0) return x;
    return step * round(x / step);
}

double qc_uniform_quantize_midrise(double x, double step) {
    if (step <= 0.0) return x;
    return step * floor(x / step) + step / 2.0;
}

double qc_uniform_optimal_step(double range, int bits) {
    int levels = 1 << bits;
    if (levels <= 0) return range;
    return range / (double)(levels - 1);
}

double qc_uniform_overload_prob(double range, double sigma) {
    /* P(|X| > range/2) = 2 * Q(range/(2*sigma)) */
    double z = range / (2.0 * sigma);
    return erfc(z / sqrt(2.0));
}

/* ================================================================
 * Logarithmic Quantizer (Elia-Mitter formulation)
 * ================================================================ */

void qc_log_quantizer_init(QCLogQuantizer *lq, double rho, double u_min) {
    if (!lq) return;
    memset(lq, 0, sizeof(QCLogQuantizer));
    lq->rho = (rho > 0.0 && rho < 1.0) ? rho : 0.5;
    lq->u_min = (u_min > 0.0) ? u_min : 0.001;
    lq->deadzone = u_min;
    lq->sector_delta = (1.0 - lq->rho) / (1.0 + lq->rho);
    lq->num_positive = 32;
    lq->num_negative = 32;
    lq->total_levels = 65;
    lq->levels_positive = malloc(lq->num_positive * sizeof(double));
    lq->levels_negative = malloc(lq->num_negative * sizeof(double));
    if (lq->levels_positive) {
        for (int i = 0; i < lq->num_positive; i++) {
            lq->levels_positive[i] = u_min / pow(rho, (double)i);
        }
    }
    if (lq->levels_negative) {
        for (int i = 0; i < lq->num_negative; i++) {
            lq->levels_negative[i] = -u_min / pow(rho, (double)i);
        }
    }
}

double qc_log_quantize(const QCLogQuantizer *lq, double x) {
    if (!lq) return 0.0;
    double abs_x = fabs(x);
    if (abs_x < lq->deadzone) return 0.0;

    /* Find the closest positive level */
    int best_idx = 0;
    double best_dist = INFINITY;
    for (int i = 0; i < lq->num_positive; i++) {
        double dist = fabs(abs_x - lq->levels_positive[i]);
        if (dist < best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }
    if (best_idx < lq->num_positive) {
        double level = lq->levels_positive[best_idx];
        return (x > 0) ? level : -level;
    }
    return 0.0;
}

double qc_log_sector_delta(double rho) {
    if (rho <= 0.0 || rho >= 1.0) return 0.0;
    return (1.0 - rho) / (1.0 + rho);
}

double qc_log_quantizer_density(const QCLogQuantizer *lq) {
    if (!lq) return 0.0;
    return 2.0 * (1.0 - lq->rho) / (1.0 + lq->rho);
}

int qc_log_quantizer_find_level(const QCLogQuantizer *lq, double x,
                                 int *is_positive, int *level_index) {
    if (!lq || !is_positive || !level_index) return -1;
    double abs_x = fabs(x);
    if (abs_x < lq->deadzone) {
        *is_positive = 0;
        *level_index = -1;
        return 0;
    }
    *is_positive = (x > 0) ? 1 : 0;
    int best_idx = 0;
    double best_dist = INFINITY;
    for (int i = 0; i < lq->num_positive; i++) {
        double dist = fabs(abs_x - lq->levels_positive[i]);
        if (dist < best_dist) { best_dist = dist; best_idx = i; }
    }
    *level_index = best_idx;
    return 0;
}

void qc_log_quantizer_free(QCLogQuantizer *lq) {
    if (!lq) return;
    free(lq->levels_positive);
    free(lq->levels_negative);
    lq->levels_positive = NULL;
    lq->levels_negative = NULL;
}

int qc_log_quantizer_build_levels(QCLogQuantizer *lq, int num_levels) {
    if (!lq || num_levels < 3 || num_levels % 2 == 0) return -1;
    int half = (num_levels - 1) / 2;
    lq->num_positive = half;
    lq->num_negative = half;
    lq->total_levels = num_levels;
    free(lq->levels_positive);
    free(lq->levels_negative);
    lq->levels_positive = malloc(half * sizeof(double));
    lq->levels_negative = malloc(half * sizeof(double));
    if (!lq->levels_positive || !lq->levels_negative) return -1;
    for (int i = 0; i < half; i++) {
        lq->levels_positive[i] = lq->u_min / pow(lq->rho, (double)i);
        lq->levels_negative[i] = -lq->levels_positive[i];
    }
    return 0;
}

double qc_log_rho_from_delta(double delta) {
    if (delta <= 0.0 || delta >= 1.0) return 0.5;
    return (1.0 - delta) / (1.0 + delta);
}

int qc_log_num_levels_from_range(double u_min, double u_max, double rho) {
    if (u_min <= 0.0 || u_max <= u_min || rho <= 0.0 || rho >= 1.0) return 0;
    double n = log(u_max / u_min) / log(1.0 / rho);
    return (int)floor(n) + 1;
}

double qc_log_deadzone_from_umin(double u_min) {
    return (u_min > 0.0) ? u_min : 0.001;
}

double qc_log_level_from_index(double u_min, double rho, int index) {
    if (index < 0) return 0.0;
    return u_min / pow(rho, (double)index);
}

int qc_log_index_from_value(double x, double u_min, double rho) {
    double abs_x = fabs(x);
    if (abs_x < u_min) return -1;
    double idx = log(abs_x / u_min) / log(1.0 / rho);
    return (int)round(idx);
}

/* ================================================================
 * Vector Quantizer
 * ================================================================ */

int qc_vec_codebook_init(QCVectorCodebook *codebook, int K, int dim) {
    if (!codebook || K <= 0 || dim <= 0 || dim > QC_MAX_STATE_DIM) return -1;
    memset(codebook, 0, sizeof(QCVectorCodebook));
    codebook->num_codewords = K;
    codebook->dimension = dim;
    codebook->entries = calloc(K, sizeof(QCCodebookEntry));
    codebook->codebook_flat = calloc(K * dim, sizeof(double));
    if (!codebook->entries || !codebook->codebook_flat) {
        free(codebook->entries);
        free(codebook->codebook_flat);
        return -1;
    }
    for (int k = 0; k < K; k++) {
        codebook->entries[k].index = k;
        codebook->entries[k].dimension = dim;
        codebook->entries[k].vector = &codebook->codebook_flat[k * dim];
        codebook->entries[k].distortion = 0.0;
        codebook->entries[k].usage_count = 0;
    }
    codebook->trained = 0;
    return 0;
}

int qc_vec_codebook_train_lbg(QCVectorCodebook *codebook,
                               const double *data, int num_vectors,
                               int max_iter, double epsilon) {
    if (!codebook || !data || num_vectors <= 0 || codebook->num_codewords <= 0) return -1;
    int K = codebook->num_codewords, dim = codebook->dimension;
    double prev_distortion = INFINITY;
    int *assignments = calloc(num_vectors, sizeof(int));
    double *centroid_sum = calloc(K * dim, sizeof(double));
    int *centroid_count = calloc(K, sizeof(int));
    if (!assignments || !centroid_sum || !centroid_count) {
        free(assignments); free(centroid_sum); free(centroid_count); return -1;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        /* Assignment step */
        double total_dist = 0.0;
        memset(centroid_sum, 0, K * dim * sizeof(double));
        memset(centroid_count, 0, K * sizeof(int));
        for (int n = 0; n < num_vectors; n++) {
            const double *x = &data[n * dim];
            int best_k = 0;
            double best_dist = INFINITY;
            for (int k = 0; k < K; k++) {
                double dist = 0.0;
                for (int d = 0; d < dim; d++) {
                    double diff = x[d] - codebook->codebook_flat[k * dim + d];
                    dist += diff * diff;
                }
                if (dist < best_dist) { best_dist = dist; best_k = k; }
            }
            assignments[n] = best_k;
            total_dist += best_dist;
            centroid_count[best_k]++;
            for (int d = 0; d < dim; d++) {
                centroid_sum[best_k * dim + d] += x[d];
            }
        }
        total_dist /= num_vectors;

        /* Update step */
        for (int k = 0; k < K; k++) {
            if (centroid_count[k] > 0) {
                for (int d = 0; d < dim; d++) {
                    codebook->codebook_flat[k * dim + d] =
                        centroid_sum[k * dim + d] / centroid_count[k];
                }
            }
        }

        if (fabs(prev_distortion - total_dist) < epsilon) break;
        prev_distortion = total_dist;
    }

    codebook->total_distortion = prev_distortion;
    codebook->trained = 1;

    free(assignments); free(centroid_sum); free(centroid_count);
    return 0;
}

int qc_vec_quantize(const QCVectorCodebook *codebook, const double *x,
                     int *index, double *distortion) {
    if (!codebook || !x || !index) return -1;
    int K = codebook->num_codewords, dim = codebook->dimension;
    int best_k = 0;
    double best_dist = INFINITY;
    for (int k = 0; k < K; k++) {
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = x[d] - codebook->codebook_flat[k * dim + d];
            dist += diff * diff;
        }
        if (dist < best_dist) { best_dist = dist; best_k = k; }
    }
    *index = best_k;
    if (distortion) *distortion = best_dist;
    return 0;
}

int qc_vec_reconstruct(const QCVectorCodebook *codebook, int index,
                        double *x_reconstructed) {
    if (!codebook || !x_reconstructed || index < 0 || index >= codebook->num_codewords) return -1;
    int dim = codebook->dimension;
    for (int d = 0; d < dim; d++) {
        x_reconstructed[d] = codebook->codebook_flat[index * dim + d];
    }
    return 0;
}

void qc_vec_codebook_free(QCVectorCodebook *codebook) {
    if (!codebook) return;
    free(codebook->entries);
    free(codebook->codebook_flat);
    codebook->entries = NULL;
    codebook->codebook_flat = NULL;
}

/* ================================================================
 * Dither Quantizer (TPDF)
 * ================================================================ */

double qc_dither_tpdf_sample(void) {
    /* TPDF: sum of two independent uniform random variables */
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    return u1 + u2 - 1.0;
}

double qc_dither_quantize_subtractive(double x, double step) {
    double dither = qc_dither_tpdf_sample() * step;
    double x_dithered = x + dither;
    double xq = qc_uniform_quantize_midtread(x_dithered, step);
    return xq - dither;
}

double qc_dither_quantize_nonsubtractive(double x, double step) {
    double dither = qc_dither_tpdf_sample() * step;
    return qc_uniform_quantize_midtread(x + dither, step);
}

/* ================================================================
 * Quantizer Performance Metrics
 * ================================================================ */

void qc_quantizer_metrics_compute(const QCQuantizer *q,
                                   const double *signal,
                                   const double *quantized,
                                   int length,
                                   QCQuantizerMetrics *metrics) {
    if (!signal || !quantized || !metrics || length <= 0) return;
    memset(metrics, 0, sizeof(QCQuantizerMetrics));
    double sum_sq = 0.0, sum_err_sq = 0.0, sum_err_abs = 0.0;
    int overloads = 0;
    for (int i = 0; i < length; i++) {
        sum_sq += signal[i] * signal[i];
        double err = signal[i] - quantized[i];
        sum_err_sq += err * err;
        sum_err_abs += fabs(err);
        if (q && fabs(signal[i]) > q->range_max) overloads++;
    }
    double sig_power = sum_sq / length;
    double noise_power = sum_err_sq / length;
    metrics->snr_db = (noise_power > 0) ? 10.0 * log10(sig_power / noise_power) : INFINITY;
    metrics->sqnr_db = metrics->snr_db;
    metrics->mse = noise_power;
    metrics->mae = sum_err_abs / length;
    metrics->enob = qc_enob_from_snr(metrics->snr_db);
    metrics->overload_percent = 100.0 * overloads / length;
    metrics->num_samples = length;
}

double qc_enob_from_snr(double snr_db) {
    if (snr_db <= 0) return 0;
    return (snr_db - 1.76) / 6.02;
}

double qc_snr_from_enob(double enob) {
    return 6.02 * enob + 1.76;
}

double qc_uniform_sqnr_theoretical(int bits, double signal_std,
                                    double range_half) {
    double step = 2.0 * range_half / (double)((1 << bits) - 1);
    double noise_var = step * step / 12.0;
    double sig_var = signal_std * signal_std;
    if (noise_var <= 0) return INFINITY;
    return 10.0 * log10(sig_var / noise_var);
}

/* ================================================================
 * Lloyd-Max Optimal Quantizer for Gaussian Source
 * ================================================================ */

double qc_lloyd_max_gaussian(int bits, double *levels,
                              double *boundaries, int max_iter) {
    int N = 1 << bits;
    if (N < 2 || !levels || !boundaries || max_iter <= 0) return -1.0;

    /* Initialize levels uniformly over [-4, 4] (covers ~99.99% of Gaussian) */
    double range = 8.0;
    for (int i = 0; i < N; i++) {
        levels[i] = -range/2.0 + range * (i + 0.5) / N;
        if (i < N - 1) boundaries[i] = -range/2.0 + range * (i + 1.0) / N;
    }

    double prev_mse = INFINITY;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Update boundaries: midpoint between adjacent levels */
        for (int i = 0; i < N - 1; i++) {
            boundaries[i] = (levels[i] + levels[i + 1]) / 2.0;
        }

        /* Update levels: conditional expectation E[X | X in cell_i] */
        double mse = 0.0;
        for (int i = 0; i < N; i++) {
            double lo = (i == 0) ? -INFINITY : boundaries[i - 1];
            double hi = (i == N - 1) ? INFINITY : boundaries[i];

            /* Approximate E[X | lo < X < hi] using numerical integration */
            double num = 0.0, den = 0.0;
            int steps = 100;
            double dx = (hi - lo) / steps;
            if (isinf(lo)) { lo = -8.0; dx = (hi - lo) / steps; }
            if (isinf(hi)) { hi = 8.0; dx = (hi - lo) / steps; }

            for (int s = 0; s < steps; s++) {
                double x = lo + (s + 0.5) * dx;
                double pdf = exp(-x * x / 2.0) / sqrt(2.0 * M_PI);
                num += x * pdf * dx;
                den += pdf * dx;
            }
            if (den > 1e-15) levels[i] = num / den;

            for (int s = 0; s < steps; s++) {
                double x = lo + (s + 0.5) * dx;
                double pdf = exp(-x * x / 2.0) / sqrt(2.0 * M_PI);
                mse += (x - levels[i]) * (x - levels[i]) * pdf * dx;
            }
        }

        if (fabs(prev_mse - mse) < 1e-6) break;
        prev_mse = mse;
    }
    return prev_mse;
}

/* ================================================================
 * Quantization Noise Analysis
 * ================================================================ */

double qc_quantization_noise_variance(double step) {
    return step * step / 12.0;
}

double qc_quantization_error_autocorr(const double *signal,
                                       const double *quantized,
                                       int length, int lag) {
    if (!signal || !quantized || length <= 0 || lag < 0 || lag >= length) return 0.0;
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < length - lag; i++) {
        double e1 = signal[i] - quantized[i];
        double e2 = signal[i + lag] - quantized[i + lag];
        sum += e1 * e2;
        count++;
    }
    return (count > 0) ? sum / count : 0.0;
}

int qc_quantization_error_is_white(const double *signal,
                                    const double *quantized,
                                    int length, int max_lag,
                                    double significance) {
    if (!signal || !quantized || length <= 0 || max_lag <= 0) return 0;
    double var0 = qc_quantization_error_autocorr(signal, quantized, length, 0);
    if (var0 <= 0.0) return 1;
    for (int lag = 1; lag <= max_lag && lag < length; lag++) {
        double corr = qc_quantization_error_autocorr(signal, quantized, length, lag) / var0;
        if (fabs(corr) > significance) return 0;
    }
    return 1;
}

double qc_signal_crest_factor(const double *signal, int length) {
    if (!signal || length <= 0) return 0.0;
    double peak = 0.0, rms = 0.0;
    for (int i = 0; i < length; i++) {
        double abs_val = fabs(signal[i]);
        if (abs_val > peak) peak = abs_val;
        rms += signal[i] * signal[i];
    }
    rms = sqrt(rms / length);
    return (rms > 0.0) ? peak / rms : 0.0;
}

double qc_non_uniform_optimal_compander(double x, double sigma, int bits) {
    /* Mu-law companding for non-uniform quantization */
    double mu = 255.0;
    double y = copysign(log(1.0 + mu * fabs(x) / sigma) / log(1.0 + mu), x);
    double step = 2.0 / (double)((1 << bits) - 1);
    double yq = step * round(y / step);
    return yq;
}

double qc_compander_expand(double y, double sigma, int bits) {
    (void)bits; /* reserved for mu-law parameterization */
    double mu = 255.0;
    double abs_y = fabs(y);
    double x = sigma / mu * (pow(1.0 + mu, abs_y) - 1.0);
    return copysign(x, y);
}
