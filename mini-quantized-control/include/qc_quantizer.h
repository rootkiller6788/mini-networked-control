/**
 * @file    qc_quantizer.h
 * @brief   Quantizer types, structures and low-level operations
 *
 * Detailed quantizer implementations:
 *   - Uniform quantizer (mid-rise and mid-tread)
 *   - Logarithmic quantizer (Elia-Mitter formulation)
 *   - Vector quantizer (lattice-based LBG)
 *   - Dither quantizer (subtractive and non-subtractive)
 *
 * Key references:
 *   - Widrow, B. & Kollar, I. (2008). Quantization Noise. Cambridge.
 *   - Gray, R.M. & Neuhoff, D.L. (1998). Quantization. IEEE TIT.
 *   - Gersho, A. & Gray, R.M. (1992). Vector Quantization. Kluwer.
 */

#ifndef QC_QUANTIZER_H
#define QC_QUANTIZER_H

#include "quantized_control.h"
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    QC_UNIFORM_MID_TREAD = 0,
    QC_UNIFORM_MID_RISE  = 1
} QCUniformSubtype;

typedef enum {
    QC_DITHER_SUBTRACTIVE     = 0,
    QC_DITHER_NON_SUBTRACTIVE = 1
} QCDitherType;

typedef struct {
    double      snr_db;
    double      sqnr_db;
    double      sfdr_db;
    double      enob;
    double      overload_percent;
    double      mse;
    double      mae;
    int         num_samples;
} QCQuantizerMetrics;

typedef struct {
    int         index;
    double     *vector;
    int         dimension;
    double      distortion;
    int         usage_count;
} QCCodebookEntry;

typedef struct {
    int             num_codewords;
    int             dimension;
    QCCodebookEntry *entries;
    double         *codebook_flat;
    double          total_distortion;
    int             trained;
} QCVectorCodebook;

typedef struct {
    QCQuantizerType type;
    int             bits;
    double          range_min;
    double          range_max;
    QCUniformSubtype uniform_subtype;
    double          log_rho;
    double          log_delta;
    QCDitherType    dither_type;
    double          dither_variance;
    int             vector_dim;
    QCOverloadStrategy overload;
} QCQuantizerConfig;

/* Uniform quantizer operations */
double qc_uniform_quantize_midtread(double x, double step);
double qc_uniform_quantize_midrise(double x, double step);
double qc_uniform_optimal_step(double range, int bits);
double qc_uniform_overload_prob(double range, double sigma);

/* Logarithmic quantizer operations */
double qc_log_rho_from_delta(double delta);
int qc_log_num_levels_from_range(double u_min, double u_max, double rho);
double qc_log_deadzone_from_umin(double u_min);
double qc_log_level_from_index(double u_min, double rho, int index);
int qc_log_index_from_value(double x, double u_min, double rho);

/* Vector quantizer operations */
int qc_vec_codebook_init(QCVectorCodebook *codebook, int K, int dim);
int qc_vec_codebook_train_lbg(QCVectorCodebook *codebook,
                               const double *data, int num_vectors,
                               int max_iter, double epsilon);
int qc_vec_quantize(const QCVectorCodebook *codebook, const double *x,
                     int *index, double *distortion);
int qc_vec_reconstruct(const QCVectorCodebook *codebook, int index,
                        double *x_reconstructed);
void qc_vec_codebook_free(QCVectorCodebook *codebook);

/* Dither quantizer operations */
double qc_dither_tpdf_sample(void);
double qc_dither_quantize_subtractive(double x, double step);
double qc_dither_quantize_nonsubtractive(double x, double step);

/* Performance analysis */
void qc_quantizer_metrics_compute(const QCQuantizer *q,
                                   const double *signal,
                                   const double *quantized,
                                   int length,
                                   QCQuantizerMetrics *metrics);
double qc_enob_from_snr(double snr_db);
double qc_snr_from_enob(double enob);
double qc_uniform_sqnr_theoretical(int bits, double signal_std,
                                    double range_half);
double qc_lloyd_max_gaussian(int bits, double *levels,
                              double *boundaries, int max_iter);
double qc_quantization_noise_variance(double step);
double qc_quantization_error_autocorr(const double *signal,
                                       const double *quantized,
                                       int length, int lag);
int qc_quantization_error_is_white(const double *signal,
                                    const double *quantized,
                                    int length, int max_lag,
                                    double significance);
double qc_signal_crest_factor(const double *signal, int length);
double qc_non_uniform_optimal_compander(double x, double sigma, int bits);
double qc_compander_expand(double y, double sigma, int bits);

#ifdef __cplusplus
}
#endif

#endif /* QC_QUANTIZER_H */
