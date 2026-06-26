/**
 * @file    quantized_control.h
 * @brief   Core Quantized Control Theory definitions and API
 *
 * Quantized control theory addresses control systems where measurements
 * and/or control signals are subject to magnitude discretization (quantization).
 * This is fundamental to networked control systems (NCS) where communication
 * bandwidth is limited and data must be represented with finite precision.
 *
 * Key pioneers: Brockett (2000), Liberzon (2003), Elia & Mitter (2001),
 * Nair & Evans (2004), Tatikonda & Mitter (2004), Fu & Xie (2005).
 *
 * Course mappings:
 *   - MIT 6.241J Dynamic Systems and Control (quantization effects)
 *   - Stanford AA 203 Optimal and Learning-Based Control
 *   - Berkeley EE 222 Nonlinear Systems (quantized nonlinear control)
 *   - Caltech CDS 110 Introduction to Control (data-rate limited control)
 *   - ETH 227-0216 System Identification (quantized measurements)
 *   - Cambridge 4F3 Nonlinear Control (quantizer design)
 *   - Oxford B4 Predictive Control (limited-precision MPC)
 *   - CMU 24-677 Nonlinear Control (quantized feedback)
 *   - Princeton MAE 546 Optimal Control (information-constrained control)
 *
 * Key references:
 *   - Elia, N. & Mitter, S.K. (2001). Stabilization of linear systems with
 *     limited information. IEEE Transactions on Automatic Control.
 *   - Nair, G.N. & Evans, R.J. (2004). Stabilizability of stochastic linear
 *     systems with finite feedback data rates. SIAM J. Control Optim.
 *   - Liberzon, D. (2003). Switching in Systems and Control. Birkhauser.
 *   - Tatikonda, S. & Mitter, S. (2004). Control under communication
 *     constraints. IEEE Transactions on Automatic Control.
 *   - Fu, M. & Xie, L. (2005). The sector bound approach to quantized
 *     feedback control. IEEE Transactions on Automatic Control.
 *   - Brockett, R.W. & Liberzon, D. (2000). Quantized feedback stabilization
 *     of linear systems. IEEE Transactions on Automatic Control.
 */

#ifndef QUANTIZED_CONTROL_H
#define QUANTIZED_CONTROL_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QC_MAX_STATE_DIM           32
#define QC_MAX_INPUT_DIM           16
#define QC_MAX_OUTPUT_DIM          16
#define QC_MAX_BITS                64
#define QC_MAX_QUANTIZER_REGIONS  4096
#define QC_EPSILON                 1e-9
#define QC_PI                      3.14159265358979323846
#define QC_MAX_ITER                10000

/* ===================================================================
 * L1: Core Definitions
 * =================================================================== */

typedef enum {
    QC_QTYPE_UNIFORM        = 0,
    QC_QTYPE_LOGARITHMIC    = 1,
    QC_QTYPE_DYNAMIC        = 2,
    QC_QTYPE_VECTOR         = 3,
    QC_QTYPE_RANDOM_DITHER  = 4,
    QC_QTYPE_FIXED_RATE     = 5,
    QC_QTYPE_VARIABLE_RATE  = 6
} QCQuantizerType;

typedef enum {
    QC_OVERLOAD_SATURATE   = 0,
    QC_OVERLOAD_ZOOM_OUT   = 1,
    QC_OVERLOAD_MODULO     = 2,
    QC_OVERLOAD_EXTEND     = 3
} QCOverloadStrategy;

typedef enum {
    QC_ENC_FIXED_LENGTH  = 0,
    QC_ENC_HUFFMAN       = 1,
    QC_ENC_DIFFERENTIAL  = 2,
    QC_ENC_ADAPTIVE      = 3,
    QC_ENC_ARITHMETIC    = 4
} QCEncodingScheme;

typedef enum {
    QC_STABLE_ASYMPTOTIC     = 0,
    QC_STABLE_PRACTICAL      = 1,
    QC_STABLE_MARGINAL       = 2,
    QC_STABLE_UNSTABLE       = 3,
    QC_STABLE_ULTIMATE_BOUND = 4
} QCStabilityStatus;

typedef struct {
    int     dim;
    double  x[QC_MAX_STATE_DIM];
    double  q_error_bound[QC_MAX_STATE_DIM];
    double  t;
} QCState;

typedef struct {
    int     dim;
    double  u[QC_MAX_INPUT_DIM];
    double  u_quantized[QC_MAX_INPUT_DIM];
    int     bits_per_component;
    double  resolution;
    double  q_error_norm;
} QCControl;

typedef struct {
    int     dim;
    double  y[QC_MAX_OUTPUT_DIM];
    double  y_quantized[QC_MAX_OUTPUT_DIM];
    int     bits_used;
    double  snr;
    double  noise_power;
} QCOutput;

typedef struct {
    int     index;
    int     dimension;
    double  center[QC_MAX_STATE_DIM];
    double  half_width[QC_MAX_STATE_DIM];
    int     code;
    double  volume;
} QCQuantizationCell;

typedef struct {
    QCQuantizerType     type;
    int                 bits;
    int                 num_levels;
    double              range_min;
    double              range_max;
    double              step;
    double              rho;
    double              delta;
    double             *levels;
    int                 levels_len;
    QCOverloadStrategy  overload;
    double              zoom_factor;
    int                 is_logarithmic;
    double              saturation_value;
    double              resolution_normalized;
} QCQuantizer;

typedef struct {
    double              rho;
    double              u_min;
    int                 num_positive;
    int                 num_negative;
    double             *levels_positive;
    double             *levels_negative;
    double              deadzone;
    double              sector_delta;
    int                 total_levels;
} QCLogQuantizer;

typedef struct {
    double              mu;
    double              mu_min;
    double              mu_max;
    double              rho_in;
    double              rho_out;
    double              M;
    int                 bits;
    double              switching_threshold;
    int                 zoom_level;
    int                 max_zoom_level;
    int                 zoom_in_count;
    int                 zoom_out_count;
    double              dwell_time;
    double              time_since_last_zoom;
} QCDynamicQuantizer;

typedef struct {
    QCEncodingScheme    scheme;
    int                 bits_per_symbol;
    int                 buffer_size;
    uint8_t            *buffer;
    int                 buffer_pos;
    int                 total_bits_encoded;
    int                 bit_offset;
    double             *symbol_probs;
    int                *huffman_codes;
    int                *huffman_lengths;
    int                 alphabet_size;
    double              prev_value;
    int                 diff_initialized;
} QCEncoder;

typedef struct {
    QCEncodingScheme    scheme;
    int                 bits_per_symbol;
    const uint8_t      *buffer;
    int                 buffer_size;
    int                 buffer_pos;
    int                 total_bits_decoded;
    int                 bit_offset;
    int                *huffman_codes;
    int                *huffman_lengths;
    double             *codebook;
    int                 alphabet_size;
    double              prev_value;
    int                 diff_initialized;
} QCDecoder;

typedef struct {
    double              total_rate_bps;
    int                 num_channels;
    int                *channel_bits;
    double             *channel_rate;
    double             *eig_magnitudes;
    int                 num_unstable;
    double              min_rate;
    double              rate_margin;
    double              channel_capacity;
    double              entropy_rate;
} QCDataRate;

typedef struct {
    int                 state_dim;
    int                 input_dim;
    int                 output_dim;
    double              A[QC_MAX_STATE_DIM * QC_MAX_STATE_DIM];
    double              B[QC_MAX_STATE_DIM * QC_MAX_INPUT_DIM];
    double              C[QC_MAX_OUTPUT_DIM * QC_MAX_STATE_DIM];
    double              D[QC_MAX_OUTPUT_DIM * QC_MAX_INPUT_DIM];
    QCQuantizer         input_quantizer;
    QCQuantizer         output_quantizer;
    QCEncoder           input_encoder;
    QCDecoder            input_decoder;
    QCEncoder           output_encoder;
    QCDecoder            output_decoder;
    QCDataRate          data_rate;
    double              sector_delta;
    double              sector_upper;
    double              sampling_period;
    QCStabilityStatus   stability;
} QCSystem;

typedef struct {
    int                 converged;
    int                 steps;
    double              final_error;
    double              max_quantization_error;
    double              avg_bit_rate;
    double              stability_margin;
    double             *state_trajectory;
    double             *input_trajectory;
    double             *output_trajectory;
    double             *time_points;
    int                 trajectory_len;
    int                 trajectory_capacity;
    QCStabilityStatus   final_status;
    double              lyapunov_value;
} QCSimulationResult;

typedef struct {
    int                 state_dim;
    double              eigenvalues_real[QC_MAX_STATE_DIM];
    double              eigenvalues_imag[QC_MAX_STATE_DIM];
    int                 num_eigenvalues;
    double              min_rate_bits;
    double              actual_rate;
    int                 is_stabilizable;
    double              rate_gap;
    double              instability_sum;
} QCDataRateTheorem;

typedef struct {
    double              delta;
    int                 is_sector_bounded;
    double              L2_gain_bound;
    int                 nominally_stable;
    int                 robustly_stable;
    double              feasibility_margin;
    double              circle_criterion_margin;
} QCSectorBoundResult;

typedef struct {
    int                 zoom_sequence_len;
    int                *zoom_directions;
    double             *lyapunov_values;
    double             *state_norms;
    double              convergence_rate;
    int                 bounded;
    int                 converges_to_zero;
    int                 num_zoom_outs;
    int                 num_zoom_ins;
    double              max_state_norm;
} QCZoomVerification;

/* API Functions */
void qc_system_init(QCSystem *sys);
int qc_system_configure(QCSystem *sys, int nx, int nu, int ny);
int qc_system_set_A(QCSystem *sys, const double *A);
int qc_system_set_B(QCSystem *sys, const double *B);
int qc_system_set_C(QCSystem *sys, const double *C);
int qc_system_set_D(QCSystem *sys, const double *D);
int qc_system_validate(const QCSystem *sys);
void qc_system_print(const QCSystem *sys);
size_t qc_system_memory_estimate(const QCSystem *sys);

void qc_quantizer_init(QCQuantizer *q, QCQuantizerType type, int bits);
int qc_quantizer_configure_range(QCQuantizer *q, double min_val, double max_val);
int qc_quantizer_set_log_params(QCQuantizer *q, double rho, double delta);
double qc_quantize_scalar(const QCQuantizer *q, double x);
void qc_quantize_vector(const QCQuantizer *q, const double *x, int dim, double *xq);
double qc_quantization_error(const QCQuantizer *q, double x);
int qc_quantizer_num_levels(const QCQuantizer *q);
double qc_quantizer_snr(const QCQuantizer *q, double signal_power);
void qc_quantizer_free(QCQuantizer *q);
const char* qc_quantizer_type_name(QCQuantizerType t);
int qc_quantizer_find_cell(const QCQuantizer *q, double x, QCQuantizationCell *cell);
double qc_quantizer_max_error(const QCQuantizer *q);

void qc_log_quantizer_init(QCLogQuantizer *lq, double rho, double u_min);
double qc_log_quantize(const QCLogQuantizer *lq, double x);
double qc_log_sector_delta(double rho);
double qc_log_quantizer_density(const QCLogQuantizer *lq);
int qc_log_quantizer_find_level(const QCLogQuantizer *lq, double x, int *is_positive, int *level_index);
void qc_log_quantizer_free(QCLogQuantizer *lq);
int qc_log_quantizer_build_levels(QCLogQuantizer *lq, int num_levels);

void qc_dyn_quantizer_init(QCDynamicQuantizer *dq, int bits, double M);
double qc_dyn_quantize(QCDynamicQuantizer *dq, double x);
void qc_dyn_zoom_in(QCDynamicQuantizer *dq);
void qc_dyn_zoom_out(QCDynamicQuantizer *dq);
int qc_dyn_should_zoom(const QCDynamicQuantizer *dq, double x);
double qc_dyn_get_range(const QCDynamicQuantizer *dq);
void qc_dyn_reset(QCDynamicQuantizer *dq);
double qc_dyn_get_resolution(const QCDynamicQuantizer *dq);

void qc_encoder_init(QCEncoder *enc, QCEncodingScheme scheme, int bits_per_sym);
int qc_encoder_encode_scalar(QCEncoder *enc, const QCQuantizer *q, double value);
int qc_encoder_encode_vector(QCEncoder *enc, const QCQuantizer *q, const double *values, int dim);
int qc_encoder_build_huffman(QCEncoder *enc, const double *probs, int alphabet_sz);
void qc_encoder_reset(QCEncoder *enc);
void qc_encoder_free(QCEncoder *enc);
size_t qc_encoder_bytes_used(const QCEncoder *enc);
int qc_encoder_total_bits(const QCEncoder *enc);

void qc_decoder_init(QCDecoder *dec, QCEncodingScheme scheme, int bits_per_sym);
double qc_decoder_decode_scalar(QCDecoder *dec, const QCQuantizer *q);
int qc_decoder_decode_vector(QCDecoder *dec, const QCQuantizer *q, double *values, int dim);
int qc_decoder_set_codebook(QCDecoder *dec, const int *codes, const int *lengths, const double *levels, int alphabet_sz);
void qc_decoder_reset(QCDecoder *dec);
void qc_decoder_free(QCDecoder *dec);

void qc_data_rate_init(QCDataRate *dr, double total_rate, int num_ch);
int qc_data_rate_compute_min(QCDataRate *dr, const double *A, int n);
int qc_data_rate_allocate(QCDataRate *dr);
int qc_data_rate_check_stabilizability(const QCDataRate *dr);
double qc_data_rate_entropy_power(const double *signal, int len, double var);
int qc_data_rate_optimal_allocation(QCDataRate *dr, const double *ch_vars, int num_ch);
double qc_data_rate_theoretical_min(const double *A, int n);
double qc_rate_distortion_gaussian(double variance, double distortion);
int qc_data_rate_channel_feasible(double capacity, double required_rate);

void qc_sector_bound_init(QCSectorBoundResult *sr);
int qc_sector_bound_analyze(const QCSystem *sys, QCSectorBoundResult *sr);
double qc_sector_bound_delta_for_rho(double rho);
int qc_sector_bound_feasibility(const QCSystem *sys, double *feasibility_margin);
int qc_circle_criterion_check(const QCSystem *sys, double delta, double *margin);

void qc_zoom_verify_init(QCZoomVerification *zv);
int qc_zoom_verify_step(QCZoomVerification *zv, double V_now, double V_prev, double state_norm);
int qc_zoom_verify_conclusion(QCZoomVerification *zv);
void qc_zoom_verify_free(QCZoomVerification *zv);

void qc_sim_result_init(QCSimulationResult *res);
void qc_sim_result_free(QCSimulationResult *res);
int qc_simulate_closed_loop(QCSystem *sys, const double *x0, double t0, double tf, double dt, double (*controller)(const double*, int, double, double*), QCSimulationResult *res);
int qc_simulate_quantized_lqr(QCSystem *sys, const double *x0, const double *K, double t0, double tf, double dt, QCSimulationResult *res);
int qc_simulate_output_feedback(QCSystem *sys, const double *x0, const double *K, const double *L, double t0, double tf, double dt, QCSimulationResult *res);
void qc_sim_result_print(const QCSimulationResult *res);

int qc_matrix_eigenvalues(const double *A, int n, double *eig_real, double *eig_imag);
int qc_is_schur_stable(const double *A, int n);
int qc_is_hurwitz_stable(const double *A, int n);
double qc_spectral_radius(const double *A, int n);
void qc_matrix_multiply(const double *A, const double *B, double *C, int m, int n, int p);
int qc_solve_lyapunov(const double *A, const double *Q, double *P, int n);
int qc_control_with_saturation(const QCSystem *sys, const double *x, double *u, const double *u_max);
double qc_ultimate_bound_estimate(const QCSystem *sys);
const char* qc_overload_strategy_name(QCOverloadStrategy s);
const char* qc_stability_status_name(QCStabilityStatus s);
const char* qc_encoding_scheme_name(QCEncodingScheme e);
double qc_bit_error_probability(double snr_db);

#ifdef __cplusplus
}
#endif

#endif /* QUANTIZED_CONTROL_H */
