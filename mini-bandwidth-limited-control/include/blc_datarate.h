/**
 * blc_datarate.h — Data Rate Theorem and Information-Theoretic Control Limits
 *
 * This header implements the mathematical framework connecting Shannon's
 * information theory to control-theoretic stabilizability.
 *
 * Fundamental Results:
 *
 * 1. Scalar Data Rate Theorem (Wong & Brockett, 1997, 1999):
 *    For ẋ = λx + u with λ > 0, the minimum bit rate is
 *       R > λ / ln(2)  bits/sec
 *
 *    For discrete-time x[k+1] = a·x[k] + u[k] with |a| > 1:
 *       R > log₂|a|  bits/sample
 *
 * 2. Multi-Dimensional Data Rate Theorem (Nair & Evans, 2000, 2004):
 *    For ẋ = Ax + Bu with eigenvalues λᵢ(A):
 *       R > Σ_{i: Re(λᵢ)>0} 2·Re(λᵢ) / ln(2)  bits/sec (continuous-time)
 *       R > Σ_{i: |λᵢ|>1} log₂|λᵢ|  bits/sample (discrete-time)
 *
 * 3. Quantization-Error-Induced Instability (Baillieul, 2002):
 *    Under limited data rate, quantization error acts as a disturbance.
 *    If the bit rate is below the critical threshold, no controller
 *    (linear or nonlinear) can stabilize the system.
 *
 * 4. Separation under Bandwidth Limits:
 *    Unlike classical LQG, estimation and control do NOT separate
 *    under rate constraints — the encoder must be co-designed with
 *    the controller (Tatikonda, Sahai, Mitter, 2004).
 *
 * Knowledge coverage: L4 (Fundamental Theorems), L5 (Analytical Methods)
 */

#ifndef BLC_DATARATE_H
#define BLC_DATARATE_H

#include "blc_core.h"
#include <complex.h>

/* ================================================================
 * Zoom Quantizer — Dynamic Range Adjustment
 *
 * The zoom quantizer (Brockett & Liberzon, 2000) dynamically
 * adjusts its range based on the observed state. If the state
 * approaches the quantization boundary, the quantizer "zooms out";
 * as the state converges to the origin, it "zooms in" for higher
 * precision.
 *
 * Parameters:
 *   - ρ_in  > 1: zoom-in factor (contraction)
 *   - ρ_out > 1: zoom-out factor (expansion)
 *   - M: number of quantization levels
 *   - L: initial range parameter
 *
 * Algorithm:
 *   1. Encode state x(k) relative to current range r(k)
 *   2. If |x(k)| ≤ r(k)/ρ_in:  zoom in:  r(k+1) = r(k) / ρ_in
 *   3. If |x(k)| > r(k):       zoom out: r(k+1) = ρ_out * r(k)
 *   4. Otherwise:              hold:     r(k+1) = r(k)
 *
 * This achieves asymptotic stability for linear systems with
 * R > Σ log₂|λᵢ(A)| bits per sample despite using finite
 * quantization levels.
 * ================================================================ */

/** Zoom quantizer state machine */
typedef struct {
    double  range;              /** Current quantization range [-range, +range] */
    double  rho_in;             /** Zoom-in contraction factor (> 1) */
    double  rho_out;            /** Zoom-out expansion factor (> 1) */
    int     levels;             /** Number of quantization levels M */
    double  step;               /** Current step size Δ = 2*range / levels */
    int     last_index;         /** Last transmitted quantization index */
    int     zoom_in_count;      /** Count of consecutive zoom-in operations */
    int     zoom_out_count;     /** Count of consecutive zoom-out operations */
    double  min_range;          /** Minimum achievable range */
    bool    is_zooming;         /** Currently zooming flag */
    double  convergence_rate;   /** Estimated convergence rate */
} BLCZoomQuantizer;

/** Predictive encoder state
 *
 *  Predictive coding (Nair & Evans, 2004): instead of transmitting
 *  the raw state, transmit the prediction error.
 *
 *  At the controller, an identical predictor runs:
 *    ẋ̂(t) = A·ẍ̂(t) + B·u(t)    (open-loop prediction)
 *
 *  The sensor transmits the quantized prediction error:
 *    e(k) = x(k·T_s) - ẍ̂(k·T_s)
 *
 *  This exploits the fact that prediction errors have smaller
 *  magnitude than the raw state when the model is good, reducing
 *  the required quantization range for a given bit budget.
 */
typedef struct {
    double  x_pred[BLC_MAX_STATES];    /** Predicted state ẍ̂(k) */
    double  x_hat[BLC_MAX_STATES];     /** Filtered estimate after decoding */
    double  prediction_error[BLC_MAX_STATES]; /** e(k) = x - x_pred */
    double  error_norm;                /** ||e(k)|| */
    double  prediction_gain;           /** Error variance reduction ratio */
    int     prediction_steps;          /** Horizon of prediction */
} BLCEncoder;

/** Delta modulator state
 *
 *  Delta modulation transmits only the sign of the change
 *  between successive samples (±Δ), using only 1 bit per sample.
 *  This is the extreme case of bandwidth-limited control.
 *
 *  The reconstruction follows: ẍ̂[k] = ẍ̂[k-1] + Δ·sgn(x[k] - ẍ̂[k-1])
 *
 *  Slope overload occurs when |dx/dt| > Δ/T_s.
 */
typedef struct {
    double  step_size;           /** Delta step size Δ */
    double  x_hat;              /** Current reconstructed value */
    int     bit_sequence[64];    /** Record of transmitted bits */
    int     seq_len;            /** Length of bit sequence */
    int     slope_overloads;    /** Count of slope overload events */
} BLCDeltaModulator;

/* ================================================================
 * Information-Theoretic Bounds for Control
 * ================================================================ */

/** Compute the theoretical minimum data rate for stabilization.
 *  @param plant The plant with computed eigenvalues
 *  @return bits per second (continuous-time) */
double blc_datarate_min_ct(const BLCPlant* plant);

/** Compute the theoretical minimum data rate per sample.
 *  For discrete-time systems: R_min = Σ log₂|λᵢ| for |λᵢ| > 1
 *  @param plant The plant model
 *  @return bits per sample */
double blc_datarate_min_dt(const BLCPlant* plant);

/** Compute the practical minimum data rate accounting for
 *  quantization overhead, packet headers, and encoding efficiency.
 *  @param plant The plant model
 *  @param Ts Sampling period
 *  @param efficiency Encoding efficiency factor (0 < η ≤ 1)
 *  @return practical minimum bits per second */
double blc_datarate_practical(const BLCPlant* plant, double Ts,
                               double efficiency);

/** Check if a given data rate is sufficient for stabilization.
 *  @param bit_rate Available bit rate (bits/sec)
 *  @param plant The plant model
 *  @param Ts Sampling period
 *  @return true if rate is sufficient */
bool   blc_datarate_is_sufficient(double bit_rate, const BLCPlant* plant,
                                   double Ts);

/** Compute the "bit-rate gap" — how many bits/sec below the minimum.
 *  Positive gap means insufficient rate. */
double blc_datarate_gap(double bit_rate, const BLCPlant* plant, double Ts);

/** Rate-Distortion function for control: R(D) minimal bit rate to
 *  achieve a given distortion D = E[||x - ẍ̂||²].
 *  For a scalar unstable pole λ > 0:
 *    R(D) = max(0, (1/2)·log₂(σ_w² / D) + λ/ln(2))
 *  where σ_w² is the process noise variance.
 *
 *  @param plant The plant model
 *  @param distortion Target distortion level D
 *  @param noise_var Process noise variance σ_w²
 *  @return minimal bits/sec to achieve distortion D */
double blc_rate_distortion(const BLCPlant* plant, double distortion,
                            double noise_var);

/** Distortion-Rate function (inverse of rate-distortion):
 *  Minimum achievable distortion given bit rate R.
 *  D(R) = σ_w² · 2^{-2(R - λ/ln(2))}  for R > λ/ln(2)
 *
 *  @param plant The plant model
 *  @param bit_rate Available bit rate
 *  @param noise_var Process noise variance σ_w²
 *  @return minimum achievable distortion */
double blc_distortion_rate(const BLCPlant* plant, double bit_rate,
                            double noise_var);

/* ================================================================
 * Zoom Quantizer API
 * ================================================================ */

/** Initialize zoom quantizer */
int     blc_zoom_init(BLCZoomQuantizer* zq, double initial_range,
                       double rho_in, double rho_out, int levels);

/** Encode a value using zoom quantizer, updating the range */
int     blc_zoom_encode(BLCZoomQuantizer* zq, double value);

/** Decode a zoom-quantized index */
double  blc_zoom_decode(const BLCZoomQuantizer* zq, int index);

/** Get current zoom quantizer step size */
double  blc_zoom_get_step(const BLCZoomQuantizer* zq);

/** Get zoom quantizer range */
double  blc_zoom_get_range(const BLCZoomQuantizer* zq);

/** Check if zoom quantizer is in zoom-in phase */
bool    blc_zoom_is_zooming_in(const BLCZoomQuantizer* zq);

/** Reset zoom quantizer to initial state */
void    blc_zoom_reset(BLCZoomQuantizer* zq, double new_range);

/** Get zoom statistics */
void    blc_zoom_stats(const BLCZoomQuantizer* zq, int* zoom_in,
                        int* zoom_out, double* conv_rate);

/* ================================================================
 * Predictive Encoder API
 * ================================================================ */

/** Initialize predictive encoder */
int     blc_encoder_init(BLCEncoder* enc, const double* x0_pred,
                          int prediction_horizon);

/** Encode state using prediction + quantization */
int     blc_encoder_encode(BLCEncoder* enc, const BLCQuantizer* q,
                            const double* true_state, const BLCPlant* plant,
                            double dt, BLCPacket* pkt, int* bit_count);

/** Decode received packet, update prediction */
int     blc_encoder_decode(BLCEncoder* enc, const BLCQuantizer* q,
                            const double* u, const BLCPlant* plant,
                            double dt, const BLCPacket* pkt);

/** Get current state estimate */
void    blc_encoder_get_estimate(const BLCEncoder* enc, double* x_hat);

/** Get prediction error norm */
double  blc_encoder_get_error_norm(const BLCEncoder* enc);

/* ================================================================
 * Delta Modulation API
 * ================================================================ */

/** Initialize delta modulator */
int     blc_delta_init(BLCDeltaModulator* dm, double step_size,
                        double initial_value);

/** Encode one sample with delta modulation (returns 1 bit: +Δ or -Δ) */
int     blc_delta_encode(BLCDeltaModulator* dm, double value);

/** Decode delta-modulated bit */
double  blc_delta_decode(BLCDeltaModulator* dm, int bit);

/** Get reconstructed signal */
double  blc_delta_get_value(const BLCDeltaModulator* dm);

/** Check for slope overload */
bool    blc_delta_slope_overload(const BLCDeltaModulator* dm,
                                 double input_derivative, double Ts);

/** Get slope overload count */
int     blc_delta_overload_count(const BLCDeltaModulator* dm);

int     blc_delta_reset(BLCDeltaModulator* dm, double step, double init);

/* ================================================================
 * Spectral Analysis for Data Rate Computation
 * ================================================================ */

/** Compute eigenvalues of a real matrix using QR algorithm
 *  (for state dimension up to BLC_MAX_STATES).
 *  @param A Input matrix (n×n, row-major from BLCPlant)
 *  @param n Dimension
 *  @param re_real Real parts of eigenvalues (output)
 *  @param re_imag Imaginary parts of eigenvalues (output)
 *  @return Number of iterations to convergence */
int     blc_eigenvalues(const double* A, int n, double* re_real,
                         double* re_imag);

/** Compute singular values of a real matrix (SVD via Golub-Reinsch).
 *  @param A Input m×n matrix
 *  @param m Rows
 *  @param n Columns
 *  @param s Output singular values (descending) */
int     blc_singular_values(const double* A, int m, int n, double* s);

/** Compute matrix exponential exp(A*dt) using Pade approximation.
 *  @param A Input n×n matrix
 *  @param n Dimension
 *  @param dt Time step
 *  @param expA Output exp(A·dt) */
int     blc_matrix_exp(const double* A, int n, double dt, double* expA);

/** Compute matrix logarithm via inverse scaling-and-squaring */
int     blc_matrix_log(const double* A, int n, double* logA);

/** Lyapunov equation solver: A*P + P*A' = -Q
 *  Bartels-Stewart algorithm for continuous-time Lyapunov equation.
 *  @param A System matrix (n×n)
 *  @param Q Right-hand side (n×n, symmetric)
 *  @param n Dimension
 *  @param P Solution (n×n, output) */
int     blc_lyapunov_solve(const double* A, const double* Q, int n, double* P);

/** Compute spectral radius ρ(A) = max|λᵢ(A)| */
double  blc_spectral_radius(const double* A, int n);

/** Compute spectral abscissa α(A) = max Re(λᵢ(A)) */
double  blc_spectral_abscissa(const double* A, int n);

/** Compute matrix norm ||A||₂ (using SVD, equals max singular value) */
double  blc_matrix_norm_2(const double* A, int m, int n);

#endif /* BLC_DATARATE_H */