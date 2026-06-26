#ifndef DELAY_STABILITY_H
#define DELAY_STABILITY_H

#include "time_delay_system.h"

/* ============================================================================
 * Stability Analysis Tools for Time-Delay Systems
 *
 * Reference:
 *   K. Gu, V. Kharitonov, J. Chen, "Stability of Time-Delay Systems" (2003)
 *   S.-I. Niculescu, "Delay Effects on Stability" (2001)
 *   W. Michiels & S.-I. Niculescu, "Stability and Stabilization of
 *     Time-Delay Systems" (2007)
 *   G. Stépán, "Retarded Dynamical Systems" (1989)
 *
 * Level 4 — Fundamental Laws
 * Level 5 — Stability Analysis Algorithms
 * ============================================================================ */

/* ============================================================================
 * Frequency-Domain Stability Analysis
 * ============================================================================ */

/* --- Nyquist Criterion for Delay Systems ---
 *
 * For ẋ = A x + A_d x(t-τ):
 * The characteristic quasipolynomial is
 *   Δ(s) = det(sI - A - A_d e^{-τs})
 *
 * Nyquist contour: evaluate Δ(jω) for ω ∈ [0, ∞) and count encirclements.
 * The system is stable iff all roots of Δ(s) have Re < 0.
 *
 * For scalar system with delay G(s) = G₀(s)e^{-τs}:
 *   The delay adds phase lag e^{-jωτ} without changing magnitude.
 *   Phase margin = π - ∠G₀(jω_c) - ω_c·τ  */

/* Compute the Nyquist curve points for a delay system.
 * Returns number of frequency points. */
int delay_nyquist_points(const TimeDelaySystem* sys,
                         double w_min, double w_max, int n_points,
                         double* out_omega, double* out_real, double* out_imag);

/* Find the gain margin for a scalar delay system */
double delay_gain_margin(const TimeDelaySystem* sys);

/* Find the phase margin for a scalar delay system */
double delay_phase_margin(const TimeDelaySystem* sys, double tau);

/* --- Delay Margin Computation --- */

/* Find τ* — the maximum constant delay for which the system
 * remains stable. Uses the frequency-sweeping method:
 * Find all ω such that Δ(jω, τ) = 0 has a solution for some τ.
 * τ* = min{θ/ω : Δ(jω) = 0 for some ω ≥ 0, θ ∈ [0, 2π)}. */
double delay_margin_frequency_sweep(const TimeDelaySystem* sys);

/* Find delay margin via the matrix pencil method.
 * For ẋ = A x + A_d x(t-τ), compute the generalized eigenvalues
 * of the matrix pair that determine stability crossings. */
double delay_margin_matrix_pencil(const TimeDelaySystem* sys);

/* --- Stability Crossing Detection ---
 *
 * For fixed τ, count the number of RHP roots using the
 * Walton-Marshall method based on the Rekasius substitution
 * e^{-jωτ} ≈ (1 - jωT) / (1 + jωT). */

/* Count RHP roots using Rekasius substitution */
int delay_rhp_root_count(const TimeDelaySystem* sys, double tau);

/* Check stability for a specific delay value */
bool delay_is_stable_at(const TimeDelaySystem* sys, double tau);

/* ============================================================================
 * Time-Domain Stability Criteria
 * ============================================================================ */

/* --- Matrix Measure (Logarithmic Norm) Approach ---
 *
 * μ(A) = lim_{ε→0⁺} (||I + εA|| - 1)/ε = λ_max((A+Aᵀ)/2)  for 2-norm.
 * Sufficient condition for delay-independent stability:
 *   μ(A) + ||A_d|| < 0  */

/* Compute the matrix measure μ_2(A) = λ_max((A+Aᵀ)/2) */
double matrix_measure_l2(const double* A, int n);

/* Check delay-independent stability via matrix measure */
bool matrix_measure_stability_check(const TimeDelaySystem* sys);

/* --- Halanay Inequality ---
 *
 * For scalar V̇(t) ≤ -α V(t) + β sup_{t-τ ≤ s ≤ t} V(s):
 * If α > β > 0, then V(t) → 0 exponentially.
 * This yields a sufficient condition for decay of the LK functional. */

/* Verify Halanay condition and compute decay rate */
bool halanay_check(double alpha, double beta, double* decay_rate);

/* ============================================================================
 * Padé Approximation of Time Delay
 * ============================================================================ */

/* The Padé approximation replaces e^{-τs} with a rational transfer function:
 *
 *   e^{-τs} ≈ P_{N,M}(-τs) / P_{N,M}(τs)
 *
 * where P_{N,M}(x) = Σ_{k=0}^{N} (N+M-k)! N! / ((N+M)! k! (N-k)!) x^k
 *
 * Common approximations:
 *   1st-order: e^{-τs} ≈ (1 - τs/2) / (1 + τs/2)
 *   2nd-order: e^{-τs} ≈ (1 - τs/2 + τ²s²/12) / (1 + τs/2 + τ²s²/12) */

typedef struct {
    int order;           /* Approximation order (1-10) */
    double tau;          /* Delay value */
    double* num;         /* Numerator coefficients (order+1) */
    double* den;         /* Denominator coefficients (order+1) */
    int num_size;        /* = order + 1 */
    int den_size;        /* = order + 1 */
} PadeApproximation;

/* Compute Padé coefficients for e^{-τs} approximation of given order */
PadeApproximation* pade_create(int order, double tau);
void pade_free(PadeApproximation* pade);

/* Convert Padé approximation to state-space form:
 *   ẋ_p = A_p x_p + B_p u
 *   y_p = C_p x_p + D_p u  */
void pade_to_state_space(const PadeApproximation* pade,
                         double** out_A, double** out_B,
                         double** out_C, double* out_D);

/* Apply Padé approximation to augment a system:
 * Replace delay block e^{-τs} in series with the plant model. */
typedef struct {
    double* A_aug;       /* Augmented system matrix (n+p)×(n+p) */
    double* B_aug;
    double* C_aug;
    int n_aug;           /* n + Padé order */
} PadeAugmentedSystem;

PadeAugmentedSystem* pade_augment_system(const TimeDelaySystem* sys,
                                          int pade_order);
void pade_augmented_free(PadeAugmentedSystem* pas);

/* ============================================================================
 * Frequency-Sweeping Test for Interval Delays
 * ============================================================================ */

/* For τ ∈ [τ_min, τ_max], check stability using
 * the sweeping test: count stability crossings in the frequency domain. */

/* Result of a sweeping test */
typedef struct {
    double tau_critical;     /* Smallest crossing delay */
    bool is_stable_min;      /* Stable at τ = τ_min? */
    bool is_stable_max;      /* Stable at τ = τ_max? */
    int n_crossings;         /* Number of stability crossings */
    double* crossing_taus;   /* Delay values at crossings */
    int crossing_count;
} SweepingTestResult;

SweepingTestResult* delay_sweeping_test(const TimeDelaySystem* sys,
                                         double tau_min, double tau_max,
                                         int n_freq_points);
void sweeping_test_free(SweepingTestResult* result);

/* ============================================================================
 * Exponential Stability Rate
 * ============================================================================ */

/* Compute the exponential decay rate α such that
 * ||x(t)|| ≤ M e^{-αt} ||φ||_c for all trajectories.
 * α = -max{Re(λ) : characteristic roots}. */
double delay_exponential_decay_rate(const TimeDelaySystem* sys);

/* Check if the system is exponentially stable (α > 0) */
bool delay_is_exponentially_stable(const TimeDelaySystem* sys);

#endif /* DELAY_STABILITY_H */
