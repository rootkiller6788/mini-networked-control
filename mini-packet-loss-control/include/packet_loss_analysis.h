#ifndef PACKET_LOSS_ANALYSIS_H
#define PACKET_LOSS_ANALYSIS_H

#include "packet_loss_core.h"
#include <stdbool.h>

/* ============================================================================
 * Packet Loss Analysis — Stability & Performance Under Packet Dropout
 *
 * Analytical tools for assessing stability, performance bounds, and
 * critical thresholds in networked control systems with packet loss.
 *
 * Key theoretical frameworks:
 *   - Stochastic Lyapunov stability (mean-square, almost-sure)
 *   - Jump linear systems (Markovian jump linear systems, MJLS)
 *   - Expected LQG cost analysis
 *   - Critical loss probability computation
 *
 * Key references:
 *   - Costa, Fragoso, Marques (2005):
 *     "Discrete-Time Markov Jump Linear Systems"
 *   - Seiler & Sengupta (2005): "An H∞ Approach to Networked Control"
 *   - You & Xie (2013): "Survey of Stability of Networked Control Systems"
 *   - Donkers, Heemels, van de Wouw, Hetel (2011):
 *     "Stability Analysis of Networked Control Systems Using a Switched
 *     Linear Systems Approach"
 * ============================================================================ */

/* --- Stability Definitions --- */

/**
 * Stability notions for stochastic NCS with packet loss.
 *
 * For a system x_{k+1} = A_γ(k) x_k where γ(k) ∈ {0,1} is the
 * packet arrival indicator (random), different stability concepts apply:
 *
 * - MEAN_SQUARE: E[||x_k||²] → 0 as k→∞
 * - ALMOST_SURE: P(lim_{k→∞} ||x_k|| = 0) = 1
 * - STOCHASTIC: ∀ε>0, lim_{k→∞} P(||x_k|| > ε) = 0
 * - EXPONENTIAL_MS: ∃c>0, α∈(0,1): E[||x_k||²] ≤ c α^k E[||x_0||²]
 *
 * For Bernoulli loss: Mean-Square ⇔ Exponential-MS (for linear systems).
 *
 * Reference: Kozin (1969), "A Survey of Stability of Stochastic Systems"
 */
typedef enum {
    STAB_MEAN_SQUARE = 0,
    STAB_ALMOST_SURE = 1,
    STAB_STOCHASTIC = 2,
    STAB_EXPONENTIAL_MS = 3
} StabilityType;

/**
 * Stability certificate: result of stability analysis.
 */
typedef struct {
    StabilityType type;
    bool is_stable;                      /* Stability verdict */
    double margin;                       /* Stability margin (≥0 if stable) */
    double decay_rate;                   /* Estimated decay rate */
    double critical_loss_prob;           /* Maximum tolerable loss probability */
    double lyapunov_value;               /* Lyapunov function value at test point */
    int iterations_to_check;             /* Steps simulated */
} StabilityCertificate;

/* --- Jump Linear System --- */

/**
 * Markovian Jump Linear System (MJLS) model.
 *
 * Represents an NCS where the closed-loop dynamics depend on the
 * packet delivery status, modeled as a finite-state Markov chain θ_k.
 *
 * System: x_{k+1} = A(θ_k) x_k
 *
 * where θ_k ∈ {0, 1, ..., M-1} is a Markov chain with
 * transition matrix P[i][j] = P(θ_{k+1}=j | θ_k=i).
 *
 * Each mode i corresponds to a specific delivery pattern:
 *   e.g., Mode 0 = both packets arrived,
 *         Mode 1 = sensor lost,
 *         Mode 2 = actuator lost,
 *         Mode 3 = both lost.
 *
 * MJLS Mean-Square Stability (MSS) condition (Costa et al., 2005):
 * The system is MSS iff the spectral radius of the augmented operator
 * ρ(A) < 1, where A = (P' ⊗ I) diag(A₀⊗A₀, ..., A_M-1⊗A_M-1).
 *
 * Complexity: O(M·n⁶) for exact test. Efficiently approximated via
 * coupled Lyapunov equations.
 */
typedef struct {
    int n_states;                   /* State dimension */
    int n_modes;                    /* Number of Markov modes */
    double** mode_matrices;         /* A[i]: n×n matrix for mode i (row-major) */
    double** transition_matrix;     /* P: M×M Markov transition matrix */
    double* steady_state;           /* π: Stationary distribution of θ */

    /* Lyapunov matrices per mode */
    double** lyapunov_matrices;     /* Pi: n×n Lyapunov matrix per mode */

    /* Diagnostics */
    double spectral_radius_estimate;
    bool mss_certified;
    int iterations_to_converge;
} JumpLinearSystem;

/* --- Critical Probability Analysis --- */

/**
 * Critical loss probability analysis context.
 *
 * For a given system, computes the maximum tolerable packet loss
 * probability p_c such that the system is stable for p < p_c.
 */
typedef struct {
    /* System parameters */
    double spectral_radius_A;       /* ρ(A) — open-loop spectral radius */
    double spectral_radius_Acl;     /* ρ(A-BL) — closed-loop spectral radius */

    /* Critical probabilities */
    double p_c_sensor;              /* Max sensor loss for stable estimation */
    double p_c_actuator;            /* Max actuator loss for stable control */
    double p_c_joint;               /* Joint critical probability */

    /* Bounds and approximations */
    double lower_bound_sensor;      /* 1 - 1/ρ(A)² */
    double upper_bound_sensor;      /* From observable mode analysis */
    double lower_bound_actuator;    /* From controllability analysis */

    /* Stability region: set of (p_s, p_a) where system is stable */
    double* stability_region;       /* Sampled boundary points */
    int n_region_points;

    /* Eigenvalue-based diagnostics */
    double max_eigenvalue_product;   /* Maximum eigenvalue product across modes */
    bool is_stabilizable_under_loss;
} CriticalProbabilityAnalysis;

/* ============================================================================
 * Stability Analysis API
 * ============================================================================ */

/**
 * Test mean-square stability for Bernoulli packet loss.
 *
 * For a system x_{k+1} = A_s x_k (success) or A_f x_k (failure/loss),
 * with P(success) = 1-p, P(failure) = p:
 *
 * The system is MSS iff:
 *   ρ( (1-p)·(A_s ⊗ A_s) + p·(A_f ⊗ A_f) ) < 1
 *
 * where ⊗ is the Kronecker product. This is a necessary and sufficient
 * condition (Costa et al., 2005).
 *
 * Complexity: O(n⁶) via full eigenvalue computation of n²×n² matrix.
 * For n ≤ 10, this is tractable. For larger n, use approximate method.
 */
StabilityCertificate* pl_stability_test_bernoulli(
    const double* A_success, const double* A_failure,
    int n, double loss_prob);

/**
 * Test mean-square stability for Gilbert-Elliott channel.
 *
 * Uses coupled Lyapunov equations:
 *
 * For i ∈ {0,1} (Good/Bad states):
 *   P_i - Σ_j P_ij · A_j' P_j A_j ≻ 0
 *
 * Feasibility of these LMIs (Linear Matrix Inequalities) certifies MSS.
 *
 * Approximated by iterating the Lyapunov recursion until convergence
 * or divergence.
 */
StabilityCertificate* pl_stability_test_gilbert_elliott(
    const double* A_success, const double* A_failure,
    int n, const GilbertElliottChannel* ch);

/**
 * Generic stability test using MJLS framework for arbitrary
 * Markovian packet loss models.
 */
StabilityCertificate* pl_stability_test_markov(
    JumpLinearSystem* jls);

/**
 * Estimate the worst-case decay rate (Lyapunov exponent) of the system.
 * λ = lim_{k→∞} (1/k) · E[log ||x_k||]
 *
 * For stable systems, λ < 0.
 * For marginally stable systems, λ ≈ 0.
 * For unstable systems, λ > 0.
 */
double pl_stability_lyapunov_exponent(
    const double* A_success, const double* A_failure,
    int n, double loss_prob, int n_trials, int n_steps,
    unsigned long seed);

/**
 * Monte Carlo simulation of mean-square error over time.
 * Estimates E[||x_k||²] by averaging over n_trials independent runs.
 *
 * Returns average final ||x||².
 */
double pl_stability_monte_carlo(
    const double* A_success, const double* A_failure,
    int n, double loss_prob, int n_trials, int n_steps,
    unsigned long seed);

/**
 * Compute the stochastic Lyapunov function value at a given state.
 * V(x) = x' P x where P is the solution to the coupled Lyapunov equations.
 *
 * For Bernoulli loss:
 *   P = (1-p)·A_s' P A_s + p·A_f' P A_f + I  (Lyapunov iteration)
 *
 * If a positive definite solution P exists, V(x) is a valid
 * stochastic Lyapunov function.
 */
double pl_stability_lyapunov_function(
    const double* x, int n,
    const double* A_success, const double* A_failure,
    double loss_prob, int max_iter, double tol);

/* --- Jump Linear System API --- */

/**
 * Create a Jump Linear System model.
 * @param n_states: State vector dimension
 * @param n_modes: Number of Markov chain modes (delivery patterns)
 */
JumpLinearSystem* pl_jls_create(int n_states, int n_modes);
void pl_jls_free(JumpLinearSystem* jls);

/** Set the system matrix for mode i (n×n, row-major). */
void pl_jls_set_mode_matrix(JumpLinearSystem* jls, int mode,
                             const double* A_mode);

/** Set one entry in the Markov transition matrix. */
void pl_jls_set_transition(JumpLinearSystem* jls, int from, int to,
                            double prob);

/** Compute steady-state distribution of the Markov chain. */
void pl_jls_compute_steady_state(JumpLinearSystem* jls);

/**
 * Test MSS using the coupled Lyapunov equation approach.
 *
 * Iteratively solve:
 *   P_i^{(t+1)} = Σ_j p_ij · A_j' P_j^{(t)} A_j + I
 *
 * starting from P_i^{(0)} = I. If the sequence converges, the
 * system is MSS. The spectral radius of the operator A determines
 * the convergence rate.
 *
 * Reference: Costa & Fragoso (1993), "Stability Results for
 * Discrete-Time Linear Systems with Markovian Jumping Parameters"
 */
bool pl_jls_test_mss(JumpLinearSystem* jls, int max_iter, double tol);

/** Compute the augmented operator spectral radius for MJLS MSS test. */
double pl_jls_operator_spectral_radius(JumpLinearSystem* jls);

void pl_jls_print(const JumpLinearSystem* jls);

/* --- Critical Probability Analysis API --- */

/**
 * Compute critical loss probabilities for a given LTI system.
 *
 * For sensor-side loss (estimation):
 *   γ_c ≥ 1 - 1/ρ(A)²  (necessary condition, Sinopoli 2004)
 *   γ_c ≤ depends on observability structure
 *
 * For actuator-side loss (control):
 *   p_c ≤ 1 - 1/ρ(A_cl)²  where A_cl = A - BL
 *
 * Joint critical region: stability requires both sensor and actuator
 * loss rates to be within a convex region (often triangular).
 */
CriticalProbabilityAnalysis* pl_critical_prob_analyze(
    const double* A, const double* B, const double* C,
    const double* L, int n, int m, int p);
void pl_critical_prob_free(CriticalProbabilityAnalysis* cpa);

/**
 * Sample points on the stability boundary in (p_sensor, p_actuator) space.
 * Returns number of boundary points generated.
 */
int pl_critical_prob_stability_region(CriticalProbabilityAnalysis* cpa,
                                       int n_samples);

void pl_critical_prob_print(const CriticalProbabilityAnalysis* cpa);

/* --- Eigenvalue & Matrix Analysis Utilities --- */

/**
 * Compute Kronecker product: C = A ⊗ B.
 * A: m×n, B: p×q, C: (m·p)×(n·q).
 *
 * Kronecker product property: (A ⊗ B)(C ⊗ D) = (AC) ⊗ (BD)
 * Used in MJLS stability condition: the system matrix is constructed
 * from Kronecker products of mode matrices.
 */
void pl_kronecker_product(const double* A, const double* B,
                           double* C, int m, int n, int p, int q);

/**
 * Compute eigenvalues of a real matrix using the QR algorithm.
 * For n ≤ 32, uses Francis QR with double shifts.
 * Returns number of real eigenvalues found (complex come in conjugate pairs).
 *
 * Complexity: O(n³) for the full Schur decomposition.
 */
int pl_eigenvalues(const double* A, double* real_part, double* imag_part, int n);

/**
 * Compute spectral radius: ρ(A) = max_i |λ_i(A)|.
 * Uses power iteration for the dominant eigenvalue.
 * Complexity: O(n² · iterations).
 */
double pl_spectral_radius_power(const double* A, int n, int max_iter, double tol);

/**
 * Matrix norm: Frobenius norm ||A||_F = sqrt(Σ_{i,j} a_{ij}²).
 */
double pl_matrix_norm_frobenius(const double* A, int rows, int cols);

/**
 * Check if a matrix is Schur stable (ρ(A) < 1).
 */
bool pl_is_schur_stable(const double* A, int n);

#endif /* PACKET_LOSS_ANALYSIS_H */