#ifndef LYAPUNOV_KRASOVSKII_H
#define LYAPUNOV_KRASOVSKII_H

#include "time_delay_system.h"

/* ============================================================================
 * Lyapunov-Krasovskii (LK) Functionals for Time-Delay Systems
 *
 * Reference:
 *   N. N. Krasovskii, "Stability of Motion" (1963) — original LK theorem
 *   J. K. Hale & S. M. Verduyn Lunel (1993) — modern formulation
 *   K. Gu, V. L. Kharitonov, J. Chen (2003) — discretized LK functional
 *   E. Fridman, "Introduction to Time-Delay Systems" (2014) — LMI approach
 *
 * Key distinction:
 *   Lyapunov FUNCTION:  V(t, x(t))          — finite dimensional
 *   Lyapunov-Krasovskii FUNCTIONAL: V(t, x_t) — on C([-τ,0], Rⁿ)
 *
 * Level 4 — Fundamental Theorem
 * ============================================================================ */

/* ============================================================================
 * LK Functional Types
 * ============================================================================ */

/* --- Standard LK functional form ---
 * V(t, x_t) = xᵀ(t) P x(t)                       ← point delay term
 *           + ∫_{t-τ}^{t} xᵀ(s) Q x(s) ds        ← interval integral
 *           + ∫_{-τ}^{0} ∫_{t+θ}^{t} ẋᵀ(s) R ẋ(s) ds dθ  ← derivative term
 *
 * Lyapunov-Krasovskii Theorem:
 *   If there exist P > 0, Q > 0, R > 0 such that
 *   dV/dt ≤ -ε||x(t)||² along trajectories,
 *   then the trivial solution is uniformly asymptotically stable
 *   for all constant delays 0 ≤ τ ≤ τ*. */

typedef struct {
    double* P;        /* n×n Lyapunov matrix (point term) */
    double* Q;        /* n×n matrix for integral term */
    double* R;        /* n×n matrix for derivative term */
    double* S;        /* n×n additional weighting matrix */
    int n;
} LKFunctional;

/* --- Augmented LK functional (for time-varying delay) ---
 * V(t, x_t) = V_standard(t, x_t)
 *           + ∫_{t-τ(t)}^{t} xᵀ(s) S1 x(s) ds
 *           + (τ_max - τ(t)) ∫ ...  additional slack terms
 * This allows delay-range-dependent stability analysis. */
typedef struct {
    LKFunctional base;
    double* Z1;       /* Augmented weighting matrices */
    double* Z2;
    double* N1;       /* Slack matrices for descriptor approach */
    double* N2;
    double* M1;       /* Free-weighting matrices */
    double* M2;
    double* M3;
    double tau_min;   /* Lower delay bound */
    double tau_max;   /* Upper delay bound */
    double mu;        /* Bound on |dτ/dt| */
} AugmentedLKFunctional;

/* --- Discretized LK functional (Gu-Kharitonov-Chen) ---
 *
 * V(x_t) = xᵀ(t) P x(t) + 2xᵀ(t) ∫_{-τ}^{0} Q(ξ) x(t+ξ) dξ
 *        + ∫_{-τ}^{0} ∫_{-τ}^{0} xᵀ(t+ξ) R(ξ,η) x(t+η) dη dξ
 *        + ∫_{-τ}^{0} xᵀ(t+ξ) S(ξ) x(t+ξ) dξ
 *
 * Where Q, R, S are piecewise linear matrix functions on [-τ, 0].
 * Discretized into N mesh points. */
typedef struct {
    int N;                /* Number of mesh points (N+1 intervals) */
    double tau;           /* Delay value */
    double* P;            /* n×n */
    double* Q_mesh;       /* (N+1)×n×n matrix function Q(ξ_i) */
    double* R_mesh;       /* (N+1)×(N+1)×n×n matrix function R(ξ_i, η_j) */
    double* S_mesh;       /* (N+1)×n×n matrix function S(ξ_i) */
    int n;
} DiscretizedLKFunctional;

/* ============================================================================
 * LK Functional Construction and Evaluation
 * ============================================================================ */

/* Allocate a standard LK functional */
LKFunctional* lkf_create(int n);

/* Free an LK functional */
void lkf_free(LKFunctional* lkf);

/* Set the P matrix from a given array */
void lkf_set_P(LKFunctional* lkf, const double* P_data);

/* Set the Q matrix from a given array */
void lkf_set_Q(LKFunctional* lkf, const double* Q_data);

/* Set the R matrix from a given array */
void lkf_set_R(LKFunctional* lkf, const double* R_data);

/* Generate P, Q, R as identity matrices (trivial solution for verification) */
void lkf_set_identity(LKFunctional* lkf);

/* Evaluate V(t, x_t) = xᵀPx + ∫xᵀQx + ∫∫ẋᵀRẋ
 * Requires: x_current (n×1), x_history (n×N_hist) over [-τ, 0],
 *           xdot_history (n×N_hist) */
double lkf_evaluate(const LKFunctional* lkf,
                    const double* x_current,
                    const double* x_history,
                    const double* xdot_history,
                    int n_hist, double tau, double dt);

/* Evaluate the derivative dV/dt along system trajectories.
 * For linear system: ẋ = A x + A_d x(t-τ),
 * dV/dt = 2xᵀP(Ax+A_d x_d) + xᵀQx - x_dᵀQx_d + τẋᵀRẋ - ∫_{t-τ}^{t} ẋᵀRẋ ds */
double lkf_derivative(const LKFunctional* lkf,
                      const TimeDelaySystem* sys,
                      const double* x_current,
                      const double* x_delayed,
                      const double* xdot_current,
                      const double* xdot_history,
                      int n_hist, double tau, double dt);

/* Check if P, Q, R are all positive definite (eigenvalue check).
 * Returns true if all eigenvalues > 0. */
bool lkf_is_positive_definite(const LKFunctional* lkf);

/* Check if dV/dt ≤ -ε||x||² (the negative-definite rate condition).
 * Returns the computed decay rate ε. */
double lkf_decay_rate(const LKFunctional* lkf,
                      const TimeDelaySystem* sys,
                      const double* x_current,
                      const double* x_delayed,
                      double tau);

/* --- Augmented LK --- */
AugmentedLKFunctional* alkf_create(int n, double tau_min, double tau_max, double mu);
void alkf_free(AugmentedLKFunctional* alkf);
double alkf_evaluate(const AugmentedLKFunctional* alkf,
                     const double* x_current,
                     const double* x_history,
                     const double* xdot_history,
                     int n_hist, double tau_t, double dt);

/* --- Discretized LK --- */
DiscretizedLKFunctional* dlkf_create(int n, int N, double tau);
void dlkf_free(DiscretizedLKFunctional* dlkf);
double dlkf_evaluate(const DiscretizedLKFunctional* dlkf,
                     const double* x_current,
                     const double* x_history,
                     int n_hist, double dt);

/* ============================================================================
 * Lyapunov-Razumikhin (LR) Theorem
 * ============================================================================ */

/* Razumikhin condition: V̇(x(t)) ≤ -w(||x(t)||) whenever
 *   V(x(t+θ)) ≤ p V(x(t)) for θ ∈ [-τ, 0] and p > 1.
 *
 * This theorem is often easier to apply than LK for time-varying delays.
 * It uses a standard Lyapunov FUNCTION (not functional) with
 * an additional "Razumikhin condition" on the history. */

/* Check Razumikhin condition for a given p > 1 */
bool razumikhin_condition_check(const double* x_history, int n_hist,
                                const double* x_current, int n,
                                double p_factor, double (*V_func)(const double*, int));

/* ============================================================================
 * LMI-Based Stability Check (using matrix conditions)
 * ============================================================================ */

/* Check the standard LMI condition for delay-independent stability:
 *
 * [ AᵀP + PA + Q    PA_d  ]
 * [   A_dᵀP        -Q    ]  < 0
 *
 * If feasible with P>0, Q>0, system is stable for ALL τ ≥ 0. */
bool lmi_delay_independent_check(const TimeDelaySystem* sys,
                                 double* out_P, double* out_Q);

/* Check delay-dependent LMI:
 *
 * [ AᵀP+PA+Q-τR    PA_d+τR    τAᵀR ]
 * [ A_dᵀP+τR       -Q-τR      τA_dᵀR ] < 0
 * [ τRA            τRA_d      -τR   ]
 *
 * With τ = tau_max, if feasible, system stable for τ ∈ [0, τ_max]. */
bool lmi_delay_dependent_check(const TimeDelaySystem* sys,
                                double tau_max, double* out_P,
                                double* out_Q, double* out_R);

#endif /* LYAPUNOV_KRASOVSKII_H */
