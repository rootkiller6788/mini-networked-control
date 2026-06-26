#ifndef PACKET_LOSS_ESTIMATOR_H
#define PACKET_LOSS_ESTIMATOR_H

#include "packet_loss_core.h"
#include <stdbool.h>

/* ============================================================================
 * Packet Loss Estimator — Kalman Filter with Intermittent Observations
 *
 * Implements estimation algorithms for systems where sensor measurements
 * arrive intermittently due to packet loss.
 *
 * Key references:
 *   - Sinopoli, Schenato, Franceschetti, Poolla, Jordan, Sastry (2004):
 *     "Kalman Filtering with Intermittent Observations"
 *     IEEE Trans. Automatic Control, 49(9):1453-1464.
 *   - Plarre & Bullo (2007): "On Kalman Filtering with Intermittent
 *     Observations: Critical Value for Second Moment Stability"
 *   - Mo & Sinopoli (2012): "Kalman Filtering with Intermittent
 *     Observations: Tail Distribution and Critical Value"
 *   - Kar, Sinopoli, Moura (2012): "Kalman Filtering with Intermittent
 *     Observations: Weak Convergence to a Stationary Distribution"
 * ============================================================================ */

/* --- Kalman Filter Types --- */

/**
 * Arrival process for sensor measurements.
 * γ_k = 1 if measurement y_k arrives at time k, 0 otherwise.
 */
typedef struct {
    int* arrivals;              /* Binary sequence: 1=arrived, 0=lost */
    int length;                 /* Length of sequence */
    int capacity;
    double arrival_rate;        /* Empirical fraction of arrivals */
    int max_consecutive_losses; /* Worst-case consecutive losses observed */
} ArrivalProcess;

/**
 * Standard Kalman filter state (for systems with perfect communication).
 *
 * Process:  x_{k+1} = A x_k + w_k,   w_k ~ N(0, Q)
 * Measure:   y_k   = C x_k + v_k,    v_k ~ N(0, R)
 *
 * Predict:  x̂_{k|k-1} = A x̂_{k-1|k-1}
 *           P_{k|k-1} = A P_{k-1|k-1} A' + Q
 * Update:   K_k = P_{k|k-1} C' (C P_{k|k-1} C' + R)^{-1}
 *           x̂_{k|k} = x̂_{k|k-1} + K_k (y_k - C x̂_{k|k-1})
 *           P_{k|k} = (I - K_k C) P_{k|k-1}
 */
typedef struct {
    /* System matrices (references, not owned) */
    const double* A;     /* n×n */
    const double* C;     /* p×n */
    const double* Q;     /* n×n */
    const double* R;     /* p×p */
    int n;               /* State dimension */
    int p;               /* Measurement dimension */

    /* Filter state */
    double* x_hat;       /* Current state estimate x̂_{k|k}, size n */
    double* P;           /* Estimation error covariance P_{k|k}, n×n */
    double* x_pred;      /* Predicted state x̂_{k|k-1}, size n */
    double* P_pred;      /* Predicted covariance P_{k|k-1}, n×n */
    double* K;           /* Kalman gain matrix, n×p */

    /* Scratch buffers */
    double* temp_nn;     /* n×n temporary */
    double* temp_np;     /* n×p temporary */
    double* temp_pp;     /* p×p temporary */
    double* innovation;  /* y - C x̂_pred, size p */

    int k;               /* Time step counter */
} KalmanFilter;

/**
 * Kalman filter with intermittent observations.
 *
 * Key modification (Sinopoli et al., 2004):
 * When γ_k = 1 (measurement arrived):
 *   Standard Kalman update with y_k.
 * When γ_k = 0 (measurement lost):
 *   Skip the update step — only time-update:
 *   x̂_{k|k} = x̂_{k|k-1},  P_{k|k} = P_{k|k-1}.
 *
 * The estimation error covariance evolves as a random process:
 *   P_{k+1} = A P_k A' + Q - γ_k A P_k C' (C P_k C' + R)^{-1} C P_k A'
 *
 * This is a random algebraic Riccati iteration.
 *
 * Critical probability γ_c: For γ_k Bernoulli i.i.d. with P(γ_k=1) = γ,
 * there exists a threshold γ_c ∈ [0,1) such that:
 *   - If γ > γ_c: E[P_k] is bounded ∀k (second-moment stable)
 *   - If γ ≤ γ_c: E[P_k] diverges as k→∞
 *
 * Bounds on γ_c (Sinopoli 2004, Theorem 2):
 *   1 - 1/ρ(A)² ≤ γ_c ≤ 1 - 1/(max_i |λ_i(A)|)² · (worst observable mode)
 */
typedef struct {
    /* Base Kalman filter state */
    KalmanFilter kf;

    /* Intermittent observation tracking */
    ArrivalProcess arrivals;
    int consecutive_losses;     /* Current run of lost measurements */
    int max_consecutive_losses; /* Maximum consecutive losses ever observed */

    /* Expected covariance (Monte Carlo estimate) */
    double* E_P;                /* Expected error covariance, n×n */
    double trace_E_P;           /* Trace of expected covariance */

    /* Critical probability estimation */
    double arrival_probability; /* γ = P(γ_k = 1) */
    double critical_gamma_lower; /* Lower bound on γ_c */
    double critical_gamma_upper; /* Upper bound on γ_c */

    /* Stability flags */
    bool is_mean_stable;        /* Is E[P_k] bounded? */
    bool is_covariance_stable;  /* Is E[||P_k||²] bounded? */
} IntermittentKalmanFilter;

/* --- Bounded-Error Estimation Under Loss --- */

/**
 * Set-membership / bounded-error estimator for systems with packet loss.
 *
 * Instead of stochastic noise assumptions, this approach assumes:
 *   ||w_k||_∞ ≤ W_max,  ||v_k||_∞ ≤ V_max
 *
 * The estimate is a set (ellipsoid/bounding box) guaranteed to contain
 * the true state.
 *
 * When a measurement is lost, the uncertainty set expands according to
 * the process dynamics. This can be used to compute guaranteed
 * stability regions.
 *
 * Reference: Schweppe (1968), "Recursive State Estimation: Unknown but
 * Bounded Errors and System Inputs"
 *            Bertsekas & Rhodes (1971), "Recursive State Estimation
 * for a Set-Membership Description of Uncertainty"
 */
typedef struct {
    double* center;         /* Center of bounding ellipsoid, size n */
    double* shape;          /* Shape matrix (inverse of ellipsoid), n×n */
    double volume;          /* Current ellipsoid volume */
    double* temp_nn;        /* n×n temporary */

    /* System reference */
    const double* A;
    const double* C;
    int n;
    int p;
    double W_max;           /* Process disturbance bound */
    double V_max;           /* Measurement noise bound */
    double rho;             /* Forgetting factor for old data */
} SetMembershipEstimator;

/* --- Mode-Dependent Kalman Filter --- */

/**
 * Kalman filter for Markov-modulated packet loss (Gilbert-Elliott).
 *
 * The optimal estimator for Markovian loss is a bank of Kalman filters
 * — one per channel state — with mode probabilities. This is the
 * Interacting Multiple Model (IMM) approach.
 *
 * Reference: Blom & Bar-Shalom (1988), "The Interacting Multiple Model
 * Algorithm for Systems with Markovian Switching Coefficients"
 */
typedef struct {
    int n_modes;                    /* Number of channel modes/states */
    KalmanFilter** mode_filters;    /* One Kalman filter per mode */
    double* mode_probabilities;     /* P(channel_state = i | data) */
    double** mode_transitions;      /* Modulator transition matrix */

    /* Combined estimate (probability-weighted average) */
    double* combined_estimate;      /* n-vector */
    double* combined_covariance;    /* n×n */

    int n;                          /* State dimension */
    int p;                          /* Measurement dimension */
    int current_mode;               /* Current most likely mode */
} ModeDependentKalmanFilter;

/* ============================================================================
 * Core Kalman Filter API
 * ============================================================================ */

/**
 * Create a standard Kalman filter.
 * A, C, Q, R are row-major matrices.
 * x0 and P0 are initial estimate and covariance.
 */
KalmanFilter* pl_kf_create(const double* A, const double* C,
                            const double* Q, const double* R,
                            int n, int p,
                            const double* x0, const double* P0);
void pl_kf_free(KalmanFilter* kf);

/**
 * Standard Kalman predict step: x̂_{k|k-1} = A x̂_{k-1|k-1}
 *                                P_{k|k-1} = A P_{k-1|k-1} A' + Q
 */
void pl_kf_predict(KalmanFilter* kf);

/**
 * Standard Kalman update step with measurement y:
 *   K = P_pred C' (C P_pred C' + R)^{-1}
 *   x̂ = x̂_pred + K (y - C x̂_pred)
 *   P = (I - K C) P_pred
 */
void pl_kf_update(KalmanFilter* kf, const double* y);

/**
 * Combined predict + update step. Returns updated x̂.
 */
const double* pl_kf_step(KalmanFilter* kf, const double* y);

/** Get current state estimate (x̂_{k|k}). */
const double* pl_kf_get_estimate(const KalmanFilter* kf);

/** Get current error covariance P_{k|k}. */
const double* pl_kf_get_covariance(const KalmanFilter* kf);

/** Get Kalman gain K. */
const double* pl_kf_get_gain(const KalmanFilter* kf);

/** Compute trace of P as scalar performance metric. */
double pl_kf_trace_P(const KalmanFilter* kf);

void pl_kf_reset(KalmanFilter* kf, const double* x0, const double* P0);
void pl_kf_print(const KalmanFilter* kf);

/* ============================================================================
 * Intermittent Kalman Filter API
 * ============================================================================ */

/**
 * Create Kalman filter with intermittent observations.
 */
IntermittentKalmanFilter* pl_ikf_create(const double* A, const double* C,
                                          const double* Q, const double* R,
                                          int n, int p,
                                          const double* x0, const double* P0,
                                          double arrival_prob);
void pl_ikf_free(IntermittentKalmanFilter* ikf);

/**
 * Predict step (same as standard KF). Always executed.
 */
void pl_ikf_predict(IntermittentKalmanFilter* ikf);

/**
 * Conditional update step:
 *   If arrived (γ_k = 1): standard Kalman update with y
 *   If lost (γ_k = 0):   skip update, propagate prediction
 */
void pl_ikf_update(IntermittentKalmanFilter* ikf, const double* y, bool arrived);

/**
 * Combined predict + conditional update.
 */
const double* pl_ikf_step(IntermittentKalmanFilter* ikf, const double* y, bool arrived);

/**
 * Compute the expected error covariance E[P_k] using the modified
 * Riccati equation (Sinopoli 2004, Eq. 9):
 *
 * E[P_{k+1}] = A E[P_k] A' + Q
 *              - γ · E[ A P_k C' (C P_k C' + R)^{-1} C P_k A' ]
 *
 * Note: The expectation of the update term is approximated by
 * evaluating it at E[P_k] (a form of certainty equivalence).
 *
 * Complexity: O(n^3 + n²p + np² + p³).
 */
void pl_ikf_expected_covariance(IntermittentKalmanFilter* ikf);

/**
 * Estimate the critical arrival probability γ_c.
 *
 * Lower bound (Sinopoli 2004):
 *   γ_c ≥ 1 - 1/ρ(A)²
 *
 * Upper bound (Plarre & Bullo 2007):
 *   γ_c ≤ 1 - 1/(max_{unstable modes} |λ_u(A)|²)
 *
 * This function computes both bounds given the system matrix A.
 */
void pl_ikf_critical_probability(IntermittentKalmanFilter* ikf);

/**
 * Check if the filter is mean-stable under current arrival probability.
 * Returns true if E[||P_k||] appears bounded (diagnostic check).
 */
bool pl_ikf_is_stable(const IntermittentKalmanFilter* ikf);

/**
 * Get trace of expected covariance as scalar stability metric.
 * Bounded trace → stable estimation.
 */
double pl_ikf_trace_expected(const IntermittentKalmanFilter* ikf);

const double* pl_ikf_get_estimate(const IntermittentKalmanFilter* ikf);
void pl_ikf_reset(IntermittentKalmanFilter* ikf, const double* x0, const double* P0);
void pl_ikf_print(const IntermittentKalmanFilter* ikf);

/* ============================================================================
 * Set-Membership Estimator API
 * ============================================================================ */

/**
 * Create bounded-error set-membership estimator.
 * @param A, C: System matrices (references)
 * @param n, p: Dimensions
 * @param center0: Initial estimate center
 * @param shape0: Initial ellipsoid shape (inverse covariance scaling)
 * @param W_max, V_max: Noise bounds (infinity norm)
 */
SetMembershipEstimator* pl_sme_create(const double* A, const double* C,
                                        int n, int p,
                                        const double* center0,
                                        const double* shape0,
                                        double W_max, double V_max);
void pl_sme_free(SetMembershipEstimator* sme);

/**
 * Time update: ellipsoid expands according to A and W_max.
 * F_{k+1} = A F_k A' + (1/β)Q  (for some β ∈ (0,1) chosen to minimize volume)
 */
void pl_sme_predict(SetMembershipEstimator* sme);

/**
 * Measurement update (when measurement Available):
 * Intersection of predicted ellipsoid with measurement strip:
 *   {x : |y_k - C x| ≤ V_max}
 *
 * Uses the ellipsoid intersection formula (Schweppe 1968):
 *   F_{k+1|k+1} = F_{k+1|k} + ρ C' C
 *   Center updated via weighted least squares.
 */
void pl_sme_update_available(SetMembershipEstimator* sme, const double* y);

/**
 * Measurement update (when measurement LOST):
 * No measurement information — the ellipsoid simply expands
 * with an additional term to account for the missing information.
 * The volume grows by factor dependent on A and W_max.
 */
void pl_sme_update_lost(SetMembershipEstimator* sme);

/**
 * Combined step with arrival flag.
 */
void pl_sme_step(SetMembershipEstimator* sme, const double* y, bool arrived);

const double* pl_sme_get_center(const SetMembershipEstimator* sme);
double pl_sme_get_volume(const SetMembershipEstimator* sme);
void pl_sme_print(const SetMembershipEstimator* sme);

/* ============================================================================
 * Mode-Dependent Kalman Filter API
 * ============================================================================ */

/**
 * Create an IMM-style estimator for Markovian packet loss.
 * @param n_modes: Number of channel states
 * @param transitions: K×K Markov transition matrix (row-major, row-stochastic)
 * @param arrival_rates: Probability of measurement arrival in each mode
 */
ModeDependentKalmanFilter* pl_mdkf_create(const double* A, const double* C,
                                            const double* Q, const double* R,
                                            int n, int p, int n_modes,
                                            const double* transitions,
                                            const double* arrival_rates,
                                            const double* x0, const double* P0);
void pl_mdkf_free(ModeDependentKalmanFilter* mdkf);

/**
 * IMM-style step: mixing → mode-matched filtering → mode probability update → combination.
 */
void pl_mdkf_step(ModeDependentKalmanFilter* mdkf,
                  const double* y, bool arrived, int channel_state);

const double* pl_mdkf_get_estimate(const ModeDependentKalmanFilter* mdkf);
int pl_mdkf_get_mode(const ModeDependentKalmanFilter* mdkf);
void pl_mdkf_print(const ModeDependentKalmanFilter* mdkf);

#endif /* PACKET_LOSS_ESTIMATOR_H */