/**
 * blc_control.h — Bandwidth-Limited Controller Design and Analysis
 *
 * Designs controllers that explicitly account for communication constraints:
 * quantization, bit-rate limits, packet loss, and latency.
 *
 * Core design methodologies:
 *
 * 1. Quantized State Feedback (Delchamps, 1990):
 *    u[k] = K·Q(x[k]) where Q is a finite-level quantizer.
 *    The controller sees only the quantized state.
 *
 * 2. Certainty-Equivalence with Quantization (Tatikonda, 2000):
 *    Separate estimation and control with quantized measurements,
 *    accounting for quantization in the Kalman gain computation.
 *
 * 3. Robust Control for Quantization Errors (Fu & Xie, 2005):
 *    Treat quantization error as bounded disturbance using
 *    H∞ or sector-bound methods.
 *
 * 4. Model Predictive Control with Rate Constraints:
 *    Include bit-rate budget as a constraint in the MPC optimization.
 *
 * 5. Minimum Attention Control (Brockett, 1997):
 *    Design controllers that require minimal information from sensors,
 *    exploiting open-loop dynamics when possible.
 *
 * Knowledge coverage: L5 (Control Design), L6 (Canonical Problems)
 */

#ifndef BLC_CONTROL_H
#define BLC_CONTROL_H

#include "blc_core.h"
#include "blc_datarate.h"
#include "blc_encoding.h"

/* ================================================================
 * LQR under Quantization
 *
 * Standard infinite-horizon LQR:
 *   J = ∫₀^∞ (x'Qx + u'Ru) dt
 *   u = -Kx, K = R^{-1}B'P, where P solves ARE: A'P + PA - PBR^{-1}B'P + Q = 0
 *
 * Under quantization, the actual control is:
 *   u = -K·Q(x) = -K·(x + ε), where ε = quantization error
 *
 * The quantization error acts as an additional disturbance:
 *   ẋ = (A - BK)x - BK·ε
 *
 * The maximum tolerable quantization error depends on the stability
 * margin of the closed-loop system.
 * ================================================================ */

/** Bandwidth-limited LQR controller state */
typedef struct {
    int     n_states;
    int     n_inputs;
    double  Q_diag[BLC_MAX_STATES];     /** State cost (diagonal) */
    double  R_diag[BLC_MAX_STATES];     /** Input cost (diagonal) */
    double  K[BLC_MAX_STATES][BLC_MAX_STATES]; /** LQR gain matrix */
    double  P[BLC_MAX_STATES][BLC_MAX_STATES]; /** Riccati solution */
    double  Ac[BLC_MAX_STATES][BLC_MAX_STATES];/** Closed-loop A - BK */
    double  cost;                       /** Last computed cost */
    double  cost_with_quantization;     /** Cost including quantization effect */
    double  quantization_penalty;       /** Additional cost due to quantization */
    int     are_iterations;             /** Newton iterations for ARE */
    bool    converged;                  /** ARE convergence flag */
} BLCLQRController;

/** Minimum attention controller (Brockett, 1997).
 *  Minimizes the information rate from sensor to controller
 *  while maintaining a guaranteed level of performance.
 *
 *  Concept: the controller "looks" at sensor data only when
 *  necessary — otherwise, it predicts using the plant model.
 *
 *  When the prediction error exceeds a threshold, the controller
 *  requests a measurement update.
 */
typedef struct {
    double  attention_rate;             /** Current attention rate (Hz) */
    double  min_attention_rate;         /** Theoretical minimum attention */
    double  performance_bound;          /** Guaranteed performance γ */
    double  error_threshold;            /** Error threshold for attention */
    int     look_interval;              /** Steps between mandatory looks */
    int     steps_since_last_look;      /** Current steps without measurement */
    double  prediction_error;           /** Current prediction error norm */
    double  total_looks;                /** Cumulative sensor accesses */
    double  average_attention;          /** Running average attention rate */
} BLCMinAttentionController;

/** Robust controller treating quantization as sector-bounded nonlinearity.
 *  Sector bound: the quantization error ε satisfies
 *    ε'(ε - δx) ≤ 0   for some δ ∈ (0, 1)
 *  where δ = 1/(2N) for uniform quantization with N levels per unit range.
 *
 *  Using the Circle Criterion or Popov criterion, we can find
 *  the maximum tolerable quantization coarseness.
 */
typedef struct {
    double  sector_delta;               /** Sector bound parameter δ */
    double  max_tolerable_delta;        /** Maximum δ for stability */
    bool    is_absolutely_stable;       /** Absolute stability verdict */
    double  gain_margin;                /** Gain margin under quantization */
    double  phase_margin;               /** Phase margin under quantization */
    double  h_inf_norm;                 /** H∞ norm of quantization→output */
} BLCRobustQuantCtrl;

/** Bandwidth-limited Model Predictive Control
 *
 *  At each step k, solve:
 *    min_{u_{k},...,u_{k+N-1}} Σ_{i=0}^{N-1} (||x_{k+i}||_Q + ||u_{k+i}||_R)
 *    s.t.  x_{k+i+1} = f(x_{k+i}, u_{k+i})
 *          |u_{k+i}| ≤ u_max
 *          R(u) ≤ R_budget   (bit-rate constraint)
 *
 *  The rate constraint is nonlinear and non-convex.
 *  We use a relaxation: limit the number of non-zero control
 *  inputs as a proxy for bit-rate.
 */
typedef struct {
    int     horizon;                    /** Prediction horizon N */
    int     max_nonzero_inputs;         /** Proxy for bit-rate constraint */
    double  u_opt[BLC_MAX_STATES*16];   /** Optimal control sequence */
    double  x_pred[BLC_MAX_STATES*16];  /** Predicted state trajectory */
    double  cost_opt;                   /** Optimal cost */
    int     active_constraints;         /** Number of active constraints */
    int     qp_iterations;              /** Quadratic program iterations */
} BLCMPCController;

/* ================================================================
 * LQR Controller API
 * ================================================================ */

/** Initialize bandwidth-limited LQR controller.
 *  @param lqr Controller structure
 *  @param n_states State dimension
 *  @param n_inputs Input dimension
 *  @param Q_diag Diagonal of state cost matrix
 *  @param R_diag Diagonal of input cost matrix
 *  @return 0 on success */
int     blc_lqr_init(BLCLQRController* lqr, int n_states, int n_inputs,
                      const double* Q_diag, const double* R_diag);

/** Solve continuous-time Algebraic Riccati Equation (CARE):
 *    A'P + PA - PBR^{-1}B'P + Q = 0
 *  Using Newton-Kleinman iteration.
 *  @param lqr Controller state
 *  @param A System matrix (n×n, row-major)
 *  @param B Input matrix (n×m, row-major)
 *  @param max_iter Maximum Newton iterations
 *  @param tol Convergence tolerance
 *  @return Number of iterations */
int     blc_lqr_solve(BLCLQRController* lqr, const double* A,
                       const double* B, int max_iter, double tol);

/** Compute the closed-loop eigenvalues of A - BK */
int     blc_lqr_closed_loop_poles(const BLCLQRController* lqr,
                                   const double* A, const double* B,
                                   double* poles_real, double* poles_imag);

/** Compute LQR cost with quantization effect.
 *  J_quant = J_nominal + tr(S·Σ_ε)
 *  where Σ_ε = diag(Δ²/12) is the quantization error covariance
 *  and S solves a Lyapunov equation.
 *  @param lqr Controller
 *  @param A Open-loop system matrix
 *  @param B Input matrix
 *  @param quant_errors Per-state quantization error standard deviations
 *  @return Total cost */
double  blc_lqr_cost_with_quantization(const BLCLQRController* lqr,
                                        const double* A, const double* B,
                                        const double* quant_errors);

/** Compute the quantization penalty: J_quant - J_nominal */
double  blc_lqr_quantization_penalty(const BLCLQRController* lqr,
                                      const double* A, const double* B,
                                      const double* quant_errors);

/** Get LQR gain matrix */
void    blc_lqr_get_gain(const BLCLQRController* lqr, double* K_flat);

/** Get Riccati solution matrix */
void    blc_lqr_get_riccati(const BLCLQRController* lqr, double* P_flat);

/** Free LQR controller */
void    blc_lqr_free(BLCLQRController* lqr);

/* ================================================================
 * Quantized State Feedback API
 * ================================================================ */

/** Quantized state feedback control law.
 *  u = K · Q(x)
 *  @param sys Bandwidth-limited system
 *  @param lqr LQR controller
 *  @param x Current state
 *  @param u Output control vector
 *  @param quant_error Output quantization error vector (optional, can be NULL) */
void    blc_quantized_state_feedback(BLCSystem* sys,
                                      const BLCLQRController* lqr,
                                      const double* x, double* u,
                                      double* quant_error);

/** Compute maximum tolerable quantization step for stability.
 *  Based on the small-gain theorem: the quantization error is
 *  a bounded disturbance, and stability requires:
 *    ||Δ|| < (1 / ||(sI - A_c)^{-1}BK||_∞)
 *  @param A_c Closed-loop matrix A-BK
 *  @param B Input matrix
 *  @param K LQR gain
 *  @param n State dimension
 *  @return Maximum allowable quantization step Δ_max */
double  blc_max_quantization_step(const double* A_c, const double* B,
                                   const double* K, int n);

/** Compute required quantization levels for guaranteed stability.
 *  N_req = ceil(2·U / Δ_max), where U is the operating range.
 *  @param range_hi Upper bound of operating range
 *  @param range_lo Lower bound of operating range
 *  @param delta_max Maximum allowable step
 *  @return Minimum number of levels */
int     blc_required_quantization_levels(double range_hi, double range_lo,
                                          double delta_max);

/** Compute bit rate from quantization levels.
 *  R = log₂(N₁) + log₂(N₂) + ... + log₂(N_n)  bits/sample
 *  @param levels Array of quantization levels per state component
 *  @param n_states Number of states
 *  @return bits per sample */
double  blc_bitrate_from_levels(const int* levels, int n_states);

/* ================================================================
 * Minimum Attention Controller API
 * ================================================================ */

/** Initialize minimum attention controller */
int     blc_minattn_init(BLCMinAttentionController* mac,
                          double error_threshold, int look_interval);

/** Check if measurement is needed now */
bool    blc_minattn_should_look(BLCMinAttentionController* mac);

/** Update with new measurement */
void    blc_minattn_update(BLCMinAttentionController* mac,
                            double prediction_error);

/** Compute minimum attention rate from Data Rate Theorem */
double  blc_minattn_theoretical_rate(const BLCPlant* plant);

/** Get current attention statistics */
void    blc_minattn_stats(const BLCMinAttentionController* mac,
                           double* avg_rate, double* total_looks);

/* ================================================================
 * Robust Control under Quantization API
 * ================================================================ */

/** Check absolute stability under quantization using Circle Criterion.
 *  The quantizer satisfies the sector condition:
 *    φ(σ)' (φ(σ) - δσ) ≤ 0
 *  where φ(σ) = Q(σ) - σ is the quantization error function.
 *
 *  Circle Criterion: The closed-loop system is absolutely stable
 *  if the Nyquist plot of G(s) = K(sI - A + BK)^{-1}B
 *  lies strictly to the right of the vertical line through -1/δ.
 *
 *  @param rqc Controller structure (populated with system data)
 *  @param A_c Closed-loop matrix A-BK
 *  @param B Input matrix
 *  @param K LQR gain
 *  @param n State dimension
 *  @return 0 if stable, 1 if unstable */
int     blc_robust_circle_criterion(BLCRobustQuantCtrl* rqc,
                                     const double* A_c, const double* B,
                                     const double* K, int n);

/** Compute H∞ norm of quantization-to-output transfer function */
double  blc_robust_hinf_norm(const double* A_c, const double* B,
                              const double* C, int n);

/* ================================================================
 * Bandwidth-Limited MPC API
 * ================================================================ */

/** Initialize bandwidth-limited MPC controller */
int     blc_mpc_init(BLCMPCController* mpc, int horizon,
                      int max_nonzero_inputs);

/** Solve one MPC step.
 *  Uses a simple projected gradient method for the QP subproblem,
 *  with sparsity-promoting L1 regularization as proxy for rate constraint.
 *  @param mpc Controller state
 *  @param x0 Current state
 *  @param A System matrix
 *  @param B Input matrix
 *  @param Q_diag Diagonal state cost
 *  @param R_diag Diagonal input cost
 *  @param u_prev Previous control (for rate penalty)
 *  @param u_next Output next control
 *  @return 0 on success */
int     blc_mpc_solve(BLCMPCController* mpc, const double* x0,
                       const double* A, const double* B,
                       const double* Q_diag, const double* R_diag,
                       const double* u_prev, double* u_next);

/** Compute effective bit-rate of MPC solution.
 *  Rate = number of non-zero control components / horizon
 *  @param mpc MPC controller after solve
 *  @return Effective bit rate proxy */
double  blc_mpc_effective_rate(const BLCMPCController* mpc);

/** Free MPC controller resources */
void    blc_mpc_free(BLCMPCController* mpc);

#endif /* BLC_CONTROL_H */