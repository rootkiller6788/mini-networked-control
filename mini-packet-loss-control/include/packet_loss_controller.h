#ifndef PACKET_LOSS_CONTROLLER_H
#define PACKET_LOSS_CONTROLLER_H

#include "packet_loss_core.h"
#include <stdbool.h>

/* ============================================================================
 * Packet Loss Controller — Control Strategies Under Packet Dropout
 *
 * Implements control strategies for networked control systems (NCS)
 * subject to packet loss in the sensor-to-controller (S→C) and
 * controller-to-actuator (C→A) channels.
 *
 * Key references:
 *   - Schenato et al. (2007): "Foundations of Control and Estimation
 *     Over Lossy Networks" — TCP/UDP-like protocol distinction,
 *     LQG optimal controller structure.
 *   - Gupta, Hassibi, Murray (2007): "Optimal LQG Control Across
 *     Packet-Dropping Links" — separation principle breakdown.
 *   - Imer, Yuksel, Basar (2006): "Optimal Control of LTI Systems
 *     over Unreliable Communication Links" — TCP-like optimal control.
 * ============================================================================ */

/* --- LTI System Model --- */

/**
 * Linear Time-Invariant (LTI) discrete-time system:
 *   x_{k+1} = A·x_k + B·u_k + w_k     (process, w_k ~ N(0, Q))
 *   y_k     = C·x_k + v_k             (measurement, v_k ~ N(0, R))
 *   z_k     = C_z·x_k                 (controlled output)
 *
 * State dimension n, input dimension m, measurement dimension p.
 */
typedef struct {
    double* A;      /* System matrix, n×n, row-major */
    double* B;      /* Input matrix, n×m, row-major */
    double* C;      /* Output matrix, p×n, row-major */
    double* Cz;     /* Controlled output matrix, q×n, row-major */
    double* Q;      /* Process noise covariance, n×n */
    double* R;      /* Measurement noise covariance, p×p */
    int n;          /* State dimension */
    int m;          /* Input dimension */
    int p;          /* Output dimension */
    int q;          /* Controlled output dimension */
} LTISystem;

/* --- Controller State --- */

/**
 * TCP-like controller state (knows which packets arrived).
 *
 * Under TCP-like protocols, the controller receives acknowledgment
 * (ACK/NACK) for each packet sent to the actuator. This means the
 * controller knows exactly which control inputs were applied.
 *
 * Key property (Schenato 2007): Under TCP-like protocols with
 * packet acknowledgments, the separation principle holds.
 * The optimal controller is:
 *   u_k = -L_k · x̂_{k|k}
 * where L_k is the standard LQR gain and x̂ is the Kalman estimate.
 */
typedef struct {
    /* System model reference */
    const LTISystem* sys;

    /* LQR gain matrix L (m×n) — optimal feedback gain */
    double* L;

    /* Kalman filter estimate of state x̂ */
    double* x_hat;              /* Current state estimate, size n */
    double* P;                  /* Estimation error covariance, n×n */

    /* Prediction for next step */
    double* x_pred;             /* Predicted state, size n */
    double* P_pred;             /* Predicted covariance, n×n */

    /* Packet delivery tracking */
    bool* sensor_acked;         /* Did sensor packet arrive? */
    bool* actuator_acked;       /* Did actuator packet arrive? */

    /* Last applied control (for hold strategies) */
    double* u_last;             /* Last control signal, size m */
    double* u_sent;             /* Most recently sent control, size m */

    /* Loss counters and statistics */
    unsigned long sensor_losses;
    unsigned long actuator_losses;
    unsigned long total_sensor_attempts;
    unsigned long total_actuator_attempts;
    double sensor_loss_rate;
    double actuator_loss_rate;

    /* Controller parameters */
    TransportProtocol protocol;
    HoldStrategy sensor_hold;    /* Strategy for sensor→controller loss */
    HoldStrategy actuator_hold;  /* Strategy for controller→actuator loss */

    /* History buffers for analysis */
    double* cost_history;        /* Running cost J_k */
    int cost_history_len;
    int cost_history_cap;
} TCPController;

/**
 * UDP-like controller state (does NOT know which packets arrived).
 *
 * Without acknowledgments, the controller does not know whether its
 * control signal was actually applied to the plant. This creates
 * a coupling between estimation and control — the separation
 * principle does NOT hold.
 *
 * The optimal UDP-like control is a nonlinear function of the
 * information state (Gupta et al., 2007).
 */
typedef struct {
    /* System model reference */
    const LTISystem* sys;

    /* LQR gain — suboptimal under UDP-like (separation fails) */
    double* L;

    /* State estimate — blindly propagated */
    double* x_hat;
    double* P;

    /* Uncertainty about what was actually applied */
    double* u_last_sent;        /* Last control signal we sent, size m */
    double* u_possible;         /* Set of possible applied controls */
    int n_possible_outcomes;    /* Number of distinct delivery patterns */

    /* Belief state about actuator delivery */
    double prob_received;       /* Belief that last packet was received */

    /* Performance tracking */
    unsigned long sensor_losses;
    unsigned long actuator_attempts;
    double* cost_history;
    int cost_history_len;
    int cost_history_cap;
} UDPController;

/* --- LQR Design --- */

/**
 * Discrete-time LQR cost:
 *   J = Σ (x_k'·Q·x_k + u_k'·R·u_k)
 *
 * Solution: Riccati equation
 *   P = A'PA - A'PB(R + B'PB)^{-1}B'PA + Q
 *
 * Optimal gain: L = (R + B'PB)^{-1} B'PA (row vector per input)
 */
typedef struct {
    double* P;          /* Solution to DARE (Discrete Algebraic Riccati Equation) */
    double* L;          /* Optimal LQR gain matrix, m×n */
    double* K_cl;       /* Closed-loop matrix: A - B·L */
    int n;
    int m;
    int iterations;     /* Number of iterations to converge */
    double tolerance;   /* Convergence tolerance */
} LQRSolution;

/**
 * Solve the Discrete-time Algebraic Riccati Equation (DARE).
 *
 * Uses iterative value iteration:
 *   P_{k+1} = A'P_k A - A'P_k B (R + B'P_k B)^{-1} B'P_k A + Q, with P_0 = Q.
 *
 * Convergence is guaranteed if (A,B) is stabilizable and (A,sqrt(Q)) is detectable.
 *
 * Complexity: O(n^3 · iterations). Each iteration requires two n×n matrix multiplies
 * and one m×m matrix inversion.
 *
 * Reference: Bertsekas, "Dynamic Programming and Optimal Control" (2012), Vol. 1, §4.1.
 */
LQRSolution* pl_lqr_solve(const LTISystem* sys, double* R_ctrl,
                           int max_iter, double tol);
void pl_lqr_free(LQRSolution* sol);

/**
 * Compute closed-loop eigenvalues of A - B·L.
 * Returns the spectral radius (max |λ|) for stability check.
 */
double pl_lqr_spectral_radius(const LQRSolution* sol, const LTISystem* sys);

/**
 * Compute LQR cost for a trajectory.
 * J = Σ_{k=0}^{T-1} (x_k'·Q·x_k + u_k'·R·u_k) + x_T'·P·x_T
 */
double pl_lqr_cost(const LTISystem* sys, const LQRSolution* sol,
                   const double* x_traj, const double* u_traj,
                   int T);

/* --- TCP-like Controller API --- */

/**
 * Create a TCP-like controller for the given LTI system.
 *
 * The controller maintains knowledge of:
 *  - Full state estimate via Kalman filter (with intermittent measurements)
 *  - Which control packets were successfully delivered (via ACKs)
 *
 * Protocol: If sensor packet lost → hold strategy applied to measurement.
 *           If actuator ACK not received → hold strategy applied to control.
 */
TCPController* pl_tcp_controller_create(const LTISystem* sys,
                                         const LQRSolution* lqr);
void pl_tcp_controller_free(TCPController* ctrl);

/**
 * Set the hold strategies for sensor and actuator loss.
 * - sensor_hold: What to do when sensor→controller packet is lost
 * - actuator_hold: What to do when controller→actuator packet is lost
 */
void pl_tcp_controller_set_hold(TCPController* ctrl,
                                 HoldStrategy sensor_hold,
                                 HoldStrategy actuator_hold);

/**
 * Process a sensor measurement with possible packet loss.
 *
 * If measurement arrives (status = PACKET_RECEIVED):
 *   - Standard Kalman update: x̂ = x̂_pred + K(y - C·x̂_pred)
 *
 * If measurement lost (status ≠ PACKET_RECEIVED):
 *   - Apply sensor hold strategy:
 *     ZERO_INPUT:  x̂ = 0 (full reset, very conservative)
 *     ZERO_ORDER:  x̂ = x̂_pred (open-loop prediction)
 *     PREDICTIVE:  x̂ = open-loop + noise model adjustment
 *
 * Returns the updated state estimate x̂.
 */
const double* pl_tcp_controller_update_estimate(TCPController* ctrl,
                                                  const double* y,
                                                  PacketStatus sensor_status);

/**
 * Compute control signal considering possible actuator loss.
 *
 * The control signal is computed as: u = -L · x̂ (LQR feedback).
 *
 * If actuator_status ≠ PACKET_RECEIVED, the hold strategy is applied:
 *   ZERO_INPUT:  Return zero vector
 *   ZERO_ORDER:  Return u_last (previously applied control)
 *   PREDICTIVE:  Predicted optimal control based on state prediction
 *
 * Returns pointer to the control signal vector (size m).
 */
const double* pl_tcp_controller_compute_control(TCPController* ctrl,
                                                  PacketStatus actuator_status);

/**
 * Apply one complete control step:
 *   1. Update estimate from measurement (with possible sensor loss)
 *   2. Compute control (with possible actuator loss)
 *   3. Store results in history
 *
 * @param y: measurement vector (NULL if no measurement generated)
 * @param sensor_status: delivery status of sensor→controller packet
 * @param actuator_status: delivery status of controller→actuator packet
 * @param u_out: output buffer for applied control (size m)
 */
void pl_tcp_controller_step(TCPController* ctrl,
                             const double* y,
                             PacketStatus sensor_status,
                             PacketStatus actuator_status,
                             double* u_out);

void pl_tcp_controller_reset(TCPController* ctrl);
void pl_tcp_controller_print(const TCPController* ctrl);

/* --- UDP-like Controller API --- */

UDPController* pl_udp_controller_create(const LTISystem* sys,
                                         const LQRSolution* lqr);
void pl_udp_controller_free(UDPController* ctrl);

/**
 * UDP-like measurement update: same as TCP-like for estimation,
 * but maintains additional uncertainty about actuator delivery.
 */
const double* pl_udp_controller_update_estimate(UDPController* ctrl,
                                                  const double* y,
                                                  PacketStatus sensor_status);

/**
 * UDP-like control: controller does NOT know if packet reached actuator.
 * Uses belief-state approach: maintains probability distribution over
 * possible applied controls and chooses action to minimize expected cost.
 */
const double* pl_udp_controller_compute_control(UDPController* ctrl,
                                                  double loss_probability);

void pl_udp_controller_step(UDPController* ctrl,
                             const double* y,
                             PacketStatus sensor_status,
                             double actuator_loss_prob,
                             double* u_out);

void pl_udp_controller_reset(UDPController* ctrl);
void pl_udp_controller_print(const UDPController* ctrl);

/* --- LTI System Utilities --- */

/**
 * Create an LTI system model.
 * Allocates memory for all matrices (n, m, p dimensions).
 * Matrices are initialized to zero; caller must fill A, B, C, Q, R.
 */
LTISystem* pl_lti_create(int n, int m, int p, int q);
void pl_lti_free(LTISystem* sys);

/**
 * Compute spectral radius of A: ρ(A) = max|λ_i(A)|.
 * Uses power iteration — does NOT compute full eigendecomposition.
 *
 * Complexity: O(n^2 · iterations).
 */
double pl_lti_spectral_radius(const LTISystem* sys);

/**
 * Compute controllability matrix: C = [B AB A²B ... A^{n-1}B].
 * Returns the rank via Gaussian elimination (non-square matrix rank).
 *
 * The pair (A,B) is controllable iff rank(C) = n.
 */
double* pl_lti_controllability_matrix(const LTISystem* sys, int* rank_out);

/**
 * Compute observability matrix: O = [C; CA; CA²; ...; CA^{n-1}]'.
 * Returns the rank via Gaussian elimination.
 *
 * The pair (A,C) is observable iff rank(O) = n.
 */
double* pl_lti_observability_matrix(const LTISystem* sys, int* rank_out);

/**
 * Check if (A,B) is stabilizable: all uncontrollable eigenvalues are stable.
 * Stabilizability is weaker than controllability — unstable modes must
 * be controllable, stable modes may be uncontrollable.
 *
 * Returns true if system is stabilizable.
 */
bool pl_lti_is_stabilizable(const LTISystem* sys);

/**
 * Check if (A,C) is detectable: all unobservable eigenvalues are stable.
 * Dual of stabilizability — unstable modes must be observable.
 */
bool pl_lti_is_detectable(const LTISystem* sys);

/**
 * Set the system matrix A and compute its spectral radius.
 * A should be in row-major order: A[i*n + j] = A_{i,j}.
 */
void pl_lti_set_A(LTISystem* sys, const double* A_data, int n);

/** Set input matrix B (n×m, row-major). */
void pl_lti_set_B(LTISystem* sys, const double* B_data);

/** Set output matrix C (p×n, row-major). */
void pl_lti_set_C(LTISystem* sys, const double* C_data);

/** Set process noise covariance Q (n×n, row-major, must be symmetric PSD). */
void pl_lti_set_Q(LTISystem* sys, const double* Q_data);

/** Set measurement noise covariance R (p×p, row-major, symmetric PD). */
void pl_lti_set_R(LTISystem* sys, const double* R_data);

/**
 * Simulate one step of the LTI system open-loop.
 * x_next = A·x + B·u + w (where w ~ N(0,Q) if noise enabled).
 */
void pl_lti_step(const LTISystem* sys, const double* x,
                 const double* u, bool add_noise,
                 unsigned long* rng_state, double* x_next);

/**
 * Generate measurement: y = C·x + v (where v ~ N(0,R) if noise enabled).
 */
void pl_lti_measure(const LTISystem* sys, const double* x,
                    bool add_noise, unsigned long* rng_state, double* y);

/* --- Matrix Operations for Control --- */

/**
 * Matrix multiplication: C = A × B.
 * A: m×k, B: k×n, C: m×n (row-major).
 * Complexity: O(m·k·n).
 */
void pl_mat_mul(const double* A, const double* B, double* C,
                int m, int k, int n);

/**
 * Matrix transpose.
 * Complexity: O(rows·cols).
 */
void pl_mat_transpose(const double* A, double* AT, int rows, int cols);

/**
 * Matrix inversion via Gaussian elimination with partial pivoting.
 * In-place: A is overwritten with its inverse.
 * Returns true on success, false if singular (within ε tolerance).
 *
 * Complexity: O(n^3). Handles n ≤ 64 efficiently.
 *
 * Reference: Golub & Van Loan, "Matrix Computations" (2013), §3.2.
 */
bool pl_mat_invert(double* A, int n);

/**
 * Cholesky decomposition: A = L·L' where A is symmetric positive definite.
 * Lower-triangular L stored in the lower triangle of A.
 * Returns true on success, false if not SPD.
 *
 * Complexity: O(n^3/3). Numerically stable.
 */
bool pl_mat_cholesky(double* A, int n);

/**
 * Solve Lyapunov equation: A·X + X·A' + Q = 0.
 * Uses Bartels-Stewart algorithm for n ≤ 32.
 * X is output (n×n).
 *
 * Reference: Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C".
 */
bool pl_solve_lyapunov(const double* A, const double* Q,
                       double* X, int n);

#endif /* PACKET_LOSS_CONTROLLER_H */