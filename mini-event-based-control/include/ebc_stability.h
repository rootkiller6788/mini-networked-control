#ifndef EBC_STABILITY_H
#define EBC_STABILITY_H
#include "ebc_core.h"

/*
 * ebc_stability.h -- Stability Analysis for Event-Based Systems (L4)
 *
 * Analyzes stability of event-triggered control systems using
 * ISS-Lyapunov theory. Key results:
 *
 * Theorem (Tabuada 2007, Theorem III.1):
 *   If there exists an ISS-Lyapunov function V such that:
 *     alpha1(|x|) <= V(x) <= alpha2(|x|)
 *     grad V * f(x, k(x+e)) <= -alpha(|x|) + gamma(|e|)
 *   and the event condition enforces
 *     gamma(|e|) <= sigma * alpha(|x|)  for some sigma in (0,1)
 *   then the closed-loop system is asymptotically stable.
 *
 * Theorem (Heemels et al. 2012, Theorem V.1 - Linear case):
 *   For dx/dt = Ax + BK(x+e), if A+BK is Hurwitz, then for any
 *   sigma in (0,1) there exists P > 0 solving
 *     (A+BK)'P + P(A+BK) = -Q
 *   such that |e| <= sigma*|x| implies dV/dt < 0.
 *   The minimum inter-event time is lower bounded by
 *     tau_min = sigma / (||A|| + sigma*(||A|| + ||BK||)) * ln(...)
 *
 * References:
 *   Tabuada (2007), IEEE TAC 52(9)
 *   Heemels, Johansson & Tabuada (2012), IEEE TAC 57(3)
 *   Sontag (2008) -- ISS theory
 *   Khalil (2002) -- Nonlinear Systems, Ch. 4
 */

/* ---- Lyapunov equation for ETC ---- */

/**
 * Solve the Lyapunov equation for event-triggered linear systems.
 *
 * Given A_cl = A + B*K (Hurwitz), solve A_cl'*P + P*A_cl = -Q
 * for P > 0. Uses the Bartels-Stewart algorithm via Schur decomposition.
 *
 * @param A       System matrix (n x n) [input]
 * @param B       Input matrix (n x m) [input]
 * @param K       Feedback gain (m x n) [input]
 * @param n       State dimension
 * @param m       Input dimension
 * @param Q       Right-hand side matrix (n x n), usually I_n [input]
 * @param P_out   Solution matrix (n x n) [output, caller-allocated]
 * @return        0 on success, -1 if A+BK is not Hurwitz
 *
 * Reference: Bartels & Stewart (1972), Comm. ACM 15(9): 820-826
 */
int ebc_lyapunov_solve(const double* A, const double* B, const double* K,
                        int n, int m, const double* Q, double* P_out);

/**
 * Compute the ISS gain gamma for linear event-triggered systems.
 *
 * For dx/dt = Ax + BK(x+e), the ISS gain is:
 *   gamma = 2 * ||P*B*K||
 * where P solves the Lyapunov equation.
 *
 * Complexity: O(n^3) for Lyapunov solve + O(n^2) for norm computation.
 */
double ebc_iss_gain_linear(const double* A, const double* B,
                            const double* K, int n, int m);

/**
 * Compute the critical sigma for event-triggered stability.
 *
 * sigma_crit = lambda_min(Q) / (2 * ||P*B*K||)
 *
 * Any sigma < sigma_crit guarantees asymptotic stability
 * under the event condition |e| <= sigma * |x|.
 *
 * Reference: Tabuada (2007), Eq. (14)
 */
double ebc_critical_sigma(const double* A, const double* B,
                           const double* K, int n, int m);

/**
 * Compute the guaranteed minimum inter-event time.
 *
 * For linear systems with mixed threshold:
 *   tau_min = (1 / ||A_cl||) * ln(1 + sigma * ||A_cl|| * xi)
 *   where xi = epsilon / (sigma * |x_k| + ... )
 *
 * This provides a positive lower bound, ruling out Zeno behavior.
 */
double ebc_minimum_iet_linear(const double* A, const double* B,
                               const double* K, int n, int m,
                               double sigma, double epsilon);

/* ---- Stability certificate generation ---- */

/**
 * Generate a full stability certificate for linear ETC systems.
 *
 * Computes: Lyapunov matrix P, alpha1/alpha2 constants,
 * ISS gain gamma, critical sigma, minimum IET, convergence rate.
 *
 * @param A, B, K, n, m  System matrices
 * @param sigma          Event threshold parameter
 * @param epsilon        Absolute tolerance
 * @param cert           Output certificate [caller-allocated]
 * @return               0 on success
 */
int ebc_stability_certify_linear(const double* A, const double* B,
                                  const double* K, int n, int m,
                                  double sigma, double epsilon,
                                  EBC_StabilityCert* cert);

/** Free internal allocations in a stability certificate */
void ebc_stability_cert_free(EBC_StabilityCert* cert);

/* ---- ISS-Lyapunov verification for nonlinear ETC ---- */

/**
 * Verify ISS property for nonlinear event-triggered systems.
 *
 * Given a candidate Lyapunov function V(x) and its Lie derivative
 * along the closed-loop dynamics, check the ISS condition:
 *   dV/dt <= -alpha3 * V(x) + gamma * |e|^2
 *
 * Uses numerical sampling over a grid in state and error space.
 *
 * @param V_func        Lyapunov function V(x) [input: returns V(x)]
 * @param dV_func       Lie derivative dV/dt(x,e) [input]
 * @param n             State dimension
 * @param bounds        State space bounds: [x_min, x_max] per dimension
 * @param grid_points   Number of grid points per dimension
 * @param alpha3_out    Output: decay constant estimate
 * @param gamma_out     Output: ISS gain estimate
 * @return              true if ISS condition verified numerically
 */
bool ebc_verify_iss_nonlinear(
    double (*V_func)(const double* x, int n, void* ctx),
    double (*dV_func)(const double* x, const double* e, int n, void* ctx),
    void* ctx, int n,
    const double* bounds_lo, const double* bounds_hi,
    int grid_points,
    double* alpha3_out, double* gamma_out);

/* ---- Zeno behavior detection ---- */

/**
 * Prove absence of Zeno behavior for given event-triggering scheme.
 *
 * Zeno behavior: infinite number of events in finite time.
 * For mixed threshold with epsilon > 0, the inter-event time
 * is bounded below by tau_min > 0, ruling out Zeno.
 *
 * @param sys    Event-based system
 * @param tp     Trigger parameters
 * @param T      Time horizon for simulation
 * @param dt     Simulation step size
 * @return       true if no Zeno detected during simulation
 */
bool ebc_zeno_free_proof(const EBC_System* sys,
                          const EBC_TriggerParams* tp,
                          double T, double dt);

/**
 * Compute the theoretical minimum inter-event time bound
 * from the ISS-Lyapunov characterization.
 *
 * For mixed threshold: tau_min = epsilon / L
 * where L is the Lipschitz constant of the dynamics.
 */
double ebc_theoretical_iet_lower_bound(const EBC_System* sys,
                                        const EBC_TriggerParams* tp);

/**
 * Estimate the Lipschitz constant of the closed-loop dynamics
 * along trajectories. Used for IET lower bound computation.
 */
double ebc_lipschitz_estimate(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx, int n,
    const double* K, int m,
    const double* bounds_lo, const double* bounds_hi,
    int samples);

#endif /* EBC_STABILITY_H */
