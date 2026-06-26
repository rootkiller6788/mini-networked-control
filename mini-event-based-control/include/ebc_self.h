#ifndef EBC_SELF_H
#define EBC_SELF_H
#include "ebc_core.h"

/*
 * ebc_self.h -- Self-Triggered Control (STC) (L5: Algorithms, L8: Advanced)
 *
 * Self-triggered control computes the next control update time
 * based on the current state and a model of the plant dynamics.
 * Unlike ETC (which continuously monitors the state), STC
 * predicts when the state will deviate beyond threshold.
 *
 * The STC law (Mazo et al. 2010):
 *   1. At time t_k, measure state x(t_k)
 *   2. Compute next update time:
 *        t_{k+1} = t_k + Gamma(x(t_k))
 *      where Gamma(x) = max { tau : |x(tau) - x| <= sigma*|x| + eps }
 *   3. Schedule next update at t_{k+1}
 *   4. At t_{k+1}, go to step 1
 *
 * Key advantage: No continuous monitoring needed.
 * Key challenge: Requires accurate plant model.
 *
 * References:
 *   Velasco et al. (2003) -- Early self-triggered concept
 *   Mazo, Anta & Tabuada (2010) -- ISS self-triggered implementation
 *     "Self-triggered control for nonlinear systems", IEEE TAC 55(9)
 *   Anta & Tabuada (2010) -- Self-triggered stabilization
 *   Heemels et al. (2012) -- Survey, Section VI
 */

/* ---- Self-triggered next-time computation (L5) ---- */

/**
 * Compute next event time for linear STC.
 *
 * For dx/dt = Ax + Bu with linear feedback u = Kx(t_k):
 *   x(t) = exp(A*(t-t_k)) * x_k + G(t-t_k) * B*K * x_k
 *   where G(tau) = integral_0^tau exp(A*s) ds
 *
 * The next time is found as the max tau satisfying:
 *   |x(t_k + tau) - x_k| <= sigma * |x_k| + epsilon
 *
 * Uses bisection search for efficiency.
 *
 * @param A        System matrix (n x n)
 * @param B        Input matrix (n x m)
 * @param K        Feedback gain (m x n)
 * @param x_k      Current state at event time
 * @param n, m     Dimensions
 * @param sigma    Relative threshold
 * @param epsilon  Absolute tolerance
 * @param tau_max  Maximum allowed inter-event time
 * @param tol      Bisection tolerance
 * @return         Next inter-event time tau in (0, tau_max]
 */
double ebc_self_next_time_linear(const double* A, const double* B,
                                  const double* K,
                                  const double* x_k,
                                  int n, int m,
                                  double sigma, double epsilon,
                                  double tau_max, double tol);

/**
 * Matrix exponential computation via scaling-and-squaring.
 * exp(A*t) computed using the Higham (2005) algorithm.
 * Pade approximation with scaling for numerical stability.
 *
 * @param A     Input matrix (n x n)
 * @param n     Dimension
 * @param t     Time scaling factor
 * @param E     Output: exp(A*t) (n x n, caller-allocated)
 */
void ebc_matrix_exponential(const double* A, int n, double t, double* E);

/**
 * Compute integral of matrix exponential: G(tau) = int_0^tau exp(A*s) ds.
 * Uses the augmented matrix approach:
 *   exp([[A, I]; [0, 0]] * tau) = [[exp(A*tau), G(tau)]; [0, I]]
 */
void ebc_matrix_exp_integral(const double* A, int n, double tau, double* G);

/* ---- Self-triggered for nonlinear systems (L8) ---- */

/**
 * Compute next event time for nonlinear STC using Lipschitz bounds.
 *
 * For dx/dt = f(x, k(x_k)), if f is Lipschitz with constant L:
 *   |x(t) - x_k| <= (|f(x_k, k(x_k))| / L) * (exp(L*tau) - 1)
 *
 * The next time satisfies:
 *   (|f_k| / L) * (exp(L*tau) - 1) <= sigma * |x_k| + epsilon
 *
 * @param f            Dynamics function
 * @param ctx          Dynamics context
 * @param K            Feedback gain
 * @param m            Input dimension
 * @param x_k          Current state
 * @param n            State dimension
 * @param L            Lipschitz constant estimate
 * @param sigma, eps   Threshold parameters
 * @param tau_max      Maximum allowed IET
 * @return             Next IET
 */
double ebc_self_next_time_nonlinear(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx,
    const double* K, int m,
    const double* x_k, int n,
    double L, double sigma, double epsilon, double tau_max);

/**
 * Run a complete self-triggered simulation.
 *
 * Unlike ETC (continuous measurement), STC only evaluates the
 * plant at event times. Between events, the plant evolves
 * autonomously with constant control input.
 *
 * @param sys       System (with dynamics set)
 * @param ctrl      Controller
 * @param T         Simulation horizon
 * @param dt        Integration step size
 * @param tp        Trigger parameters
 * @param traj      Output: state trajectory (allocated internally)
 * @param traj_len  Output: trajectory length
 * @param events    Output: event times (allocated internally)
 * @param evt_len   Output: number of events
 * @return          0 on success
 */
int ebc_self_simulate(EBC_System* sys, const EBC_Controller* ctrl,
                       double T, double dt,
                       const EBC_TriggerParams* tp,
                       double** traj, int* traj_len,
                       double** events, int* evt_len);

/* ---- Comparison: STC vs ETC vs Periodic ---- */

/**
 * Compare STC, ETC and periodic control on the same system.
 * Outputs performance metrics for all three schemes.
 */
typedef struct {
    EBC_Performance stc;
    EBC_Performance etc;
    EBC_Performance periodic;
} EBC_ComparisonResult;

EBC_ComparisonResult ebc_compare_schemes(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx, int n, int m,
    const double* A, const double* K,
    const double* x0, double T, double dt,
    double sigma, double epsilon, double period_h);

#endif /* EBC_SELF_H */
