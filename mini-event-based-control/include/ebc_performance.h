#ifndef EBC_PERFORMANCE_H
#define EBC_PERFORMANCE_H
#include "ebc_core.h"
#include "ebc_self.h"

/*
 * ebc_performance.h -- Performance Analysis of Event-Based Control (L6, L7)
 *
 * Compares event-triggered, self-triggered, periodic event-triggered,
 * and traditional periodic control in terms of:
 *   - Communication reduction
 *   - Control performance (ISE, IAE, ITAE)
 *   - Energy efficiency
 *   - Robustness to disturbances and delays
 *
 * Key metric: Communication reduction ratio
 *   R = 1 - N_events / N_periodic
 * where N_periodic = T / h for a baseline periodic controller.
 *
 * References:
 *   Astrom & Bernhardsson (1999) -- First comparison paper
 *   Heemels et al. (2008) -- ETC vs periodic for linear systems
 *   Lunze & Lehmann (2010) -- State-feedback event-based control
 *   Trimpe & D'Andrea (2014) -- Laboratory comparison ETC vs periodic
 */

/* ---- Performance computation ---- */

/**
 * Compute performance from trajectory data.
 *
 * Takes the state trajectory and event times and computes
 * all standard performance metrics in one pass.
 *
 * Complexity: O(traj_len * n)
 */
EBC_Performance ebc_compute_performance(
    const double* traj, int traj_len, int n,
    const double* events, int evt_len,
    double dt, double T,
    const double* K, int m,
    double comm_energy_per_event,
    double compute_energy_per_step);

/**
 * Compute Integral of Squared Error: ISE = sum |x_i|^2 * dt
 */
double ebc_ise(const double* traj, int len, int n, double dt);

/**
 * Compute Integral of Absolute Error: IAE = sum |x_i| * dt
 */
double ebc_iae(const double* traj, int len, int n, double dt);

/**
 * Compute Integral of Time-weighted Absolute Error:
 *   ITAE = sum t_i * |x_i| * dt
 */
double ebc_itae(const double* traj, int len, int n, double dt);

/**
 * Compute Integral of Squared Control Input:
 *   ISCI = sum |u_i|^2 * dt
 */
double ebc_isci(const double* K, const double* traj,
                int len, int n, int m, double dt);

/**
 * Compute settling time (time to stay within 2% of final value).
 */
double ebc_settling_time(const double* traj, int len, int n, double dt,
                          double threshold);

/**
 * Compute maximum overshoot relative to steady-state value.
 */
double ebc_overshoot(const double* traj, int len, int n, double dt);

/**
 * Compute average inter-event time.
 */
double ebc_average_iet(const double* events, int evt_len);

/**
 * Compute maximum state deviation between events.
 * Scans trajectory to find max |x(t) - x(t_k)| where t_k is the
 * last event time before t.
 */
double ebc_max_inter_event_deviation(const double* traj, int traj_len,
                                      int n,
                                      const double* events, int evt_len,
                                      double dt);

/* ---- Multi-scenario comparison (L7: Applications) ---- */

/**
 * Run a comprehensive comparison between all triggering schemes.
 *
 * Tests the same system under:
 *   1. Continuous ETC with mixed threshold
 *   2. Self-triggered with model-based prediction
 *   3. Periodic ETC with fixed sampling period
 *   4. Traditional periodic with equivalent rate
 *
 * @return Comparison result with metrics for all schemes
 */
EBC_ComparisonResult ebc_compare_all_schemes(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx, int n, int m,
    const double* A, const double* B, const double* K,
    const double* x0, double T, double dt,
    double sigma, double epsilon, double period_h);

/**
 * Benchmark communication reduction vs performance degradation.
 *
 * Sweeps sigma from 0.01 to 0.99, computing:
 *   - Comm reduction ratio
 *   - ISE degradation vs ideal continuous control
 *
 * Outputs a Pareto frontier of the communication-performance trade-off.
 */
typedef struct {
    int    n_points;
    double* sigma_values;
    double* comm_reduction;
    double* ise_degradation;
    double* iae_degradation;
    double* min_iet;
} EBC_ParetoFrontier;

EBC_ParetoFrontier ebc_pareto_frontier(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx, int n, int m,
    const double* A, const double* B, const double* K,
    const double* x0, double T, double dt,
    double epsilon, int n_sigma_points);

/** Free the Pareto frontier data */
void ebc_pareto_free(EBC_ParetoFrontier* pf);

/**
 * Robustness analysis: test performance under additive disturbance.
 *
 * Disturbance model: dx/dt = f(x, u) + w(t)
 * where |w(t)| <= w_bound.
 *
 * Tests whether event-triggered control maintains stability
 * and acceptable performance under bounded disturbances.
 */
typedef struct {
    double w_bound;         /* disturbance bound */
    double ise_mean;        /* mean ISE across noise realizations */
    double ise_std;         /* std dev of ISE */
    double comm_reduction;  /* average communication reduction */
    bool   stable;          /* true if stable for all realizations */
    int    zeno_count;      /* number of realizations with Zeno */
    int    n_trials;        /* number of Monte Carlo trials */
} EBC_RobustnessResult;

EBC_RobustnessResult ebc_robustness_analysis(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx, int n, int m,
    const double* K, const double* x0,
    double T, double dt, double sigma, double epsilon,
    double w_bound, int n_trials);

#endif /* EBC_PERFORMANCE_H */
