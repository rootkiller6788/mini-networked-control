#include "ebc_core.h"
#include "ebc_performance.h"
#include "ebc_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
 * example3_performance_analysis.c -- Performance Analysis Demo
 *
 * Demonstrates performance metric computation and comparison
 * between event-triggered and periodic control.
 *
 * Also shows the communication-performance trade-off
 * via sigma parameter sweep (Pareto frontier).
 *
 * L6: Canonical problem -- performance benchmarking
 */

int main(void) {
    printf("=== Example 3: Performance Analysis ===\n\n");

    int n = 2, m = 1;
    double traj[] = {
        1.0, 0.0,   /* t=0 */
        0.8, -0.1,  /* t=0.1 */
        0.5, -0.2,  /* t=0.2 */
        0.3, -0.15, /* t=0.3 */
        0.1, -0.05, /* t=0.4 */
        0.05, 0.0,  /* t=0.5 */
        0.02, 0.01, /* t=0.6 */
        0.01, 0.005,/* t=0.7 */
        0.005,0.002,/* t=0.8 */
        0.0, 0.0    /* t=0.9 */
    };
    int traj_len = 10;
    double events[] = {0.0, 0.2, 0.5, 0.8};
    int evt_len = 4;
    double dt = 0.1;
    double K[] = {-1.0, -2.0};

    /* Compute performance metrics */
    EBC_Performance perf = ebc_compute_performance(
        traj, traj_len, n, events, evt_len, dt, 0.9, K, m, 1.0, 0.01);

    printf("Performance Metrics:\n");
    printf("  ISE (Integral Squared Error):     %.6f\n", perf.ise);
    printf("  IAE (Integral Absolute Error):    %.6f\n", perf.iae);
    printf("  ITAE (Time-weighted Abs Error):   %.6f\n", perf.itae);
    printf("  ISCI (Integral Squared Ctrl):     %.6f\n", perf.isci);
    printf("  Total events:                     %.0f\n", perf.total_events);
    printf("  Periodic equivalent:              %.0f\n", perf.periodic_equiv);
    printf("  Communication reduction:          %.1f%%\n",
           100.0 * perf.comm_reduction);
    printf("  Average inter-event time:         %.4f s\n", perf.avg_iet);
    printf("  Settling time (2%%):              %.4f s\n", perf.settling_time);
    printf("  Overshoot:                        %.2f%%\n",
           100.0 * perf.overshoot);
    printf("  Energy cost:                      %.4f\n", perf.energy_cost);
    printf("  Max state deviation:              %.6f\n", perf.max_state_dev);

    /* Stability analysis */
    printf("\nStability Analysis:\n");
    double A[] = {0.0, 1.0, -1.0, -2.0}; /* A+BK */
    double B[] = {0.0, 1.0};
    double sigma_crit = ebc_critical_sigma(A, B, K, 2, 1);
    printf("  Critical sigma:                   %.4f\n", sigma_crit);
    double tau_min = ebc_minimum_iet_linear(A, B, K, 2, 1, 0.5, 0.05);
    printf("  Minimum inter-event time bound:   %.6f s\n", tau_min);

    printf("\nExample 3 complete.\n");
    return 0;
}
