#include "ebc_core.h"
#include "ebc_trigger.h"
#include "ebc_stability.h"
#include "ebc_performance.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
 * example1_etc_stabilization.c -- Event-Triggered Stabilization Demo
 *
 * Demonstrates event-triggered control for a linear system:
 *   dx/dt = A*x + B*u,  u = K*x(t_k)
 *
 * The system is a simple 2D double integrator:
 *   A = [0 1; 0 0]  (double integrator)
 *   B = [0; 1]
 *   K = [-1 -2]      (pole placement at -1, -1: s^2 + 3s + 1)
 *
 * This example shows:
 *   - System creation
 *   - Trigger parameter configuration
 *   - Event-based simulation
 *   - Performance comparison vs periodic
 *
 * L6: Canonical problem -- stabilization with limited communication
 */

static void double_integrator_dynamics(double t, const double* x,
    const double* u, int n, double* dx, void* ctx) {
    (void)t; (void)ctx;
    /* dx = A*x + B*u */
    dx[0] = x[1];           /* x1' = x2 */
    dx[1] = 0.0 + u[0];     /* x2' = u (double integrator) */
}

int main(void) {
    printf("=== Example 1: Event-Triggered Stabilization ===\n\n");

    int n = 2, m = 1;
    double x0[] = {1.0, 0.0};
    double K[] = {-1.0, -2.0};  /* feedback gain */

    /* Create ETC system */
    EBC_System* sys = ebc_system_create(n, m, EBC_CONTINUOUS_ETC);
    ebc_system_set_state(sys, x0);

    /* Configure trigger */
    EBC_TriggerParams tp = ebc_trigger_mixed(0.5, 0.05);
    printf("Trigger: sigma=%.3f, epsilon=%.3f\n", tp.sigma, tp.epsilon);

    /* Configure controller */
    EBC_Controller ctrl;
    memset(&ctrl, 0, sizeof(ctrl));
    ctrl.K = K;
    ctrl.n = n;
    ctrl.m = m;

    /* Simulate */
    double T = 5.0, dt = 0.001;
    double *traj = NULL, *events = NULL;
    int traj_len = 0, evt_len = 0;

    int ret = ebc_simulate(sys, &ctrl, T, dt, &tp,
                            &traj, &traj_len, &events, &evt_len);

    if (ret == 0) {
        printf("Simulation complete: %d steps, %d events\n", traj_len, evt_len);
        printf("Events triggered at times: ");
        for (int i = 0; i < evt_len && i < 10; i++)
            printf("%.3f ", events[i]);
        printf("\n");
        printf("Total events: %d (periodic equivalent: %.0f)\n",
               sys->event_count, T / dt);
        printf("Communication reduction: %.1f%%\n",
               100.0 * (1.0 - sys->event_count / (T / dt)));
        printf("Minimum IET: %.6f s\n", sys->min_iet);
        printf("Average IET: %.6f s\n", ebc_average_iet(events, evt_len));

        /* Print trajectory */
        printf("\nTrajectory (x1, x2):\n");
        int step = (traj_len > 50) ? traj_len / 50 : 1;
        for (int k = 0; k < traj_len; k += step) {
            printf("  t=%.3f: [%.4f, %.4f]\n",
                   k * dt, traj[k * n], traj[k * n + 1]);
        }
        free(traj);
        free(events);
    }

    ebc_system_free(sys);
    printf("\nExample 1 complete.\n");
    return 0;
}
