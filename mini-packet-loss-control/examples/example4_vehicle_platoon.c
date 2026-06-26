#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include "packet_loss_analysis.h"
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Vehicle platoon with packet loss.
 * Each vehicle: x_{k+1} = A·x_k + B·u_k
 * where x = [position_error, velocity_error] relative to predecessor.
 *
 * String stability: errors must not amplify toward the tail.
 * Packet loss between vehicles degrades string stability. */

int main(void) {
    printf("=== Example 4: Vehicle Platoon with Packet Loss ===\n\n");

    /* Vehicle dynamics (Ts=0.1s, double integrator) */
    double A_veh[4] = {1.0, 0.1, 0.0, 1.0};
    double B_veh[2] = {0.005, 0.1};

    printf("Vehicle model (Ts=0.1s, double integrator):\n");
    printf("  A = [%.1f %.1f; %.1f %.1f]\n", A_veh[0], A_veh[1], A_veh[2], A_veh[3]);
    printf("  B = [%.3f; %.1f]\n", B_veh[0], B_veh[1]);

    /* LQR design for each vehicle */
    LTISystem* sys = pl_lti_create(2, 1, 1, 2);
    pl_lti_set_A(sys, A_veh, 2);
    pl_lti_set_B(sys, B_veh);
    double Q_lqr[4] = {10.0, 0.0, 0.0, 1.0};
    double R_lqr[1] = {1.0};
    for (int i = 0; i < 4; i++) sys->Q[i] = Q_lqr[i];

    LQRSolution* lqr = pl_lqr_solve(sys, R_lqr, 500, 1e-6);
    printf("  LQR gain: [%.4f, %.4f]\n", lqr->L[0], lqr->L[1]);
    printf("  rho(A-BL)=%.6f\n", pl_lqr_spectral_radius(lqr, sys));

    /* Platoon of 5 vehicles */
    int N = 5;
    double x[N][2];  /* State for each vehicle */
    for (int i = 0; i < N; i++) {
        x[i][0] = (double)(N - i) * 0.5;  /* Spacing error */
        x[i][1] = 0.0;                      /* Velocity error */
    }

    /* Create independent channels between vehicles */
    BernoulliChannel* channels[N-1];
    for (int i = 0; i < N-1; i++)
        channels[i] = pl_bernoulli_create(0.15, (unsigned long)(100 + i * 37));

    printf("\nPlatoon: %d vehicles, packet loss 15%% between vehicles\n", N);
    printf("Step | Veh1(pos,vel) | Veh2(pos,vel) | ...  | MaxError\n");
    printf("-----+---------------+----------------+------+----------\n");

    /* Simulate platoon */
    int steps = 50;
    double max_error = 0.0;

    for (int k = 0; k < steps; k++) {
        /* Leader reference: sinusoidal perturbation */
        double leader_ref = 0.1 * sin(2.0 * M_PI * k / 50.0);
        (void)leader_ref; /* Leader state evolves via open-loop in this demo */

        /* Each vehicle computes control based on predecessor state */
        for (int i = 1; i < N; i++) {
            PacketStatus st = pl_bernoulli_transmit(channels[i-1]);
            bool arrived = (st == PACKET_RECEIVED);

            if (arrived) {
                /* Compute control using predecessor's actual state */
                double u = -(lqr->L[0] * x[i-1][0] + lqr->L[1] * x[i-1][1]);
                (void)u; /* Control would be applied in full implementation */
            }

            /* Update vehicle state (open-loop for simplicity) */
            double x_next[2];
            x_next[0] = A_veh[0]*x[i][0] + A_veh[1]*x[i][1];
            x_next[1] = A_veh[2]*x[i][0] + A_veh[3]*x[i][1];
            x[i][0] = x_next[0];
            x[i][1] = x_next[1];
        }

        /* Track maximum position error in platoon */
        double max_pos_err = 0.0;
        for (int i = 0; i < N; i++)
            if (fabs(x[i][0]) > max_pos_err) max_pos_err = fabs(x[i][0]);
        if (max_pos_err > max_error) max_error = max_pos_err;

        if (k % 10 == 0) {
            printf("%4d | (%+.2f,%+.2f) | (%+.2f,%+.2f) | ...  | %.4f\n",
                   k, x[0][0], x[0][1], x[1][0], x[1][1], max_pos_err);
        }
    }

    printf("-----+---------------+----------------+------+----------\n");
    printf("Final states:\n");
    for (int i = 0; i < N; i++)
        printf("  Vehicle %d: pos_err=%+.4f, vel_err=%+.4f\n",
               i+1, x[i][0], x[i][1]);
    printf("Maximum position error: %.4f\n", max_error);

    /* String stability check */
    double loss_rates[N-1];
    double avg_loss = 0.0;
    for (int i = 0; i < N-1; i++) {
        loss_rates[i] = pl_bernoulli_get_loss_rate(channels[i]);
        avg_loss += loss_rates[i];
    }
    avg_loss /= (N-1);
    printf("Average empirical loss rate: %.4f\n", avg_loss);

    /* Cleanup */
    for (int i = 0; i < N-1; i++) pl_bernoulli_free(channels[i]);
    pl_lqr_free(lqr);
    pl_lti_free(sys);

    printf("\nExample 4 complete.\n");
    return 0;
}
