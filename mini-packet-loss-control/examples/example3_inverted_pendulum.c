#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include "packet_loss_analysis.h"
#include "packet_loss_predictor.h"
#include <stdio.h>
#include <math.h>

/* Discretized inverted pendulum on a cart (linearized, Ts=0.01s)
 * States: [angle, angular_velocity]
 * M=0.5kg, m=0.2kg, l=0.3m, g=9.81
 *
 * Continuous: A_c = [0, 1; (M+m)g/(Ml), 0], B_c = [0; -1/(Ml)]
 * Discrete (Euler): A = I + Ts·A_c, B = Ts·B_c */

int main(void) {
    printf("=== Example 3: Inverted Pendulum over Lossy Network ===\n\n");

    /* Discretized inverted pendulum */
    double A[4] = {1.0000, 0.0100, 0.4578, 1.0000};
    double B[2] = {0.0000, -0.0667};

    printf("Inverted Pendulum Model (linearized, upright, Ts=0.01s):\n");
    printf("  A = [%.4f %.4f; %.4f %.4f]\n", A[0], A[1], A[2], A[3]);
    printf("  B = [%.4f; %.4f]\n", B[0], B[1]);
    printf("  Open-loop eigenvalues: unstable (pendulum falls)\n");

    /* LQR design */
    LTISystem* sys = pl_lti_create(2, 1, 1, 2);
    pl_lti_set_A(sys, A, 2);
    pl_lti_set_B(sys, B);
    double Q_lqr[4] = {100.0, 0.0, 0.0, 1.0};
    double R_lqr[1] = {0.1};
    for (int i = 0; i < 4; i++) sys->Q[i] = Q_lqr[i];

    LQRSolution* lqr = pl_lqr_solve(sys, R_lqr, 1000, 1e-8);
    printf("  LQR solution: %d iterations, rho(A-BL)=%.6f\n",
           lqr->iterations, pl_lqr_spectral_radius(lqr, sys));
    printf("  LQR gain: [%.4f, %.4f]\n", lqr->L[0], lqr->L[1]);

    /* Create Packetized Predictive Control with horizon H=5 */
    double x_init[2] = {0.2, 0.0}; /* 0.2 rad initial angle (~11.5 deg) */
    PacketizedPredictiveControl* ppc = pl_ppc_create(5, A, B, lqr->L, 2, 1, x_init);
    printf("\nPacketized Predictive Control (horizon H=5):\n");

    /* Simulate: introduce varying packet loss */
    BernoulliChannel* ch = pl_bernoulli_create(0.3, 99);
    double x[2] = {0.2, 0.0};
    double x_hist[100]; (void)x_hist;
    int hist_idx = 0; (void)hist_idx;

    printf("\nSimulating 100 steps with 30%% packet loss...\n");
    printf("Step | Angle   | AngVel  | Control  | Status\n");
    printf("-----+---------+---------+----------+--------\n");

    for (int k = 0; k < 100; k++) {
        PacketStatus st = pl_bernoulli_transmit(ch);
        bool arrived = (st == PACKET_RECEIVED);

        pl_ppc_update_state(ppc, x);
        const double* u = pl_ppc_consume(ppc, arrived);

        /* Simulate plant: x_next = A·x + B·u */
        double x_next[2];
        x_next[0] = A[0]*x[0] + A[1]*x[1] + B[0]*u[0];
        x_next[1] = A[2]*x[0] + A[3]*x[1] + B[1]*u[0];

        if (k % 20 == 0) {
            printf("%4d | %+.4f | %+.4f | %+.4f | %s\n",
                   k, x[0], x[1], u[0],
                   arrived ? "received" : "LOST");
        }

        x[0] = x_next[0]; x[1] = x_next[1];
        if (hist_idx < 100) x_hist[hist_idx++] = x[0];
    }

    printf("-----+---------+---------+----------+--------\n");
    printf("Final state: angle=%.4f, angvel=%.4f\n", x[0], x[1]);
    printf("Max consecutive losses: %d (tolerated: %d)\n",
           ppc->max_consecutive_losses_seen, ppc->horizon);
    printf("Buffer exhausted: %s\n", ppc->buffer_exhausted ? "YES" : "NO");
    printf("Avg buffer utilization: %.3f\n", ppc->avg_buffer_utilization);

    /* Stability analysis */
    printf("\nStability Analysis:\n");
    double A_cl[4]; (void)A_cl;
    double BL[2][2] = {{0}};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int l = 0; l < 1; l++)
                BL[i][j] += B[i] * lqr->L[l*2+j];
    for (int i = 0; i < 4; i++) A_cl[i] = A[i] - ((double*)BL)[i];

    CriticalProbabilityAnalysis* cpa = pl_critical_prob_analyze(
        A, B, NULL, lqr->L, 2, 1, 1);
    printf("  p_critical(sensor):  %.4f\n", cpa->p_c_sensor);
    printf("  p_critical(actuator): %.4f\n", cpa->p_c_actuator);
    printf("  p_critical(joint):   %.4f\n", cpa->p_c_joint);

    pl_critical_prob_free(cpa);
    pl_ppc_free(ppc);
    pl_lqr_free(lqr);
    pl_lti_free(sys);
    pl_bernoulli_free(ch);

    printf("\nExample 3 complete.\n");
    return 0;
}
