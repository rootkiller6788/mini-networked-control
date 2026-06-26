#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include "packet_loss_estimator.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 2: Gilbert-Elliott Markov Channel & Kalman Filter ===\n\n");

    /* Create Gilbert-Elliott channel: 5% Good→Bad, 30% Bad→Good,
     * 1% loss in Good, 40% loss in Bad */
    GilbertElliottChannel* ge = pl_gilbert_elliott_create(0.05, 0.30, 0.01, 0.40, 42);
    printf("Gilbert-Elliott channel:\n");
    printf("  P(Good->Bad)=%.2f, P(Bad->Good)=%.2f\n", ge->p_gb, ge->p_bg);
    printf("  Loss(Good)=%.2f, Loss(Bad)=%.2f\n", ge->loss_rate_good, ge->loss_rate_bad);
    printf("  Steady-state: pi_G=%.3f, pi_B=%.3f, avg loss=%.3f\n",
           ge->steady_p_good, ge->steady_p_bad, ge->steady_loss_rate);
    printf("  Burstiness: %.3f (>1 means bursty)\n",
           pl_gilbert_elliott_burstiness(ge));

    /* Create an intermittent Kalman filter for a 2-state system */
    double A[4] = {1.0, 0.05, 0.0, 0.95};
    double C[2] = {1.0, 0.0};
    double Q[4] = {0.001, 0.0, 0.0, 0.001};
    double R[1] = {0.01};
    double x0[2] = {0.0, 0.0};
    double P0[4] = {1.0, 0.0, 0.0, 1.0};

    IntermittentKalmanFilter* ikf = pl_ikf_create(A, C, Q, R, 2, 1, x0, P0, 0.85);

    printf("\nIntermittent Kalman Filter:\n");
    printf("  Arrival probability gamma=%.2f\n", ikf->arrival_probability);
    printf("  Critical gamma bounds: [%.4f, %.4f]\n",
           ikf->critical_gamma_lower, ikf->critical_gamma_upper);
    printf("  Stable: %s\n", pl_ikf_is_stable(ikf) ? "YES" : "NO");

    /* Simulate 200 time steps */
    printf("\nSimulating 200 steps...\n");
    double y_true[1] = {1.0};
    for (int k = 0; k < 200; k++) {
        PacketStatus st = pl_gilbert_elliott_transmit(ge);
        bool arrived = (st == PACKET_RECEIVED);
        pl_ikf_step(ikf, y_true, arrived);
    }

    /* Print results */
    pl_ikf_expected_covariance(ikf);
    printf("\nResults after 200 steps:\n");
    printf("  Max consecutive losses: %d\n", ikf->max_consecutive_losses);
    printf("  Trace(E[P]): %.6f\n", pl_ikf_trace_expected(ikf));
    printf("  State estimate: [%.4f]\n", pl_ikf_get_estimate(ikf)[0]);
    printf("  Empirical loss rate: %.3f\n",
           (double)ge->total_lost / (double)ge->total_sent);
    printf("  Channel state: %s\n", ge->current_state == 0 ? "Good" : "Bad");

    pl_ikf_free(ikf);
    pl_gilbert_elliott_free(ge);
    printf("\nExample 2 complete.\n");
    return 0;
}
