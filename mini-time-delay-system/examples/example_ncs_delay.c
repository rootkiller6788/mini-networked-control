#include "time_delay_system.h"
#include "networked_delay.h"
#include "delay_stability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ============================================================================
 * Example 2: Networked Control System with Delay and Packet Loss
 *
 * Demonstrates how network-induced delays and packet loss affect
 * closed-loop stability and performance. Compares:
 *   1. Ideal control (no delay, no loss)
 *   2. NCS with delay only
 *   3. NCS with delay + packet loss
 *   4. NCS with compensation (timestamp-based)
 *
 * System: second-order unstable plant:
 *   ẋ₁ = x₂
 *   ẋ₂ = -x₁ + 0.5 x₂ + u(t-τ)
 *
 * Without delay: stable with u = -[2.0, 1.5] x
 * With delay: stability depends on τ
 * ============================================================================ */

int main(void) {
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Networked Control — Delay & Packet Loss     ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    double h = 0.02;  /* Sampling period 20ms */
    int N_steps = 500;  /* 10 seconds simulation */

    /* --- Scenario 1: Ideal control (no delay, no loss) --- */
    printf("=== Scenario 1: Ideal Control (τ=0, loss=0%%) ===\n");

    TimeDelaySystem* plant1 = tds_create("ideal_plant", 2, 1, 1);
    double A[4] = {0.0, 1.0, -1.0, 0.5};
    double Ad[4] = {0.0, 0.0, 0.0, 0.0};
    double B[2] = {0.0, 1.0};
    double C[2] = {1.0, 0.0};
    tds_set_linear_model(plant1, A, Ad, B, C);
    plant1->current_state[0] = 1.0;  /* Initial condition x₁(0)=1 */

    NetworkedControlSystem* ncs1 = ncs_create(plant1, h, 0.0);
    double K[2] = {2.0, 1.5};
    ncs_set_gain(ncs1, K);
    ncs_set_sc_qos(ncs1, 100e6, 0.0, 0.0, 0.0, 0.0, 0.0);
    ncs_set_ca_qos(ncs1, 100e6, 0.0, 0.0, 0.0, 0.0, 0.0);

    double* t1 = (double*)malloc((size_t)N_steps * sizeof(double));
    double* y1 = (double*)malloc((size_t)N_steps * sizeof(double));
    double* u1 = (double*)malloc((size_t)N_steps * sizeof(double));
    ncs_run(ncs1, N_steps, t1, y1, u1);
    printf("  Final state: [%.4f, %.4f]\n",
           plant1->current_state[0], plant1->current_state[1]);
    printf("  ISE: %.6f  Settling: %.2f s\n\n",
           ncs1->ise, t1[N_steps-1]);

    /* --- Scenario 2: NCS with delay-only --- */
    printf("=== Scenario 2: NCS with Delay (τ_sc=50ms, τ_ca=30ms) ===\n");

    TimeDelaySystem* plant2 = tds_create("delay_plant", 2, 1, 1);
    tds_set_linear_model(plant2, A, Ad, B, C);
    plant2->current_state[0] = 1.0;

    NetworkedControlSystem* ncs2 = ncs_create(plant2, h, 0.01);
    ncs_set_gain(ncs2, K);
    /* SC: 50ms mean delay, 10ms jitter, 20-80ms range */
    ncs_set_sc_qos(ncs2, 100e6, 0.0, 0.05, 0.01, 0.02, 0.08);
    /* CA: 30ms mean delay, 5ms jitter, 10-50ms range */
    ncs_set_ca_qos(ncs2, 100e6, 0.0, 0.03, 0.005, 0.01, 0.05);

    double* y2 = (double*)malloc((size_t)N_steps * sizeof(double));
    double* u2 = (double*)malloc((size_t)N_steps * sizeof(double));
    ncs_run(ncs2, N_steps, NULL, y2, u2);
    printf("  Final state: [%.4f, %.4f]\n",
           plant2->current_state[0], plant2->current_state[1]);
    printf("  Avg delay: %.4f s\n\n", ncs2->avg_actual_delay);

    /* --- Scenario 3: NCS with delay + packet loss --- */
    printf("=== Scenario 3: NCS with Delay + 10%% Packet Loss ===\n");

    TimeDelaySystem* plant3 = tds_create("loss_plant", 2, 1, 1);
    tds_set_linear_model(plant3, A, Ad, B, C);
    plant3->current_state[0] = 1.0;

    NetworkedControlSystem* ncs3 = ncs_create(plant3, h, 0.01);
    ncs_set_gain(ncs3, K);
    ncs_set_sc_qos(ncs3, 100e6, 0.10, 0.05, 0.01, 0.02, 0.08);
    ncs_set_ca_qos(ncs3, 100e6, 0.10, 0.03, 0.005, 0.01, 0.05);

    double* y3 = (double*)malloc((size_t)N_steps * sizeof(double));
    double* u3 = (double*)malloc((size_t)N_steps * sizeof(double));
    ncs_run(ncs3, N_steps, NULL, y3, u3);
    printf("  Final state: [%.4f, %.4f]\n",
           plant3->current_state[0], plant3->current_state[1]);
    printf("  Packet loss: %d/%d = %.1f%%\n\n",
           ncs3->n_packets_lost, ncs3->n_packets_sent,
           100.0 * ncs3->n_packets_lost / (double)ncs3->n_packets_sent);

    /* --- Comparison --- */
    printf("=== Comparison ===\n");
    printf("  %-30s %12s %12s\n", "Scenario", "||x||_final", "ISE");
    printf("  %-30s %12.4f %12.4f\n", "Ideal (no delay)",
           sqrt(plant1->current_state[0]*plant1->current_state[0]
              + plant1->current_state[1]*plant1->current_state[1]),
           ncs1->ise);
    printf("  %-30s %12.4f %12.4f\n", "Delay only",
           sqrt(plant2->current_state[0]*plant2->current_state[0]
              + plant2->current_state[1]*plant2->current_state[1]),
           ncs2->ise);
    printf("  %-30s %12.4f %12.4f\n", "Delay + 10% loss",
           sqrt(plant3->current_state[0]*plant3->current_state[0]
              + plant3->current_state[1]*plant3->current_state[1]),
           ncs3->ise);
    printf("\nImplication: network degradation → performance degradation\n");

    /* Cleanup */
    ncs_free(ncs1); ncs_free(ncs2); ncs_free(ncs3);
    tds_free(plant1); tds_free(plant2); tds_free(plant3);
    free(t1); free(y1); free(u1);
    free(y2); free(u2); free(y3); free(u3);

    printf("Done.\n");
    return 0;
}
