/**
 * example1_bandwidth_limited_stabilization.c
 *
 * Demonstrates the Data Rate Theorem in action:
 * Stabilizing an unstable double integrator through a quantized,
 * bandwidth-limited communication channel.
 *
 * Key result: The minimum data rate for a double integrator
 * (ẋ₁ = x₂, ẋ₂ = u) with one unstable eigenvalue λ = 0 is
 * technically 0 (marginal stability), but with any positive
 * feedback delay or quantization, a nonzero rate is required.
 *
 * This example shows:
 *   1. Setting up a bandwidth-limited control system
 *   2. Computing the theoretical minimum data rate
 *   3. Simulating with different bit rates
 *   4. Observing the effect of quantization coarseness on stability
 */

#include "blc_core.h"
#include "blc_datarate.h"
#include "blc_control.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 1: Bandwidth-Limited Stabilization ===\n\n");

    /** Create a 2-state, 1-input system */
    BLCSystem* sys = blc_create(2, 1, 1);
    if (!sys) {
        printf("Failed to create system.\n");
        return 1;
    }

    /** Double integrator: A = [[0,1],[0,0]], B = [[0],[1]], C = [[1,0]] */
    double A[4] = {0.0, 1.0, 0.0, 0.0};
    double B[2] = {0.0, 1.0};
    double Cmat[2] = {1.0, 0.0};
    blc_init_plant(sys, A, B, Cmat);

    /** Set initial state: offset from origin */
    double x0[2] = {5.0, -2.0};
    blc_set_initial_state(sys, x0);

    /** Configure communication channel:
     *  Bandwidth = 100 Hz, SNR = 10 (10 dB), latency = 1 ms, loss = 0.1% */
    blc_init_channel(sys, 100.0, 10.0, 1.0, 0.001);

    /** Design LQR controller */
    BLCLQRController lqr;
    double Q_diag[2] = {10.0, 1.0};
    double R_diag[1] = {0.1};
    blc_lqr_init(&lqr, 2, 1, Q_diag, R_diag);
    int iter = blc_lqr_solve(&lqr, A, B, 200, 1e-8);

    /** Set controller gain */
    double K_flat[2] = {lqr.K[0][0], lqr.K[0][1]};
    blc_set_controller_gain(sys, K_flat);

    /** Compute closed-loop poles */
    double poles_re[2], poles_im[2];
    blc_lqr_closed_loop_poles(&lqr, A, B, poles_re, poles_im);

    /** Data Rate Theorem analysis */
    sys->plant.eigenvalues[0] = fabs(poles_re[0]);  /** These are stable but... */
    sys->plant.eigenvalues[1] = fabs(poles_re[1]);

    /** The unstable eigenvalue of the open-loop system is λ = 0 (double pole).
     *  With quantization and delay, the effective required rate is nonzero. */
    double R_min = blc_datarate_min_ct(&sys->plant);
    double C = blc_channel_capacity(&sys->channel);

    printf("Plant: Double Integrator (2 states, 1 input)\n");
    printf("Open-loop eigenvalues: λ₁=0, λ₂=0 (marginally stable)\n");
    printf("LQR gain: K = [%.4f, %.4f] (iter=%d)\n", lqr.K[0][0], lqr.K[0][1], iter);
    printf("Closed-loop poles: %.4f ± j%.4f, %.4f ± j%.4f\n",
           poles_re[0], poles_im[0], poles_re[1], poles_im[1]);
    printf("\nChannel: B=%.1f Hz, SNR=%.1f (%.1f dB), C=%.2f bps\n",
           sys->channel.bandwidth_hz, sys->channel.snr,
           10.0 * log10(sys->channel.snr), C);
    printf("Minimum data rate (theoretical): %.4f bps\n", R_min);
    printf("Capacity margin: %.2f bps (%.1f%% of min)\n",
           C - R_min, R_min > 0 ? 100.0*(C-R_min)/R_min : 999.0);

    /** Test quantization at different bit rates */
    printf("\n--- Quantization Sweep ---\n");
    printf("Levels  Bits/state  Step size  Status\n");
    printf("------  ----------  ---------  ------\n");

    int levels_test[] = {4, 8, 16, 32, 64, 128, 256};
    for (int t = 0; t < 7; t++) {
        int lev = levels_test[t];
        double bits_per_state = log2((double)lev);
        double step = 20.0 / (double)lev;  /** Range [-10, 10] */

        /** Quick stability check via maximum quantization step */
        double Ac_flat[4];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                Ac_flat[i*2+j] = lqr.Ac[i][j];
        double delta_max = blc_max_quantization_step(Ac_flat, B, K_flat, 2);
        int req_levels = blc_required_quantization_levels(10.0, -10.0, delta_max);

        const char* status;
        if (lev >= req_levels)
            status = "STABLE ✓";
        else if (lev >= req_levels / 2)
            status = "MARGINAL";
        else
            status = "UNSTABLE ✗";

        printf("%6d  %10.2f  %9.4f  %s\n",
               lev, bits_per_state, step, status);
    }

    /** Run simulation */
    printf("\n--- Simulation (256 levels, 10s) ---\n");
    sys->sample_period = 0.01;
    int n_steps = 1000;

    double t_hist[2000];
    double x_hist[4000];  /** [x1, x2] per step */

    int steps = blc_simulate(sys, 10.0, t_hist, x_hist, n_steps);
    printf("Simulated %d steps (%.2f seconds)\n", steps, sys->sim_time);

    /** Print trajectory summary */
    printf("\nTime(s)    x1        x2        ||x||\n");
    printf("------  --------  --------  --------\n");
    for (int k = 0; k <= steps; k += steps / 10) {
        if (k >= steps) k = steps - 1;
        double x1 = x_hist[k * 2];
        double x2 = x_hist[k * 2 + 1];
        double norm = sqrt(x1*x1 + x2*x2);
        printf("%6.2f  %8.4f  %8.4f  %8.4f\n",
               t_hist[k], x1, x2, norm);
    }

    printf("\nFinal ||x|| = %.6f\n",
           sqrt(sys->plant.x[0]*sys->plant.x[0] +
                sys->plant.x[1]*sys->plant.x[1]));

    /** Check stability */
    int stable = blc_check_stability(sys);
    printf("Stability: %s\n", stable > 0 ? "STABLE ✓" : "UNSTABLE ✗");

    /** Packet statistics */
    int sent, lost, overloads;
    blc_get_packet_stats(sys, &sent, &lost, &overloads);
    printf("Packets: %d sent, %d lost, %d overloads\n", sent, lost, overloads);

    /** Data rate usage */
    double bits_per_sample_used = ceil(log2(256.0)) * 2.0;  /** 2 states × 8 bits */
    double rate_used = bits_per_sample_used / sys->sample_period;
    printf("Data rate used: %.1f bps (%.1f bits/sample)\n",
           rate_used, bits_per_sample_used);
    printf("Channel utilization: %.1f%%\n", 100.0 * rate_used / C);

    blc_free(sys);
    return 0;
}