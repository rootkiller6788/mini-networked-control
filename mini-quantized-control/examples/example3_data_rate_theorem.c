/**
 * @file    example3_data_rate_theorem.c
 * @brief   Example: Data Rate Theorem verification
 *
 * Demonstrates the Nair-Evans Data Rate Theorem by computing the
 * minimum data rate required to stabilize unstable LTI systems
 * over quantized feedback channels.
 *
 * L6 Canonical Problem: Minimum data rate for stabilization
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quantized_control.h"

int main(void) {
    printf("=== Example 3: Data Rate Theorem ===\n\n");
    printf("Theorem (Nair & Evans, 2004):\n");
    printf("For x_{k+1} = A x_k + B u_k with quantized state feedback,\n");
    printf("stabilization requires R > sum_i max(0, log2|lambda_i(A)|).\n\n");

    /* Test case 1: Scalar unstable system */
    printf("--- Case 1: Scalar system x_{k+1} = a * x_k + u_k ---\n");
    double a_vals[] = {1.1, 1.5, 2.0, 3.0, 5.0};
    int n_a = 5;

    printf("  a   | |lambda| | Min Rate (bits) | Verdict\n");
    printf("------|---------|-----------------|--------\n");
    for (int i = 0; i < n_a; i++) {
        double a = a_vals[i];
        double A[1] = {a};
        double min_rate = qc_data_rate_theoretical_min(A, 1);
        double available = 2.0; /* assume 2 bits/sample available */
        int stab = qc_data_rate_is_stabilizable(A, 1, available, 0.0);
        printf(" %4.1f | %7.1f | %15.4f | %s\n",
               a, fabs(a), min_rate, stab ? "STABILIZABLE" : "NOT STABILIZABLE");
    }

    /* Test case 2: 2D system with mixed stability */
    printf("\n--- Case 2: 2D system with mixed stability ---\n");
    double A2[4] = {2.0, 0.5, 0.0, 0.8};
    /* Eigenvalues: 2.0 (unstable, needs log2(2)=1 bit), 0.8 (stable, 0 bits) */
    double min_rate_2d = qc_data_rate_theoretical_min(A2, 2);
    printf("A = [[2.0, 0.5], [0, 0.8]]\n");
    printf("Minimum rate: %.4f bits/sample (expected: log2(2) = 1.0000)\n",
           min_rate_2d);
    printf("Verification: %.4f ≈ 1.0000\n", min_rate_2d);

    /* Test case 3: Complex eigenvalues */
    printf("\n--- Case 3: Oscillatory unstable system ---\n");
    double A3[4] = {0.0, -2.0, 2.0, 0.0};
    /* Rotation-scaling: eigenvalues ±j2, |lambda|=2 */
    double min_rate_3d = qc_data_rate_theoretical_min(A3, 2);
    printf("A = [[0, -2], [2, 0]] → eigenvalues = ±2j, |lambda| = 2\n");
    printf("Minimum rate: %.4f bits/sample\n", min_rate_3d);

    /* Channel capacity demonstration */
    printf("\n--- Case 4: Channel Capacity ---\n");
    printf("AWGN channel, B=1 Hz:\n");
    printf("  SNR=0dB  → C=%.4f bps\n", qc_channel_capacity_awgn(1.0, 1.0));
    printf("  SNR=10dB → C=%.4f bps\n",
           qc_channel_capacity_awgn(1.0, pow(10.0, 1.0)));
    printf("  SNR=20dB → C=%.4f bps\n",
           qc_channel_capacity_awgn(1.0, pow(10.0, 2.0)));

    printf("\nKey Insight:\n");
    printf("- Each unstable eigenvalue contributes log2|lambda| to min rate\n");
    printf("- Stable eigenvalues contribute 0 (no information needed)\n");
    printf("- If channel capacity < min rate → stabilization impossible\n");

    return 0;
}
