/**
 * @file    demo_quantization.c
 * @brief   Demo: Interactive quantized control visualization
 *
 * Demonstrates the effect of different quantizer types and bit depths
 * on a simple first-order control system. Outputs data for plotting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quantized_control.h"

int main(void) {
    printf("=== Quantized Control Demo ===\n\n");

    QCSystem sys;
    qc_system_init(&sys);
    qc_system_configure(&sys, 1, 1, 1);

    double A[1] = {-1.0};
    double B[1] = {1.0};
    double C[1] = {1.0};
    qc_system_set_A(&sys, A);
    qc_system_set_B(&sys, B);
    qc_system_set_C(&sys, C);

    int bits_list[] = {2, 4, 8, 16};
    int n_configs = 4;

    printf("Comparing quantization resolutions for first-order system\n");
    printf("System: dx/dt = -x + u, C = 1\n\n");

    for (int cfg = 0; cfg < n_configs; cfg++) {
        int bits = bits_list[cfg];
        qc_quantizer_init(&sys.input_quantizer, QC_QTYPE_UNIFORM, bits);
        qc_quantizer_configure_range(&sys.input_quantizer, -1.0, 1.0);

        double x = 1.0;
        double dt = 0.01;
        int steps = 200;

        printf("--- %d-bit Quantizer (step=%.6f) ---\n", bits, sys.input_quantizer.step);
        printf("  t      x           u (ideal)   u (quantized)\n");

        for (int k = 0; k <= steps; k += 40) {
            double t = k * dt;
            double u_ideal = -2.0 * x;
            double u_q = qc_quantize_scalar(&sys.input_quantizer, u_ideal);
            printf("  %5.2f  %10.6f  %10.6f  %10.6f\n", t, x, u_ideal, u_q);

            /* Euler step */
            x += dt * (-x + u_q);
        }

        /* Steady-state analysis */
        printf("  Final: x=%.6f | Expected steady-state near zero\n\n", x);
        qc_quantizer_free(&sys.input_quantizer);
    }

    printf("Observation: Higher bit depths → finer control → smaller steady-state error.\n");
    printf("This is the 'practical stability' property of quantized control.\n");

    return 0;
}
