/**
 * @file    example2_log_quantizer.c
 * @brief   Example: Logarithmic quantizer and sector bound analysis
 *
 * Demonstrates the logarithmic quantizer's constant relative error
 * property and the sector bound delta = (1-rho)/(1+rho).
 *
 * L6 Canonical Problem: Logarithmic quantization for control
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quantized_control.h"

int main(void) {
    printf("=== Example 2: Logarithmic Quantizer Analysis ===\n\n");

    double rho_values[] = {0.3, 0.5, 0.7, 0.9};
    int n_rho = 4;

    printf("Logarithmic Quantizer: q(x) = sign(x) * u_min * rho^{-i}\n");
    printf("Sector bound: |e| <= delta * |x|, delta = (1-rho)/(1+rho)\n\n");

    printf(" rho  | delta  | Relative Error (tested)\n");
    printf("------|--------|------------------------\n");

    for (int k = 0; k < n_rho; k++) {
        double rho = rho_values[k];
        double delta = qc_log_sector_delta(rho);

        QCLogQuantizer lq;
        qc_log_quantizer_init(&lq, rho, 0.001);

        /* Test relative error at various magnitudes */
        double test_vals[] = {0.002, 0.01, 0.1, 1.0};
        double max_rel_err = 0.0;
        for (int i = 0; i < 4; i++) {
            double x = test_vals[i];
            double xq = qc_log_quantize(&lq, x);
            if (fabs(x) > lq.deadzone) {
                double rel_err = fabs(xq - x) / fabs(x);
                if (rel_err > max_rel_err) max_rel_err = rel_err;
            }
        }

        printf(" %.2f | %.4f | %.4f\n", rho, delta, max_rel_err);
        qc_log_quantizer_free(&lq);
    }

    printf("\nInterpretation:\n");
    printf("- Smaller rho (denser quantization) → smaller delta → tighter bound\n");
    printf("- The relative error is CONSTANT regardless of signal magnitude\n");
    printf("- This enables robust stability analysis via sector bound\n");
    printf("- Application: quantized feedback stabilization (Elia & Mitter 2001)\n");

    return 0;
}
