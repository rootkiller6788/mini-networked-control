/**
 * @file    example1_uniform_quantizer.c
 * @brief   Example: Uniform quantizer performance analysis
 *
 * Demonstrates the effect of bit resolution on quantization error
 * for a sinusoidal signal. Computes SQNR and compares with
 * theoretical predictions.
 *
 * L6 Canonical Problem: Quantized sinusoidal signal analysis
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quantized_control.h"

int main(void) {
    printf("=== Example 1: Uniform Quantizer Performance ===\n\n");

    QCQuantizer q;
    int bits_list[] = {2, 4, 6, 8, 10, 12};
    int n_bits = 6;
    int N = 1000;
    double *signal = malloc(N * sizeof(double));
    double *quantized = malloc(N * sizeof(double));

    /* Generate sinusoidal test signal */
    for (int i = 0; i < N; i++) {
        signal[i] = sin(2.0 * M_PI * i / 50.0);
    }

    printf("Bits | Step      | SQNR (dB) | Theoretical (dB) | ENOB\n");
    printf("-----|-----------|-----------|------------------|------\n");

    for (int b = 0; b < n_bits; b++) {
        int bits = bits_list[b];
        qc_quantizer_init(&q, QC_QTYPE_UNIFORM, bits);
        qc_quantizer_configure_range(&q, -1.0, 1.0);

        /* Quantize signal */
        for (int i = 0; i < N; i++) {
            quantized[i] = qc_quantize_scalar(&q, signal[i]);
        }

        /* Compute metrics */
        QCQuantizerMetrics metrics;
        qc_quantizer_metrics_compute(&q, signal, quantized, N, &metrics);
        double theory_sqnr = qc_uniform_sqnr_theoretical(bits, 0.7071, 1.0);

        printf(" %3d | %9.6f | %9.2f | %14.2f | %5.2f\n",
               bits, q.step, metrics.sqnr_db, theory_sqnr, metrics.enob);

        qc_quantizer_free(&q);
    }

    printf("\nObservation: Each additional bit improves SQNR by ~6 dB.\n");
    printf("This matches the theoretical formula: SQNR = 6.02*N + 1.76 dB.\n");

    free(signal);
    free(quantized);
    return 0;
}
