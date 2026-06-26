/**
 * @file    bench_quantizer.c
 * @brief   Benchmarks for quantizer performance
 *
 * Measures throughput of quantizer operations in million samples per second.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "quantized_control.h"

static double get_time_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

static void bench_uniform_quantizer(int N) {
    QCQuantizer q;
    qc_quantizer_init(&q, QC_QTYPE_UNIFORM, 8);
    qc_quantizer_configure_range(&q, -1.0, 1.0);

    double *input = malloc(N * sizeof(double));
    double *output = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) input[i] = sin(i * 0.001);

    double t0 = get_time_sec();
    for (int i = 0; i < N; i++) output[i] = qc_quantize_scalar(&q, input[i]);
    double t1 = get_time_sec();

    double msps = N / ((t1 - t0) * 1e6);
    printf("  Uniform (8-bit): %.2f Msamples/sec\n", msps);

    free(input); free(output);
    qc_quantizer_free(&q);
}

static void bench_log_quantizer(int N) {
    QCLogQuantizer lq;
    qc_log_quantizer_init(&lq, 0.5, 0.001);

    double *input = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) input[i] = exp(sin(i * 0.001));

    double t0 = get_time_sec();
    for (int i = 0; i < N; i++) qc_log_quantize(&lq, input[i]);
    double t1 = get_time_sec();

    double msps = N / ((t1 - t0) * 1e6);
    printf("  Logarithmic: %.2f Msamples/sec\n", msps);

    free(input);
    qc_log_quantizer_free(&lq);
}

static void bench_huffman(int N) {
    double probs[4] = {0.4, 0.3, 0.2, 0.1};
    int symbols[4] = {0, 1, 2, 3};
    QCHuffmanTree tree;
    qc_huffman_build(&tree, probs, symbols, 4);

    QCBitWriter bw;
    qc_bit_writer_init(&bw, N * 4);

    double t0 = get_time_sec();
    for (int i = 0; i < N; i++) qc_huffman_encode(&tree, i % 4, &bw);
    double t1 = get_time_sec();

    double msps = N / ((t1 - t0) * 1e6);
    printf("  Huffman encode: %.2f Msymbols/sec\n", msps);

    qc_bit_writer_free(&bw);
    qc_huffman_free(&tree);
}

int main(void) {
    printf("=== Quantizer Benchmarks ===\n\n");
    int N = 1000000;
    printf("N = %d samples\n\n", N);

    bench_uniform_quantizer(N);
    bench_log_quantizer(N / 100);
    bench_huffman(N / 10);

    printf("\nDone.\n");
    return 0;
}
