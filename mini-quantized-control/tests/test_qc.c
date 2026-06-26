/**
 * @file    test_qc.c
 * @brief   Tests for mini-quantized-control module
 *
 * Tests cover:
 *   - System initialization and configuration
 *   - Uniform quantizer accuracy bounds
 *   - Logarithmic quantizer sector bound
 *   - Data rate theorem minimum rate computation
 *   - Encoder/decoder roundtrip
 *   - Closed-loop simulation with quantization
 *   - Dynamic zoom quantizer
 *   - Huffman coding
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "quantized_control.h"
#include "qc_quantizer.h"
#include "qc_data_rate.h"
#include "qc_encoder.h"

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define FAIL(msg) do { printf("FAILED: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)

/* ================================================================
 * Test 1: System Initialization
 * ================================================================ */
static void test_system_init(void) {
    TEST("system_init");
    QCSystem sys;
    qc_system_init(&sys);
    CHECK(sys.state_dim == 0, "initial state_dim should be 0");
    CHECK(sys.stability == QC_STABLE_MARGINAL, "initial stability status");
    PASS();
}

/* ================================================================
 * Test 2: System Configuration
 * ================================================================ */
static void test_system_configure(void) {
    TEST("system_configure");
    QCSystem sys;
    qc_system_init(&sys);
    int ret = qc_system_configure(&sys, 2, 1, 1);
    CHECK(ret == 0, "configure should succeed");
    CHECK(sys.state_dim == 2, "state_dim=2");
    CHECK(sys.input_dim == 1, "input_dim=1");
    CHECK(sys.output_dim == 1, "output_dim=1");
    /* Check A initialized to identity */
    CHECK(fabs(sys.A[0] - 1.0) < 1e-9, "A[0,0]=1");
    CHECK(fabs(sys.A[3] - 1.0) < 1e-9, "A[1,1]=1");
    PASS();
}

/* ================================================================
 * Test 3: Uniform Quantizer Accuracy
 * ================================================================ */
static void test_uniform_quantizer(void) {
    TEST("uniform_quantizer");
    QCQuantizer q;
    qc_quantizer_init(&q, QC_QTYPE_UNIFORM, 8);
    qc_quantizer_configure_range(&q, -1.0, 1.0);

    /* Test mid-tread quantizer */
    double xq = qc_quantize_scalar(&q, 0.3);
    double err = fabs(0.3 - xq);
    CHECK(err <= q.step / 2.0 + 1e-9, "quantization error within bound");

    /* Test zero maps to zero */
    double xq0 = qc_quantize_scalar(&q, 0.0);
    CHECK(fabs(xq0) < 1e-9, "zero maps to zero");

    /* Test saturation */
    double xq_sat = qc_quantize_scalar(&q, 2.0);
    CHECK(xq_sat <= q.range_max + 1e-9, "saturation upper bound");

    qc_quantizer_free(&q);
    PASS();
}

/* ================================================================
 * Test 4: Quantization Noise Variance (L4)
 * ================================================================ */
static void test_quantization_noise(void) {
    TEST("quantization_noise");
    double step = 0.01;
    double expected_var = step * step / 12.0;
    double computed_var = qc_quantization_noise_variance(step);
    CHECK(fabs(computed_var - expected_var) < 1e-12,
          "noise variance = step^2/12");
    PASS();
}

/* ================================================================
 * Test 5: Logarithmic Quantizer Sector Bound (L4)
 * ================================================================ */
static void test_log_quantizer_sector(void) {
    TEST("log_quantizer_sector_bound");
    double rho = 0.5;
    double delta = qc_log_sector_delta(rho);
    double expected = (1.0 - rho) / (1.0 + rho);
    CHECK(fabs(delta - expected) < 1e-12,
          "delta = (1-rho)/(1+rho)");

    /* Test quantization */
    QCLogQuantizer lq;
    qc_log_quantizer_init(&lq, 0.5, 0.01);
    double xq = qc_log_quantize(&lq, 0.05);
    CHECK(xq >= 0.0, "positive input gives positive output");
    double xq_neg = qc_log_quantize(&lq, -0.05);
    CHECK(xq_neg <= 0.0, "negative input gives negative output");
    /* Deadzone test */
    double xq_zero = qc_log_quantize(&lq, 0.001);
    CHECK(fabs(xq_zero) < 1e-9, "deadzone maps to zero");

    qc_log_quantizer_free(&lq);
    PASS();
}

/* ================================================================
 * Test 6: Data Rate Theorem (L4) — minimum rate computation
 * ================================================================ */
static void test_data_rate_theorem(void) {
    TEST("data_rate_theorem");
    /* Unstable system: A = [[1.5, 0], [0, 0.5]]
     * Unstable eigenvalue at 1.5: requires log2(1.5) ≈ 0.585 bits/sample
     * Stable eigenvalue at 0.5: requires 0 bits
     * Total minimum: ≈ 0.585 */
    double A[4] = {1.5, 0.0, 0.0, 0.5};
    double min_rate = qc_data_rate_theoretical_min(A, 2);
    double expected = log2(1.5);
    CHECK(fabs(min_rate - expected) < 1e-6, "min_rate = log2(1.5)");

    /* Check stabilizability */
    int stab = qc_data_rate_is_stabilizable(A, 2, 1.0, 0.0);
    CHECK(stab == 1, "rate 1.0 > 0.585 → stabilizable");
    int not_stab = qc_data_rate_is_stabilizable(A, 2, 0.1, 0.0);
    CHECK(not_stab == 0, "rate 0.1 < 0.585 → not stabilizable");

    /* Stable system requires 0 rate */
    double A_stable[4] = {0.5, 0.0, 0.0, 0.5};
    double min_stable = qc_data_rate_theoretical_min(A_stable, 2);
    CHECK(fabs(min_stable) < 1e-9, "stable system requires 0 rate");

    PASS();
}

/* ================================================================
 * Test 7: Data Rate Allocation
 * ================================================================ */
static void test_data_rate_allocate(void) {
    TEST("data_rate_allocate");
    QCDataRate dr;
    qc_data_rate_init(&dr, 10.0, 2);
    double A[4] = {2.0, 0.0, 0.0, 1.5};
    CHECK(qc_data_rate_compute_min(&dr, A, 2) == 0, "compute_min OK");
    CHECK(dr.num_unstable == 2, "2 unstable eigenvalues");
    CHECK(qc_data_rate_check_stabilizability(&dr) == 1, "total > min");

    free(dr.channel_bits);
    free(dr.channel_rate);
    free(dr.eig_magnitudes);
    PASS();
}

/* ================================================================
 * Test 8: Encoder/Decoder Roundtrip
 * ================================================================ */
static void test_encoder_decoder(void) {
    TEST("encoder_decoder_roundtrip");
    QCQuantizer q;
    qc_quantizer_init(&q, QC_QTYPE_UNIFORM, 8);
    qc_quantizer_configure_range(&q, -1.0, 1.0);

    QCEncoder enc;
    QCDecoder dec;
    qc_encoder_init(&enc, QC_ENC_FIXED_LENGTH, 8);
    qc_decoder_init(&dec, QC_ENC_FIXED_LENGTH, 8);

    /* Encode a value */
    double original = 0.42;
    int bits = qc_encoder_encode_scalar(&enc, &q, original);
    CHECK(bits > 0, "encoding produces bits");

    /* Share buffer with decoder */
    dec.buffer = enc.buffer;
    dec.buffer_size = enc.buffer_size;
    dec.buffer_pos = 0;
    dec.bit_offset = 0;

    /* Decode */
    double decoded = qc_decoder_decode_scalar(&dec, &q);
    double xq = qc_quantize_scalar(&q, original);
    CHECK(fabs(decoded - xq) < q.step, "decoded ≈ quantized original");

    /* Prevent double-free */
    enc.buffer = NULL;
    qc_encoder_free(&enc);
    qc_decoder_free(&dec);
    qc_quantizer_free(&q);
    PASS();
}

/* ================================================================
 * Test 9: Huffman Coding
 * ================================================================ */
static void test_huffman_coding(void) {
    TEST("huffman_coding");
    double probs[4] = {0.4, 0.3, 0.2, 0.1};
    int symbols[4] = {0, 1, 2, 3};
    QCHuffmanTree tree;
    int ret = qc_huffman_build(&tree, probs, symbols, 4);
    CHECK(ret == 0, "huffman build OK");

    /* Check that entropy <= avg_code_length < entropy + 1 */
    double h = qc_entropy_discrete(probs, 4);
    CHECK(tree.avg_code_length >= h - 1e-9, "avg_len >= entropy");
    CHECK(tree.avg_code_length < h + 1.0, "avg_len < entropy + 1 (optimality)");
    CHECK(tree.efficiency > 0.8, "efficiency > 80%");

    /* Test encode/decode */
    QCBitWriter bw;
    qc_bit_writer_init(&bw, 256);
    int bits_enc = qc_huffman_encode(&tree, 2, &bw);
    CHECK(bits_enc > 0, "encoded bits > 0");
    qc_bit_writer_flush(&bw);

    QCBitReader br;
    qc_bit_reader_init(&br, bw.data, bw.byte_pos);
    int decoded_sym = -1;
    qc_huffman_decode(&tree, &br, &decoded_sym);
    CHECK(decoded_sym == 2, "huffman roundtrip OK");

    qc_bit_writer_free(&bw);
    qc_huffman_free(&tree);
    PASS();
}

/* ================================================================
 * Test 10: Dynamic Zoom Quantizer
 * ================================================================ */
static void test_dynamic_zoom(void) {
    TEST("dynamic_zoom_quantizer");
    QCDynamicQuantizer dq;
    qc_dyn_quantizer_init(&dq, 8, 1.0);

    double xq = qc_dyn_quantize(&dq, 0.3);
    CHECK(fabs(xq) <= dq.M * dq.mu + 1e-9, "quantized within range");

    /* Test zoom out */
    dq.mu = 0.1;
    int zoom_out = qc_dyn_should_zoom(&dq, 0.5);
    CHECK(zoom_out > 0, "should zoom out when |x| > M*mu");
    qc_dyn_zoom_out(&dq);
    CHECK(dq.mu > 0.1, "mu increased after zoom out");

    /* Test zoom in */
    dq.time_since_last_zoom = dq.dwell_time;
    int zoom_in = qc_dyn_should_zoom(&dq, 0.001);
    CHECK(zoom_in < 0, "should zoom in when |x| << M*mu");
    qc_dyn_zoom_in(&dq);
    CHECK(dq.mu < 0.2, "mu decreased after zoom in");

    PASS();
}

/* ================================================================
 * Test 11: Closed-Loop Simulation
 * ================================================================ */
static double test_controller(const double *y, int ny, double t, double *u) {
    (void)t;
    u[0] = -2.0 * y[0];
    return 0.0;
}

static void test_simulation(void) {
    TEST("closed_loop_simulation");
    QCSystem sys;
    qc_system_init(&sys);
    qc_system_configure(&sys, 1, 1, 1);

    /* Stable system: A = -1, B = 1, C = 1 */
    double A[1] = {-1.0};
    double B[1] = {1.0};
    double C[1] = {1.0};
    qc_system_set_A(&sys, A);
    qc_system_set_B(&sys, B);
    qc_system_set_C(&sys, C);

    double x0[1] = {1.0};
    QCSimulationResult res;
    qc_sim_result_init(&res);

    int ret = qc_simulate_closed_loop(&sys, x0, 0.0, 2.0, 0.01,
                                       test_controller, &res);
    CHECK(ret == 0, "simulation returned OK");
    CHECK(res.steps > 0, "steps > 0");

    /* System is stable, state should converge */
    CHECK(res.final_error < 1.0, "state converges (stable system)");

    qc_sim_result_free(&res);
    PASS();
}

/* ================================================================
 * Test 12: Quantized LQR Simulation
 * ================================================================ */
static void test_quantized_lqr(void) {
    TEST("quantized_lqr");
    QCSystem sys;
    qc_system_init(&sys);
    qc_system_configure(&sys, 1, 1, 1);

    double A[1] = {-0.5};
    double B[1] = {1.0};
    double C[1] = {1.0};
    qc_system_set_A(&sys, A);
    qc_system_set_B(&sys, B);
    qc_system_set_C(&sys, C);

    /* Simple LQR gain K = 0.5 */
    double K[1] = {0.5};
    double x0[1] = {1.0};
    QCSimulationResult res;
    qc_sim_result_init(&res);

    int ret = qc_simulate_quantized_lqr(&sys, x0, K, 0.0, 2.0, 0.01, &res);
    CHECK(ret == 0, "LQR simulation OK");
    CHECK(res.final_error < 1.0, "converges with quantized LQR");

    qc_sim_result_free(&res);
    PASS();
}

/* ================================================================
 * Test 13: Sector Bound Analysis
 * ================================================================ */
static void test_sector_bound(void) {
    TEST("sector_bound_analysis");
    QCSectorBoundResult sr;
    qc_sector_bound_init(&sr);

    QCSystem sys;
    qc_system_init(&sys);
    qc_system_configure(&sys, 1, 1, 1);
    double A[1] = {-1.0};
    double B[1] = {1.0};
    double C[1] = {1.0};
    qc_system_set_A(&sys, A);
    qc_system_set_B(&sys, B);
    qc_system_set_C(&sys, C);
    sys.sector_delta = 0.1;

    int ret = qc_sector_bound_analyze(&sys, &sr);
    CHECK(ret == 0, "sector bound analysis OK");
    CHECK(sr.nominally_stable == 1, "stable system detected");

    PASS();
}

/* ================================================================
 * Test 14: Matrix Operations
 * ================================================================ */
static void test_matrix_operations(void) {
    TEST("matrix_operations");
    double A[4] = {2.0, 1.0, 1.0, 2.0};
    double er[2], ei[2];
    int ret = qc_matrix_eigenvalues(A, 2, er, ei);
    CHECK(ret == 0, "eigenvalues computed");

    /* For symmetric matrix [[2,1],[1,2]], eigenvalues: 3 and 1 */
    double mag1 = sqrt(er[0]*er[0] + ei[0]*ei[0]);
    double mag2 = sqrt(er[1]*er[1] + ei[1]*ei[1]);
    CHECK(mag1 > 0.0 && mag2 > 0.0, "valid eigenvalue magnitudes");

    /* Spectral radius */
    double rho = qc_spectral_radius(A, 2);
    CHECK(rho >= 1.0, "spectral radius >= 1 for this matrix");

    PASS();
}

/* ================================================================
 * Test 15: Channel Capacity
 * ================================================================ */
static void test_channel_capacity(void) {
    TEST("channel_capacity");
    /* AWGN: C = B*log2(1+SNR) */
    double c_awgn = qc_channel_capacity_awgn(1.0, 1.0);
    CHECK(fabs(c_awgn - 1.0) < 1e-6, "AWGN C(1,1) = 1 bps");

    /* BSC(p=0): C = 1 */
    double c_bsc0 = qc_channel_capacity_bsc(0.0);
    CHECK(fabs(c_bsc0 - 1.0) < 1e-9, "BSC(0) capacity = 1");

    /* BSC(p=0.5): C = 0 */
    double c_bsc5 = qc_channel_capacity_bsc(0.5);
    CHECK(fabs(c_bsc5) < 1e-9, "BSC(0.5) capacity = 0");

    /* Erasure channel */
    double c_erase = qc_channel_capacity_erasure(0.2);
    CHECK(fabs(c_erase - 0.8) < 1e-9, "erasure(0.2) cap = 0.8");

    PASS();
}

/* ================================================================
 * Test 16: Lloyd-Max Quantizer
 * ================================================================ */
static void test_lloyd_max(void) {
    TEST("lloyd_max_gaussian");
    double levels[8], boundaries[7];
    double mse = qc_lloyd_max_gaussian(3, levels, boundaries, 50);
    CHECK(mse > 0.0, "Lloyd-Max returns positive MSE");
    /* Levels should be sorted */
    for (int i = 0; i < 7; i++) {
        CHECK(levels[i] < levels[i+1], "levels sorted ascending");
    }
    PASS();
}

/* ================================================================
 * Test 17: Vector Quantizer
 * ================================================================ */
static void test_vector_quantizer(void) {
    TEST("vector_quantizer");
    QCVectorCodebook cb;
    int ret = qc_vec_codebook_init(&cb, 4, 2);
    CHECK(ret == 0, "codebook init OK");
    CHECK(cb.num_codewords == 4, "4 codewords");
    CHECK(cb.dimension == 2, "dim=2");

    double x[2] = {0.5, -0.3};
    int idx;
    double dist;
    ret = qc_vec_quantize(&cb, x, &idx, &dist);
    CHECK(ret == 0, "vector quantize OK");
    CHECK(idx >= 0 && idx < 4, "valid index");

    double xr[2];
    ret = qc_vec_reconstruct(&cb, idx, xr);
    CHECK(ret == 0, "reconstruct OK");

    qc_vec_codebook_free(&cb);
    PASS();
}

/* ================================================================
 * Test 18: Delta Modulation
 * ================================================================ */
static void test_delta_modulation(void) {
    TEST("delta_modulation");
    QCDeltaModulator dm;
    qc_delta_mod_init(&dm, 0.1, 0.01, 1.0);

    int b1 = qc_delta_mod_encode(&dm, 0.15);
    CHECK(b1 == 1, "increase detection");

    double dec1 = qc_delta_mod_decode(&dm, b1);
    CHECK(dec1 > 0.0, "decoded positive");

    int b2 = qc_delta_mod_encode(&dm, -0.1);
    CHECK(b2 == 0, "decrease detection");

    PASS();
}

/* ================================================================
 * Test 19: SNR and ENOB
 * ================================================================ */
static void test_snr_enob(void) {
    TEST("snr_enob");
    double enob = qc_enob_from_snr(61.96);
    CHECK(fabs(enob - 10.0) < 0.1, "ENOB ≈ 10 for 62dB SNR");

    double snr = qc_snr_from_enob(10.0);
    CHECK(fabs(snr - 61.96) < 0.1, "SNR ≈ 62dB for 10 ENOB");

    PASS();
}

/* ================================================================
 * Test 20: Zoom Verification
 * ================================================================ */
static void test_zoom_verification(void) {
    TEST("zoom_verification");
    QCZoomVerification zv;
    qc_zoom_verify_init(&zv);

    /* Simulate decreasing Lyapunov function (20 steps for clear convergence) */
    for (int i = 0; i < 20; i++) {
        double V = exp(-0.5 * i);
        double V_prev = (i > 0) ? exp(-0.5 * (i - 1)) : V;
        double norm = exp(-0.25 * i);
        qc_zoom_verify_step(&zv, V, V_prev, norm);
    }

    int converged = qc_zoom_verify_conclusion(&zv);
    CHECK(zv.bounded == 1, "state bounded");
    CHECK(zv.converges_to_zero == 1, "converges to zero");
    CHECK(converged == 1, "zoom strategy successful");

    qc_zoom_verify_free(&zv);
    PASS();
}

/* ================================================================
 * Main
 * ================================================================ */
int main(void) {
    printf("========================================\n");
    printf(" mini-quantized-control Test Suite\n");
    printf("========================================\n");

    test_system_init();
    test_system_configure();
    test_uniform_quantizer();
    test_quantization_noise();
    test_log_quantizer_sector();
    test_data_rate_theorem();
    test_data_rate_allocate();
    test_encoder_decoder();
    test_huffman_coding();
    test_dynamic_zoom();
    test_simulation();
    test_quantized_lqr();
    test_sector_bound();
    test_matrix_operations();
    test_channel_capacity();
    test_lloyd_max();
    test_vector_quantizer();
    test_delta_modulation();
    test_snr_enob();
    test_zoom_verification();

    printf("========================================\n");
    printf(" Results: %d/%d tests PASSED\n", tests_passed, tests_run);
    printf("========================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
