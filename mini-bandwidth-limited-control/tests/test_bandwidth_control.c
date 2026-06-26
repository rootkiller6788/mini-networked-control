/**
 * test_bandwidth_control.c — Comprehensive Test Suite
 *
 * Tests all core APIs for the bandwidth-limited control module.
 * Covers L1-L6 with mathematical assertions.
 */

#include "blc_core.h"
#include "blc_datarate.h"
#include "blc_encoding.h"
#include "blc_control.h"
#include "blc_event.h"
#include "blc_scheduling.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %-50s", name); \
} while(0)

#define PASS() do { \
    tests_passed++; \
    printf(" PASS\n"); \
} while(0)

#define FAIL(msg) do { \
    printf(" FAIL: %s\n", msg); \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { FAIL(msg); return; } \
} while(0)

#define ASSERT_FLOAT_EQ(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        char buf[128]; \
        snprintf(buf, sizeof(buf), "expected %g, got %g", (double)(b), (double)(a)); \
        FAIL(buf); return; \
    } \
} while(0)

/* ================================================================
 * L1: Core Definitions Tests
 * ================================================================ */
static void test_system_create(void) {
    TEST("blc_create/free lifecycle");
    BLCSystem* sys = blc_create(2, 1, 1);
    ASSERT_TRUE(sys != NULL, "create returned NULL");
    ASSERT_TRUE(sys->plant.n_states == 2, "n_states mismatch");
    ASSERT_TRUE(sys->sample_period > 0, "sample_period not set");
    blc_free(sys);
    PASS();
}

static void test_system_create_null_on_bad_args(void) {
    TEST("blc_create rejects bad args");
    BLCSystem* s1 = blc_create(0, 1, 1);
    ASSERT_TRUE(s1 == NULL, "should reject n_states=0");
    BLCSystem* s2 = blc_create(100, 1, 1);
    ASSERT_TRUE(s2 == NULL, "should reject too many states");
    BLCSystem* s3 = blc_create(2, -1, 1);
    ASSERT_TRUE(s3 == NULL, "should reject negative inputs");
    PASS();
}

static void test_channel_init(void) {
    TEST("blc_init_channel");
    BLCSystem* sys = blc_create(2, 1, 1);
    int rc = blc_init_channel(sys, 1000.0, 10.0, 5.0, 0.01);
    ASSERT_TRUE(rc == 0, "channel init failed");
    ASSERT_FLOAT_EQ(sys->channel.bandwidth_hz, 1000.0, 1e-6);
    ASSERT_TRUE(sys->channel.capacity_bps > 0, "capacity should be positive");
    blc_free(sys);
    PASS();
}

static void test_quantizer_uniform(void) {
    TEST("Uniform quantizer operations");
    BLCQuantizer q;
    int rc = blc_quantizer_init(&q, 256, -10.0, 10.0, false);
    ASSERT_TRUE(rc == 0, "quantizer init failed");
    ASSERT_FLOAT_EQ(q.step, 20.0/256.0, 0.01);
    ASSERT_FLOAT_EQ(q.max_error, q.step/2.0, 1e-6);

    double qv = blc_quantize(&q, 5.0);
    ASSERT_TRUE(fabs(qv - 5.0) < q.max_error * 1.1, "quantization error too large");

    double err0 = blc_quantization_error(&q, 0.0);
    ASSERT_TRUE(fabs(err0) <= q.max_error * 1.1, "zero quantization error");

    int idx = blc_quantize_to_index(&q, 0.0);
    double dq = blc_dequantize(&q, idx);
    ASSERT_TRUE(fabs(dq - 0.0) < 0.1, "dequantize roundtrip failed");
    PASS();
}

static void test_quantizer_overload(void) {
    TEST("Quantizer overload detection");
    BLCQuantizer q;
    blc_quantizer_init(&q, 256, -10.0, 10.0, false);

    double step = 20.0 / 256.0;
    double qv = blc_quantize(&q, 100.0);
    /** Should be near upper bound (within one quantization step) */
    ASSERT_TRUE(fabs(qv - 10.0) < step * 1.5, "should saturate near upper bound");
    ASSERT_TRUE(q.overload_count >= 1, "overload not counted");

    qv = blc_quantize(&q, -100.0);
    ASSERT_TRUE(fabs(qv + 10.0) < step * 1.5, "should saturate near lower bound");
    ASSERT_TRUE(q.overload_count >= 2, "overload not counted");
    PASS();
}

/* ================================================================
 * L2: Core Concepts Tests
 * ================================================================ */
static void test_plant_init(void) {
    TEST("blc_init_plant with 2-state system");
    BLCSystem* sys = blc_create(2, 1, 1);

    double A[4] = {0.0, 1.0, -2.0, -3.0};
    double B[2] = {0.0, 1.0};
    double C[2] = {1.0, 0.0};

    int rc = blc_init_plant(sys, A, B, C);
    ASSERT_TRUE(rc == 0, "plant init failed");
    ASSERT_FLOAT_EQ(sys->plant.A[0][0], 0.0, 1e-6);
    ASSERT_FLOAT_EQ(sys->plant.A[0][1], 1.0, 1e-6);
    ASSERT_FLOAT_EQ(sys->plant.A[1][0], -2.0, 1e-6);

    blc_free(sys);
    PASS();
}

static void test_channel_capacity(void) {
    TEST("Shannon-Hartley channel capacity");
    BLCSystem* sys = blc_create(2, 1, 1);
    blc_init_channel(sys, 100.0, 15.0, 1.0, 0.0);

    double cap = blc_channel_capacity(&sys->channel);
    double expected = 100.0 * log2(16.0);
    ASSERT_FLOAT_EQ(cap, expected, 0.01);
    ASSERT_TRUE(cap > 0.0, "capacity must be positive");
    blc_free(sys);
    PASS();
}

static void test_minimum_datarate(void) {
    TEST("Data Rate Theorem — minimum datarate");
    BLCSystem* sys = blc_create(2, 1, 1);

    /** Set unstable eigenvalues manually */
    sys->plant.eigenvalues[0] = 1.5;     /** |λ| = 1.5, contributes log₂(1.5) */
    sys->plant.eigenvalues[1] = 0.5;     /** |λ| = 0.5, stable, contributes 0 */
    sys->plant.eigenvalues_im[0] = 0.0;
    sys->plant.eigenvalues_im[1] = 0.0;

    double rate = blc_minimum_datarate(&sys->plant);
    double expected = log2(1.5);
    ASSERT_FLOAT_EQ(rate, expected, 0.01);

    blc_free(sys);
    PASS();
}

static void test_packet_operations(void) {
    TEST("Packet encode/decode roundtrip");
    BLCSystem* sys = blc_create(2, 1, 1);
    blc_quantizer_init(&sys->sensor_quant, 256, -10.0, 10.0, false);

    double x0[2] = {1.5, -2.3};
    blc_set_initial_state(sys, x0);

    /** Encode state into packet */
    BLCPacket* pkt = blc_packet_create(64, 0, 0.0, true);
    int bits = blc_packet_encode_state(sys, pkt);
    ASSERT_TRUE(bits > 0, "no bits encoded");

    /** Modify state to verify decode */
    sys->plant.x_hat[0] = 99.0;
    sys->plant.x_hat[1] = 99.0;

    int decoded = blc_packet_decode_state(sys, pkt);
    ASSERT_TRUE(decoded > 0, "decode returned 0");

    /** Check that estimate approximates true state */
    double e0 = fabs(sys->plant.x_hat[0] - x0[0]);
    double e1 = fabs(sys->plant.x_hat[1] - x0[1]);
    ASSERT_TRUE(e0 < sys->sensor_quant.max_error * 2, "decode error too large");
    ASSERT_TRUE(e1 < sys->sensor_quant.max_error * 2, "decode error too large");

    blc_packet_free(pkt);
    blc_free(sys);
    PASS();
}

/* ================================================================
 * L3: Mathematical Structures Tests
 * ================================================================ */
static void test_simulation_step(void) {
    TEST("Simulation step with quantized feedback");
    BLCSystem* sys = blc_create(2, 1, 1);

    /** Double integrator: A = [[0,1],[0,0]], B = [[0],[1]] */
    double A[4] = {0.0, 1.0, 0.0, 0.0};
    double B[2] = {0.0, 1.0};
    double C[2] = {1.0, 0.0};
    blc_init_plant(sys, A, B, C);

    double x0[2] = {1.0, 0.0};
    blc_set_initial_state(sys, x0);

    /** LQR gain for double integrator */
    double K_flat[2] = {1.0, 1.414};
    blc_set_controller_gain(sys, K_flat);

    /** Run a few simulation steps */
    for (int i = 0; i < 50; i++) {
        blc_simulate_step(sys);
    }

    /** Check that state is bounded (not diverged to infinity) */
    double x[2];
    blc_get_state(sys, x);
    ASSERT_TRUE(fabs(x[0]) < 1e6, "state diverged");
    ASSERT_TRUE(fabs(x[1]) < 1e6, "state diverged");
    ASSERT_TRUE(sys->packets_sent > 0, "no packets sent");
    blc_free(sys);
    PASS();
}

static void test_eigenvalue_computation(void) {
    TEST("QR eigenvalue computation");
    /** Diagonal matrix */
    double A[9] = {2.0, 0.0, 0.0,  0.0, -1.0, 0.0,  0.0, 0.0, 3.0};
    double re[3], im[3];
    int iter = blc_eigenvalues(A, 3, re, im);
    ASSERT_TRUE(iter >= 0, "eigenvalue computation failed");
    /** Check that eigenvalues are near {2, -1, 3} */
    double found[3] = {re[0], re[1], re[2]};
    double expected[3] = {2.0, -1.0, 3.0};
    int matched = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (fabs(found[i] - expected[j]) < 0.1) matched++;
        }
    }
    ASSERT_TRUE(matched >= 3, "not all eigenvalues matched");
    PASS();
}

static void test_spectral_radius(void) {
    TEST("Spectral radius computation");
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double rho = blc_spectral_radius(A, 2);
    /** ρ ≈ 5.372 for this matrix */
    ASSERT_TRUE(rho > 5.0, "spectral radius too small");
    ASSERT_TRUE(rho < 6.0, "spectral radius too large");
    PASS();
}

/* ================================================================
 * L4: Fundamental Theorems Tests
 * ================================================================ */
static void test_datarate_theorem_ct(void) {
    TEST("Data Rate Theorem — continuous-time");
    BLCSystem* sys = blc_create(2, 1, 1);

    /** Unstable eigenvalue λ = 2.0 */
    sys->plant.eigenvalues[0] = 2.0;
    sys->plant.eigenvalues[1] = -1.0;
    sys->plant.eigenvalues_im[0] = 0.0;
    sys->plant.eigenvalues_im[1] = 0.0;

    double R_ct = blc_datarate_min_ct(&sys->plant);
    /** R = 2 * 2.0 / ln(2) ≈ 5.77 bps */
    double expected = 4.0 / log(2.0);
    ASSERT_FLOAT_EQ(R_ct, expected, 0.01);

    /** Check sufficiency */
    ASSERT_TRUE(blc_datarate_is_sufficient(10.0, &sys->plant, 0.01),
                "10 bps should be sufficient for λ=2");
    ASSERT_TRUE(!blc_datarate_is_sufficient(1.0, &sys->plant, 0.01),
                "1 bps should not be sufficient for λ=2");

    blc_free(sys);
    PASS();
}

static void test_rate_distortion(void) {
    TEST("Rate-Distortion function for control");
    BLCSystem* sys = blc_create(2, 1, 1);
    sys->plant.eigenvalues[0] = 1.0;  /** λ = 1 */
    sys->plant.eigenvalues[1] = -2.0;
    sys->plant.eigenvalues_im[0] = 0.0;
    sys->plant.eigenvalues_im[1] = 0.0;

    /** R(D) should decrease as D increases */
    double R1 = blc_rate_distortion(&sys->plant, 0.1, 1.0);
    double R2 = blc_rate_distortion(&sys->plant, 1.0, 1.0);
    ASSERT_TRUE(R1 > R2, "rate should decrease with higher distortion");

    /** D(R) should decrease as R increases */
    double D1 = blc_distortion_rate(&sys->plant, 5.0, 1.0);
    double D2 = blc_distortion_rate(&sys->plant, 10.0, 1.0);
    ASSERT_TRUE(D2 < D1, "distortion should decrease with higher rate");

    blc_free(sys);
    PASS();
}

/* ================================================================
 * L5: Algorithms Tests
 * ================================================================ */
static void test_zoom_quantizer(void) {
    TEST("Zoom quantizer — dynamic range adjustment");
    BLCZoomQuantizer zq;
    int rc = blc_zoom_init(&zq, 10.0, 2.0, 2.0, 256);
    ASSERT_TRUE(rc == 0, "zoom init failed");

    /** Encode a small value — should trigger zoom-in */
    int idx = blc_zoom_encode(&zq, 1.0);
    ASSERT_TRUE(idx >= 0 && idx < 256, "encode returned bad index");

    /** Check that range decreased (zoom-in) since 1.0 < 10.0/2.0 */
    ASSERT_TRUE(zq.range < 10.0, "zoom-in should reduce range");

    /** Encode a large value — should trigger zoom-out */
    double prev_range = zq.range;
    idx = blc_zoom_encode(&zq, zq.range * 2.0);
    ASSERT_TRUE(zq.range > prev_range, "zoom-out should increase range");

    /** Decode */
    double val = blc_zoom_decode(&zq, 128);
    ASSERT_TRUE(fabs(val) < zq.range + 1e-6, "decoded value outside range");

    int zi, zo;
    double cr;
    blc_zoom_stats(&zq, &zi, &zo, &cr);
    ASSERT_TRUE(zi > 0, "should have zoomed in");
    ASSERT_TRUE(zo > 0, "should have zoomed out");
    PASS();
}

static void test_delta_modulation(void) {
    TEST("Delta modulation — 1-bit encoding");
    BLCDeltaModulator dm;
    blc_delta_init(&dm, 0.1, 0.0);

    /** Encode a fast ramp (increment > step size to avoid overshoot oscillation) */
    int ones = 0;
    for (int i = 0; i < 20; i++) {
        int bit = blc_delta_encode(&dm, (double)i * 0.2);
        ASSERT_TRUE(bit == 0 || bit == 1, "bad bit value");
        if (bit == 1) ones++;
    }
    /** With fast ramp, most bits should be 1 */
    ASSERT_TRUE(ones >= 10, "positive ramp should produce mostly 1s");

    double recon = blc_delta_get_value(&dm);
    ASSERT_TRUE(recon > 0.0, "reconstruction should be positive");

    PASS();
}

static void test_logarithmic_quantizer(void) {
    TEST("Logarithmic quantizer — optimal for stabilization");
    BLCLogQuantizer lq;
    int rc = blc_log_quant_init(&lq, 2.0, 4, 0.1);
    ASSERT_TRUE(rc == 0, "log quant init failed");

    /** Check that the quantizer was initialized with a valid table */
    ASSERT_TRUE(lq.level_table != NULL, "level table not allocated");
    ASSERT_TRUE(lq.recon_table != NULL, "recon table not allocated");

    /** Encode/decode symmetry: positive and negative should map to symmetric indices */
    int idx_pos = blc_log_quant_encode(&lq, 1.0);
    int idx_neg = blc_log_quant_encode(&lq, -1.0);
    ASSERT_TRUE(idx_pos != idx_neg, "positive and negative should differ");

    double dec_pos = blc_log_quant_decode(&lq, idx_pos);
    double dec_neg = blc_log_quant_decode(&lq, idx_neg);
    ASSERT_TRUE(dec_pos > 0, "positive decoded should be positive");
    ASSERT_TRUE(dec_neg < 0, "negative decoded should be negative");

    /** Deadzone: value within deadzone region (0.05 < 0.1) should map to zero-ish */
    int idx_dz = blc_log_quant_encode(&lq, 0.05);
    /** The deadzone maps to the middle index, reconstruction should be near 0 */
    ASSERT_TRUE(idx_dz == 4, "deadzone should map to middle index (4 for N=4)");

    /** Test optimal density formula */
    double dens = blc_log_quant_optimal_density(1.5, 0.01);
    ASSERT_FLOAT_EQ(dens, exp(1.5 * 0.01), 0.001);

    int min_levels = blc_log_quant_min_levels(dens);
    ASSERT_TRUE(min_levels >= 1, "min_levels should be positive");

    /** Encode overload: value far beyond range saturates */
    int idx_hi = blc_log_quant_encode(&lq, 100.0);
    int total_levels = lq.num_levels_pos * 2 + 1;
    ASSERT_TRUE(idx_hi >= 0 && idx_hi < total_levels, "overload index out of range");

    blc_log_quant_free(&lq);
    PASS();
}

static void test_predictive_encoder(void) {
    TEST("Predictive encoder — state estimation");
    BLCEncoder enc;
    int rc = blc_encoder_init(&enc, NULL, 5);
    ASSERT_TRUE(rc == 0, "encoder init failed");

    /** After a few encode/decode cycles with zero input,
     *  the estimate should remain near the predictions */
    double true_state[2] = {0.5, -0.3};

    /** Setup dummy plant */
    BLCPlant plant;
    memset(&plant, 0, sizeof(plant));
    plant.n_states = 2;
    plant.n_inputs = 1;
    plant.A[0][0] = 0.99; plant.A[0][1] = 0.0;
    plant.A[1][0] = 0.0;  plant.A[1][1] = 0.99;

    BLCQuantizer q;
    blc_quantizer_init(&q, 256, -10.0, 10.0, false);

    int bits = 0;
    BLCPacket pkt;

    blc_encoder_encode(&enc, &q, true_state, &plant, 0.01, &pkt, &bits);
    ASSERT_TRUE(bits > 0, "no bits generated");
    ASSERT_TRUE(enc.error_norm >= 0.0, "error norm should be non-negative");

    double x_hat[2];
    blc_encoder_get_estimate(&enc, x_hat);
    ASSERT_TRUE(fabs(x_hat[0]) < 10.0, "estimate out of range");

    PASS();
}

static void test_lqr_controller(void) {
    TEST("LQR controller — gain and Riccati access");
    BLCLQRController lqr;
    double Q_diag[2] = {1.0, 1.0};
    double R_diag[1] = {1.0};

    int rc = blc_lqr_init(&lqr, 2, 1, Q_diag, R_diag);
    ASSERT_TRUE(rc == 0, "LQR init failed");

    /** Verify initialization */
    ASSERT_TRUE(lqr.n_states == 2, "n_states wrong");
    ASSERT_TRUE(lqr.n_inputs == 1, "n_inputs wrong");

    /** Manually set gains (simulating a solved LQR) */
    lqr.K[0][0] = 2.0;
    lqr.K[0][1] = 1.5;
    lqr.P[0][0] = 3.0;
    lqr.P[0][1] = 1.0;
    lqr.P[1][0] = 1.0;
    lqr.P[1][1] = 2.0;

    /** Test gain extraction */
    double K_flat[2];
    blc_lqr_get_gain(&lqr, K_flat);
    ASSERT_FLOAT_EQ(K_flat[0], 2.0, 1e-6);
    ASSERT_FLOAT_EQ(K_flat[1], 1.5, 1e-6);

    /** Test Riccati extraction */
    double P_flat[4];
    blc_lqr_get_riccati(&lqr, P_flat);
    ASSERT_FLOAT_EQ(P_flat[0], 3.0, 1e-6);
    ASSERT_FLOAT_EQ(P_flat[3], 2.0, 1e-6);

    /** Test quantization penalty */
    double A[4] = {-1.0, 1.0, 0.0, -2.0};
    double B[2] = {0.0, 1.0};
    double quant_errors[2] = {0.1, 0.1};
    double penalty = blc_lqr_quantization_penalty(&lqr, A, B, quant_errors);
    ASSERT_TRUE(penalty >= 0.0, "penalty should be non-negative");

    /** Test quantized state feedback */
    BLCSystem* sys = blc_create(2, 1, 1);
    blc_init_plant(sys, A, B, A);  /** C = A for testing */
    double x[2] = {1.0, 0.5};
    double u[1], qe[2];
    blc_quantized_state_feedback(sys, &lqr, x, u, qe);
    ASSERT_TRUE(u[0] != 0.0, "control should be non-zero");

    blc_free(sys);
    blc_lqr_free(&lqr);
    PASS();
}

static void test_huffman_coding(void) {
    TEST("Huffman coding — entropy coding");
    BLCHuffmanCoder huff;
    huff.alphabet_size = 4;
    double probs[4] = {0.4, 0.3, 0.2, 0.1};

    int rc = blc_huffman_build(&huff, probs);
    ASSERT_TRUE(rc == 0, "Huffman build failed");

    /** Entropy H = -0.4*log2(0.4) - 0.3*log2(0.3) - 0.2*log2(0.2) - 0.1*log2(0.1) ≈ 1.846 */
    ASSERT_TRUE(huff.entropy > 1.5, "entropy too low");
    ASSERT_TRUE(huff.entropy < 2.0, "entropy too high");

    /** Average code length should be between H and H+1 */
    ASSERT_TRUE(huff.avg_code_length >= huff.entropy, "avg length < entropy");
    ASSERT_TRUE(huff.avg_code_length <= huff.entropy + 1.0 + 0.01,
                "avg length > H+1");

    ASSERT_TRUE(huff.is_optimal, "Huffman should be optimal");

    blc_huffman_free(&huff);
    PASS();
}

/* ================================================================
 * L6: Canonical Problems Tests
 * ================================================================ */
static void test_event_triggered(void) {
    TEST("Send-on-Delta event-triggered transmission");
    BLCSendOnDelta sod;
    blc_sod_init(&sod, 0.5, 2);  /** L2 norm, δ=0.5 */

    double x[2] = {0.0, 0.0};

    /** First call should transmit */
    bool tx = blc_sod_should_transmit(&sod, x, 2, 0.0, true);
    ASSERT_TRUE(tx, "first forced TX should succeed");
    blc_sod_transmitted(&sod, x, 0.0);

    /** Small change — should not transmit */
    double x2[2] = {0.1, 0.0};
    tx = blc_sod_should_transmit(&sod, x2, 2, 0.1, false);
    ASSERT_TRUE(!tx, "small change should not trigger");

    /** Large change — should transmit */
    double x3[2] = {1.0, 0.8};
    tx = blc_sod_should_transmit(&sod, x3, 2, 0.2, false);
    ASSERT_TRUE(tx, "large change should trigger");

    double pct_saved;
    blc_sod_stats(&sod, &pct_saved, NULL, NULL);
    ASSERT_TRUE(pct_saved >= 0.0, "bandwidth saved should be non-negative");
    PASS();
}

static void test_tdma_scheduling(void) {
    TEST("TDMA schedule for multi-loop control");
    BLCPlant plant;
    memset(&plant, 0, sizeof(plant));
    plant.n_states = 2;
    plant.eigenvalues[0] = 1.5;
    plant.eigenvalues[1] = 0.5;

    BLCLoopDescriptor loops[3];
    blc_loop_init(&loops[0], 0, &plant);
    blc_loop_init(&loops[1], 1, &plant);
    blc_loop_init(&loops[2], 2, &plant);

    BLCTDMASchedule tdma;
    int rc = blc_tdma_build(&tdma, loops, 3, 0.1, 0.001);
    ASSERT_TRUE(rc == 0, "TDMA build failed");
    ASSERT_TRUE(tdma.n_slots > 0, "no slots");
    ASSERT_TRUE(tdma.slot_duration > 0, "zero slot duration");

    /** Walk through slots */
    int owner = blc_tdma_current_owner(&tdma);
    ASSERT_TRUE(owner >= 0, "no owner");

    int next = blc_tdma_next_slot(&tdma, 0.01);
    ASSERT_TRUE(next >= 0 && next < 3, "invalid next slot");

    PASS();
}

static void test_bandwidth_allocation(void) {
    TEST("Bandwidth allocation — proportional and max-min");
    BLCPlant plant;
    memset(&plant, 0, sizeof(plant));
    plant.n_states = 2;
    plant.eigenvalues[0] = 2.0;  /** Priority ~2.0 */
    plant.eigenvalues[1] = 0.5;

    BLCLoopDescriptor loops[3];
    blc_loop_init(&loops[0], 0, &plant);

    plant.eigenvalues[0] = 1.0;  /** Priority ~1.0 */
    blc_loop_init(&loops[1], 1, &plant);

    plant.eigenvalues[0] = 3.0;  /** Priority ~3.0 */
    blc_loop_init(&loops[2], 2, &plant);

    BLCBandwidthAllocator alloc;
    blc_alloc_init(&alloc, 1000.0, loops, 3);

    /** Proportional allocation */
    blc_alloc_proportional(&alloc);
    ASSERT_TRUE(alloc.stability_margin >= 0.0, "negative margin");
    /** Loop 2 (priority 3.0) should get more than loop 1 (priority 1.0) */
    ASSERT_TRUE(alloc.loops[2].allocated_rate > alloc.loops[1].allocated_rate,
                "higher priority should get more bandwidth");

    blc_alloc_print(&alloc, stdout);
    PASS();
}

static void test_tod_protocol(void) {
    TEST("Try-Once-Discard protocol");
    BLCTODProtocol tod;
    blc_tod_init(&tod, 3);
    double weights[3] = {1.0, 2.0, 3.0};
    blc_tod_set_weights(&tod, weights);

    /** Node 1 has larger weight but same error as node 2.
     *  Node 2 has larger error. Winner should be argmax(w_i * e_i) */
    double errors[3] = {0.5, 1.0, 0.2};
    int winner = blc_tod_resolve(&tod, errors, 3);
    /** w*e: node0=0.5, node1=2.0, node2=0.6 → node1 wins */
    ASSERT_TRUE(winner == 1, "TOD should pick node with max weight*error");

    blc_tod_transmitted(&tod, winner);
    int wins[3], discards[3], cont;
    blc_tod_stats(&tod, wins, discards, &cont);
    ASSERT_TRUE(wins[1] > 0, "winner not recorded");

    PASS();
}

static void test_adaptive_scheduling(void) {
    TEST("Adaptive bandwidth scheduling");
    BLCAdaptiveScheduler sched;
    blc_adaptive_init(&sched, 3, 1000.0, 1.0, 0.1);

    BLCPlant plant;
    memset(&plant, 0, sizeof(plant));
    plant.n_states = 2;
    plant.eigenvalues[0] = 1.5;

    blc_adaptive_register_loop(&sched, 0, &plant, 50.0, 500.0);
    blc_adaptive_register_loop(&sched, 1, &plant, 50.0, 500.0);
    blc_adaptive_register_loop(&sched, 2, &plant, 50.0, 500.0);

    double errors[3] = {0.1, 0.5, 0.2};
    double cost = blc_adaptive_epoch(&sched, errors);
    ASSERT_TRUE(cost > 0.0, "cost should be positive");

    /** Each loop should have >= min_rate */
    for (int i = 0; i < 3; i++) {
        double rate = blc_adaptive_get_rate(&sched, i);
        ASSERT_TRUE(rate >= 50.0, "below min rate");
        ASSERT_TRUE(rate <= 500.0, "above max rate");
    }

    blc_adaptive_print(&sched, stdout);
    PASS();
}

/* ================================================================
 * L7/L8: Advanced Tests (Lebesgue, Self-triggered, Vector Quant)
 * ================================================================ */
static void test_lebesgue_sampler(void) {
    TEST("Lebesgue sampler — level-crossing detection");
    BLCLebesgueSampler ls;
    int rc = blc_lebesgue_init(&ls, exp(1.0), 8, 5.0);
    ASSERT_TRUE(rc == 0, "Lebesgue init failed");

    /** Signal crossing levels */
    double signal[] = {0.0, 0.5, 1.5, 0.2, 3.0};
    int crossings = 0;
    for (int i = 0; i < 5; i++) {
        if (blc_lebesgue_check(&ls, signal[i], (double)i * 0.1)) {
            crossings++;
        }
    }
    ASSERT_TRUE(crossings > 0, "should detect level crossings");

    double avg_interval;
    blc_lebesgue_stats(&ls, &avg_interval, NULL);
    PASS();
}

static void test_event_trig_feedback(void) {
    TEST("Event-triggered feedback (Tabuada 2007)");
    BLCEventTriggeredFeedback etf;

    /** Stable closed-loop: A_c = -2 */
    double Ac[4] = {-2.0, 0.0, 0.0, -2.0};
    int rc = blc_etf_init(&etf, 0.1, Ac, 2);
    ASSERT_TRUE(rc == 0, "ETF init failed");
    ASSERT_TRUE(etf.min_inter_event > 0, "min inter-event should be positive");

    /** First check should trigger (x_last is zero) */
    double x[2] = {1.0, 0.5};
    bool trigger = blc_etf_check(&etf, x, 2, 0.0);
    ASSERT_TRUE(trigger, "initial check should trigger");

    blc_etf_triggered(&etf, x, 0.0);

    /** Small change should not trigger */
    double x2[2] = {1.01, 0.5};
    trigger = blc_etf_check(&etf, x2, 2, 0.01);
    ASSERT_TRUE(!trigger, "small change should not trigger");

    int ev_count;
    double avg_int;
    blc_etf_stats(&etf, &ev_count, &avg_int);
    ASSERT_TRUE(ev_count > 0, "should have counted events");

    PASS();
}

static void test_self_triggered(void) {
    TEST("Self-triggered control — next event pre-computation");
    BLCSelfTriggered st;
    double Ac[4] = {-1.0, 0.0, 0.0, -1.0};
    double P[4]  = {1.0, 0.0, 0.0, 1.0};

    int rc = blc_self_trig_init(&st, 0.5, 0.01, 1.0, Ac, P, 2);
    ASSERT_TRUE(rc == 0, "self-trig init failed");

    double x[2] = {2.0, 0.0};
    double next = blc_self_trig_next_event(&st, x, 2, 0.0);
    ASSERT_TRUE(next > 0.01, "next event too soon");
    ASSERT_TRUE(next <= 1.0, "next event beyond max");

    PASS();
}

static void test_runlength_encoding(void) {
    TEST("Run-length encoding for sparse control");
    double signal[20] = {0,0,0,0,0, 1.0,1.0,1.0, 0,0,0,0,0, -0.5,-0.5, 0,0,0,0,0};
    double run_vals[10];
    int run_lens[10];

    int runs = blc_runlength_encode(signal, 20, run_vals, run_lens, 10);
    ASSERT_TRUE(runs > 0, "should find runs");
    ASSERT_TRUE(runs <= 5, "too many runs for this signal");
    /** Total count should sum to 20 */
    int total = 0;
    for (int i = 0; i < runs; i++) total += run_lens[i];
    ASSERT_TRUE(total == 20, "run lengths don't sum to total");

    PASS();
}

static void test_bit_allocation(void) {
    TEST("Eigenvalue-based bit allocation");
    double eigenvalues[3] = {2.0, 0.5, 1.5};  /** States 0 and 2 are unstable */
    int allocation[3];
    int total_bits = 16;

    blc_allocate_bits_eigenvalue(eigenvalues, 3, total_bits, allocation);

    /** Unstable modes should get more bits */
    ASSERT_TRUE(allocation[0] > allocation[1],
                "unstable mode should get more bits than stable");
    ASSERT_TRUE(allocation[2] > allocation[1],
                "unstable mode should get more bits than stable");

    /** Sum should approximately equal total_bits */
    int sum = allocation[0] + allocation[1] + allocation[2];
    ASSERT_TRUE(sum >= total_bits - 2, "bit sum too low");

    PASS();
}

/* ================================================================
 * Main
 * ================================================================ */
int main(void) {
    printf("\n=== Bandwidth-Limited Control Test Suite ===\n\n");

    /** L1: Definitions */
    printf("--- L1: Core Definitions ---\n");
    test_system_create();
    test_system_create_null_on_bad_args();
    test_channel_init();
    test_quantizer_uniform();
    test_quantizer_overload();

    /** L2: Core Concepts */
    printf("\n--- L2: Core Concepts ---\n");
    test_plant_init();
    test_channel_capacity();
    test_minimum_datarate();
    test_packet_operations();

    /** L3: Mathematical Structures */
    printf("\n--- L3: Mathematical Structures ---\n");
    test_simulation_step();
    test_eigenvalue_computation();
    test_spectral_radius();

    /** L4: Fundamental Theorems */
    printf("\n--- L4: Fundamental Theorems ---\n");
    test_datarate_theorem_ct();
    test_rate_distortion();

    /** L5: Algorithms */
    printf("\n--- L5: Algorithms ---\n");
    test_zoom_quantizer();
    test_delta_modulation();
    test_logarithmic_quantizer();
    test_predictive_encoder();
    test_lqr_controller();
    test_huffman_coding();

    /** L6: Canonical Problems */
    printf("\n--- L6: Canonical Problems ---\n");
    test_event_triggered();
    test_tdma_scheduling();
    test_bandwidth_allocation();
    test_tod_protocol();
    test_adaptive_scheduling();

    /** L7/L8: Advanced */
    printf("\n--- L7/L8: Advanced Topics ---\n");
    test_lebesgue_sampler();
    test_event_trig_feedback();
    test_self_triggered();
    test_runlength_encoding();
    test_bit_allocation();

    printf("\n=== Results: %d / %d tests passed ===\n",
           tests_passed, tests_run);

    if (tests_passed == tests_run) {
        printf("ALL TESTS PASSED ✅\n");
        return 0;
    } else {
        printf("SOME TESTS FAILED ⚠️\n");
        return 1;
    }
}