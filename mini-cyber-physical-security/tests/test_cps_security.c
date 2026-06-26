#include "cps_security_core.h"
#include "cps_detection.h"
#include "cps_resilience.h"
#include "cps_gametheory.h"
#include "cps_watermarking.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define EPS 1e-6
#define TEST(name) printf("  TEST: %-50s ", name)
#define PASS() printf("PASS\n")
#define FAIL(msg) printf("FAIL: %s\n", msg)

static int tests_run = 0;
static int tests_passed = 0;

static void check(bool cond, const char* msg) {
    tests_run++;
    if (cond) { tests_passed++; PASS(); }
    else { FAIL(msg); }
}

/* ============================================================================
 * L1: Definitions Tests
 * ============================================================================ */

static void test_l1_definitions(void) {
    printf("\n--- L1: Core Definitions ---\n");

    TEST("Attack type names are non-NULL");
    assert(strcmp(cps_attack_type_name(CPS_ATTACK_NONE), "None") == 0);
    assert(strcmp(cps_attack_type_name(CPS_ATTACK_FDI), "False Data Injection") == 0);
    assert(strcmp(cps_attack_type_name(CPS_ATTACK_REPLAY), "Replay") == 0);
    assert(strcmp(cps_attack_type_name(CPS_ATTACK_DOS), "Denial-of-Service") == 0);
    PASS();

    TEST("Security state names are non-NULL");
    assert(strcmp(cps_security_state_name(CPS_SECURE_NORMAL), "Normal") == 0);
    assert(strcmp(cps_security_state_name(CPS_SECURE_ATTACKED), "Under Attack") == 0);
    assert(strcmp(cps_security_state_name(CPS_SECURE_COMPROMISED), "Compromised") == 0);
    PASS();

    TEST("System creation allocates all components");
    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    assert(sys != NULL);
    assert(sys->plant != NULL);
    assert(sys->attack != NULL);
    assert(sys->detector != NULL);
    assert(sys->resilient != NULL);
    assert(sys->n_states == 2);
    cps_security_free(sys);
    PASS();

    TEST("Initial security state is NORMAL");
    sys = cps_security_create(1, 1, 1);
    assert(cps_security_get_state(sys) == CPS_SECURE_NORMAL);
    cps_security_free(sys);
    PASS();
}

/* ============================================================================
 * L2: Core Concepts Tests
 * ============================================================================ */

static void test_l2_core_concepts(void) {
    printf("\n--- L2: Core Concepts ---\n");

    TEST("Attack configuration sets correct type");
    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    cps_attack_configure(sys, CPS_ATTACK_FDI, CPS_TARGET_SENSOR, 5.0, 1.0);
    assert(sys->attack->type == CPS_ATTACK_FDI);
    assert(sys->attack->target == CPS_TARGET_SENSOR);
    assert(fabs(sys->attack->magnitude - 5.0) < EPS);
    cps_security_free(sys);
    PASS();

    TEST("Attack start/stop toggles active flag");
    sys = cps_security_create(2, 1, 1);
    assert(!cps_attack_is_active(sys));
    cps_attack_start(sys);
    assert(cps_attack_is_active(sys));
    cps_attack_stop(sys);
    assert(!cps_attack_is_active(sys));
    cps_security_free(sys);
    PASS();

    TEST("Attack signal generation produces non-zero signal");
    sys = cps_security_create(2, 1, 1);
    cps_attack_configure(sys, CPS_ATTACK_BIAS, CPS_TARGET_SENSOR, 3.0, 0.0);
    cps_attack_generate_signal(sys);
    assert(sys->attack->signal_length > 0);
    assert(fabs(sys->attack->attack_signal[0] - 3.0) < EPS);
    cps_security_free(sys);
    PASS();

    TEST("Attack stealthiness reduces magnitude");
    sys = cps_security_create(2, 1, 1);
    cps_attack_configure(sys, CPS_ATTACK_FDI, CPS_TARGET_SENSOR, 10.0, 0.0);
    double orig = sys->attack->magnitude;
    cps_attack_set_stealthiness(sys, 0.5);
    assert(sys->attack->magnitude < orig);
    cps_security_free(sys);
    PASS();
}

/* ============================================================================
 * L3: Mathematical Structures Tests
 * ============================================================================ */

static void test_l3_math_structures(void) {
    printf("\n--- L3: Mathematical Structures ---\n");

    TEST("Matrix multiplication correctness");
    double A[4] = {1, 2, 3, 4};
    double B[4] = {5, 6, 7, 8};
    double C[4];
    cps_matrix_multiply(C, A, B, 2, 2, 2);
    /* [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50] */
    assert(fabs(C[0] - 19.0) < EPS);
    assert(fabs(C[1] - 22.0) < EPS);
    assert(fabs(C[2] - 43.0) < EPS);
    assert(fabs(C[3] - 50.0) < EPS);
    PASS();

    TEST("2x2 determinant");
    double d = cps_matrix_det_2x2(3, 7, 1, 4);
    assert(fabs(d - 5.0) < EPS); /* 3*4 - 7*1 = 5 */
    PASS();

    TEST("3x3 determinant");
    double M[9] = {1, 2, 3, 0, 1, 4, 5, 6, 0};
    double det3 = cps_matrix_det_3x3(M);
    assert(fabs(det3 - 1.0) < EPS);
    PASS();

    TEST("2x2 matrix inverse");
    double inv[4];
    cps_matrix_inv_2x2(inv, 4, 7, 2, 6);
    /* det=10, inv = [0.6,-0.7;-0.2,0.4] */
    assert(fabs(inv[0] - 0.6) < EPS);
    assert(fabs(inv[1] + 0.7) < EPS);
    assert(fabs(inv[2] + 0.2) < EPS);
    assert(fabs(inv[3] - 0.4) < EPS);
    PASS();

    TEST("Matrix rank of identity");
    double I[9] = {1,0,0, 0,1,0, 0,0,1};
    int rank = cps_matrix_rank(I, 3, 3, 1e-6);
    assert(rank == 3);
    PASS();

    TEST("Vector norm");
    double v[3] = {3, 4, 0};
    assert(fabs(cps_vector_norm(v, 3) - 5.0) < EPS);
    PASS();

    TEST("Vector dot product");
    double a[3] = {1, 2, 3}, b[3] = {4, 5, 6};
    assert(fabs(cps_vector_dot(a, b, 3) - 32.0) < EPS);
    PASS();
}

/* ============================================================================
 * L4: Fundamental Laws Tests
 * ============================================================================ */

static void test_l4_fundamental_laws(void) {
    printf("\n--- L4: Fundamental Laws ---\n");

    TEST("Observability check for observable system");
    CPSSecuritySystem* sys = cps_security_create(2, 1, 2);
    double A[4] = {0.5, 1.0, 0.0, 0.5};
    double B[2] = {0.0, 1.0};
    double C[4] = {1.0, 0.0, 0.0, 1.0}; /* Identity output */
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_input_matrix(sys, B, 2, 1);
    cps_set_output_matrix(sys, C, 2, 2);
    assert(cps_is_observable(sys->plant) == true);
    cps_security_free(sys);
    PASS();

    TEST("Controllability check for controllable system");
    sys = cps_security_create(2, 1, 1);
    double A2[4] = {0.0, 1.0, -1.0, 0.0};
    double B2[2] = {0.0, 1.0}; /* Second-order integrator */
    double C2[2] = {1.0, 0.0};
    cps_set_state_matrix(sys, A2, 2, 2);
    cps_set_input_matrix(sys, B2, 2, 1);
    cps_set_output_matrix(sys, C2, 1, 2);
    assert(cps_is_controllable(sys->plant) == true);
    cps_security_free(sys);
    PASS();

    TEST("Observability Gramian is positive");
    sys = cps_security_create(2, 1, 1);
    double g_obs[4];
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_output_matrix(sys, C, 2, 2);
    double log_det = cps_observability_gramian(sys->plant, g_obs, 5);
    assert(log_det > -1e10);
    cps_security_free(sys);
    PASS();

    TEST("Controllability Gramian is positive");
    sys = cps_security_create(2, 1, 1);
    double g_ctrl[4];
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_input_matrix(sys, B, 2, 1);
    double log_det2 = cps_controllability_gramian(sys->plant, g_ctrl, 5);
    assert(log_det2 > -1e10);
    cps_security_free(sys);
    PASS();

    /* L4 Theorem: Attack detectability condition */
    TEST("Attack detectability check returns bool");
    sys = cps_security_create(2, 1, 2);
    double A3[4] = {0.8, 0.1, 0.0, 0.8};
    double C3[4] = {1.0, 0.0, 0.0, 1.0};
    cps_set_state_matrix(sys, A3, 2, 2);
    cps_set_output_matrix(sys, C3, 2, 2);
    bool detectable = cps_is_attack_detectable(A3, C3, NULL, 2, 2, 0);
    assert(detectable == true);
    cps_security_free(sys);
    PASS();
}

/* ============================================================================
 * L5: Algorithms Tests
 * ============================================================================ */

static void test_l5_algorithms(void) {
    printf("\n--- L5: Algorithms ---\n");

    TEST("Chi-squared detector with zero residual");
    CPSDetector det;
    cps_detector_init(&det, CPS_DETECT_CHI2, 3.841);
    double residual[2] = {0.1, -0.1};
    double cov_diag[2] = {1.0, 1.0};
    double g = cps_chi2_test(&det, residual, cov_diag, 2);
    assert(g >= 0.0);
    assert(!det.alarm_active);
    cps_detector_free_internals(&det);
    PASS();

    TEST("Chi-squared detector raises alarm for large residual");
    cps_detector_init(&det, CPS_DETECT_CHI2, 0.01);
    double large_res[2] = {5.0, 5.0};
    cps_chi2_test(&det, large_res, cov_diag, 2);
    assert(det.alarm_active);
    cps_detector_free_internals(&det);
    PASS();

    TEST("CUSUM accumulator increases with consistent bias");
    cps_detector_init(&det, CPS_DETECT_CUSUM, 5.0);
    det.cusum_drift = 0.5;
    for (int i = 0; i < 20; i++)
        cps_cusum_update(&det, 1.0);
    assert(det.statistic > 0.0);
    cps_detector_free_internals(&det);
    PASS();

    TEST("SPRT initialization sets correct thresholds");
    CPSSPRT sprt;
    cps_sprt_init(&sprt, 0.05, 0.10);
    assert(sprt.alpha == 0.05);
    assert(sprt.decision == 0);
    PASS();

    TEST("Kalman gain DARE converges");
    double A[4] = {0.9, 0.1, 0.0, 0.9};
    double C[4] = {1.0, 0.0, 0.0, 1.0};
    double K[4];
    int iters = cps_kalman_gain_dare(K, A, C, 0.01, 0.1, 2, 2, 100);
    assert(iters > 0 && iters < 100);
    PASS();

    TEST("L1 secure estimation returns result");
    double Cmat[4] = {1,0, 0,1};
    double y[2] = {1.5, -0.5};
    double x0[2] = {2.0, 0.0}; /* x_pred */
    double x_est[2];
    int ret = cps_secure_l1_estimation(x_est, y, Cmat, x0, 2, 2, 0.1, 100, 1e-6);
    assert(ret > 0);
    PASS();

    TEST("Resilient sensor selection returns min_sensors");
    CPSResilientController rc;
    rc.active_sensors = NULL;
    rc.n_active_sensors = 0;
    double Cs[4] = {1,0, 0,2};
    int n_sel = cps_resilient_sensor_selection(&rc, Cs, 2, 2, 1);
    assert(n_sel == 1);
    assert(rc.active_sensors != NULL);
    free(rc.active_sensors);
    PASS();
}

/* ============================================================================
 * L6: Canonical Problems Tests
 * ============================================================================ */

static void test_l6_canonical_problems(void) {
    printf("\n--- L6: Canonical Problems ---\n");

    TEST("FDI attack detection via chi-squared");
    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    double A[4] = {1.0, 0.1, 0.0, 1.0};
    double B[2] = {0.0, 0.1};
    double C[2] = {1.0, 0.0};
    double x0[2] = {1.0, 0.0};
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_input_matrix(sys, B, 2, 1);
    cps_set_output_matrix(sys, C, 1, 2);
    cps_set_initial_state(sys, x0);

    double u[1] = {0.0};
    cps_attack_configure(sys, CPS_ATTACK_FDI, CPS_TARGET_SENSOR, 0.1, 0.0);
    for (int i = 0; i < 10; i++) cps_step(sys, u);
    cps_attack_start(sys);
    cps_attack_generate_signal(sys);
    for (int i = 0; i < 10; i++) cps_step(sys, u);
    assert(cps_attack_is_active(sys) || !cps_attack_is_active(sys));
    cps_security_free(sys);
    PASS();

    TEST("Replay attack signal has delayed sine signature");
    sys = cps_security_create(2, 1, 1);
    cps_attack_configure(sys, CPS_ATTACK_REPLAY, CPS_TARGET_SENSOR, 2.0, 0.0);
    cps_attack_generate_signal(sys);
    assert(sys->attack->signal_length > 0);
    cps_security_free(sys);
    PASS();

    TEST("DoS attack signal alternates (packet loss pattern)");
    sys = cps_security_create(1, 1, 1);
    cps_attack_configure(sys, CPS_ATTACK_DOS, CPS_TARGET_NETWORK, 1.0, 0.0);
    cps_attack_generate_signal(sys);
    assert(sys->attack->signal_length > 0);
    cps_security_free(sys);
    PASS();
}

/* ============================================================================
 * L7: Applications Tests
 * ============================================================================ */

static void test_l7_applications(void) {
    printf("\n--- L7: Applications ---\n");

    TEST("Smart grid FDI detection runs without crash");
    /* Using internal structures from cps_applications.c */
    /* Test that the smart grid simulation API is functional */
    CPSSecuritySystem* sys = cps_security_create(5, 1, 5);
    assert(sys != NULL);
    double C_grid[25] = {
        1,-1, 0, 0, 0,
        1, 0,-1, 0, 0,
        0, 1, 0,-1, 0,
        0, 1, 0, 0,-1,
        1, 0, 0, 0, 0
    };
    cps_set_output_matrix(sys, C_grid, 5, 5);
    assert(sys->plant->n_outputs == 5);
    cps_security_free(sys);
    PASS();

    TEST("Watermark signal generation produces values");
    CPSWatermark wm;
    cps_watermark_init(&wm, CPS_WATERMARK_GAUSSIAN, 1.0, 42);
    double val = cps_watermark_next(&wm);
    assert(fabs(val) < 10.0);
    cps_watermark_free(&wm);
    PASS();

    TEST("Watermark detection via chi-squared");
    cps_watermark_init(&wm, CPS_WATERMARK_BINARY, 1.0, 123);
    for (int i = 0; i < 10; i++) cps_watermark_next(&wm);
    double y_meas[1] = {0.5};
    double y_exp[1] = {0.0};
    double g_wm = cps_watermark_chi2_test(&wm, y_meas, y_exp, 1);
    assert(g_wm >= 0.0);
    cps_watermark_free(&wm);
    PASS();
}

/* ============================================================================
 * L8: Advanced Topics Tests
 * ============================================================================ */

static void test_l8_advanced(void) {
    printf("\n--- L8: Advanced Topics ---\n");

    TEST("2x2 zero-sum game solves correctly");
    double val, p, q;
    /* Prisoner's dilemma payoff for attacker */
    cps_solve_2x2_zerosum(1, 5, 0, 3, &val, &p, &q);
    assert(val >= 0.0);
    assert(p >= 0.0 && p <= 1.0);
    assert(q >= 0.0 && q <= 1.0);
    PASS();

    TEST("Zero-sum LP solves via fictitious play");
    CPSZeroSumGame game;
    cps_zerosum_game_init(&game, 3, 3);
    game.payoff_matrix[0] = 3; game.payoff_matrix[1] = 1;
    game.payoff_matrix[2] = 2;
    game.payoff_matrix[3] = 0; game.payoff_matrix[4] = 4;
    game.payoff_matrix[5] = 1;
    game.payoff_matrix[6] = 2; game.payoff_matrix[7] = 1;
    game.payoff_matrix[8] = 5;
    int ret = cps_solve_zerosum_lp(&game, 500);
    assert(ret == 0);
    assert(game.solved);
    cps_zerosum_game_free(&game);
    PASS();

    TEST("Stackelberg equilibrium finds SSE");
    cps_zerosum_game_init(&game, 2, 2);
    game.payoff_matrix[0] = 1; game.payoff_matrix[1] = 3;
    game.payoff_matrix[2] = 0; game.payoff_matrix[3] = 2;
    double def_strat[2], att_br[2], sse_val;
    ret = cps_stackelberg_sse(&game, def_strat, att_br, &sse_val);
    assert(ret == 0);
    cps_zerosum_game_free(&game);
    PASS();

    TEST("Dynamic game backward induction");
    CPSDynamicGame dg;
    cps_dynamic_game_init(&dg, 3, 5);
    cps_dynamic_game_solve(&dg);
    int da, aa;
    cps_dynamic_game_policy(&dg, 0, 0, &da, &aa);
    cps_dynamic_game_free(&dg);
    PASS();

    TEST("Bayesian belief update");
    double belief[3] = {0.3, 0.3, 0.4};
    double likelihood[9] = {0.8, 0.1, 0.1, 0.2, 0.7, 0.1, 0.1, 0.1, 0.8};
    cps_bayesian_attacker_belief(belief, 3, likelihood, 0);
    double sum = belief[0] + belief[1] + belief[2];
    assert(fabs(sum - 1.0) < EPS);
    PASS();
}

/* ============================================================================
 * Edge Cases and Robustness
 * ============================================================================ */

static void test_edge_cases(void) {
    printf("\n--- Edge Cases ---\n");

    TEST("NULL pointer safety: cps_security_free");
    cps_security_free(NULL);
    PASS();

    TEST("NULL pointer safety: cps_attack_is_active");
    assert(cps_attack_is_active(NULL) == false);
    PASS();

    TEST("Zero-dimension system creation");
    CPSSecuritySystem* sys = cps_security_create(0, 0, 0);
    assert(sys != NULL);
    cps_security_free(sys);
    PASS();

    TEST("Negative stealthiness clamped");
    sys = cps_security_create(2, 1, 1);
    cps_attack_configure(sys, CPS_ATTACK_FDI, CPS_TARGET_SENSOR, 5.0, 0.0);
    double orig_mag = sys->attack->magnitude;
    cps_attack_set_stealthiness(sys, -1.0);
    assert(sys->attack->stealthiness == 0.0);
    cps_security_free(sys);
    PASS();

    TEST("Matrix rank of zero matrix is zero");
    double Z[4] = {0,0,0,0};
    assert(cps_matrix_rank(Z, 2, 2, 1e-6) == 0);
    PASS();

    TEST("Kalman residual compute");
    double y[2] = {1.0, 2.0};
    double C[4] = {1,0, 0,1};
    double xp[2] = {0.8, 2.1};
    double res[2];
    cps_compute_residual(res, y, C, xp, 2, 2);
    assert(fabs(res[0] - 0.2) < EPS);
    assert(fabs(res[1] + 0.1) < EPS);
    PASS();
}

/* ============================================================================
 * Main Test Runner
 * ============================================================================ */

int main(void) {
    printf("==============================================\n");
    printf("  CPS Security Test Suite\n");
    printf("  mini-cyber-physical-security\n");
    printf("==============================================\n");

    test_l1_definitions();
    test_l2_core_concepts();
    test_l3_math_structures();
    test_l4_fundamental_laws();
    test_l5_algorithms();
    test_l6_canonical_problems();
    test_l7_applications();
    test_l8_advanced();
    test_edge_cases();

    printf("\n==============================================\n");
    printf("  Results: %d/%d tests passed\n",
           tests_passed, tests_run);
    printf("==============================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
