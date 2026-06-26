#include "ebc_core.h"
#include "ebc_trigger.h"
#include "ebc_stability.h"
#include "ebc_self.h"
#include "ebc_periodic.h"
#include "ebc_performance.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
 * test_ebc.c -- Comprehensive test suite for mini-event-based-control
 *
 * Tests cover:
 *   - System lifecycle (create, set_state, reset, free)
 *   - Event detection (all trigger types)
 *   - Lyapunov equation solver
 *   - ISS gain and critical sigma computation
 *   - Minimum inter-event time bound
 *   - Stability certificate generation
 *   - Matrix exponential computation
 *   - PETC state machine
 *   - Performance metric computation
 *   - Edge cases (NULL, zero-dim, degenerate matrices)
 */

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s: ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)

/* ================================================================
 * Test 1: System lifecycle (L1)
 * ================================================================ */
static void test_system_lifecycle(void) {
    TEST("system_create_valid");
    EBC_System* sys = ebc_system_create(3, 2, EBC_CONTINUOUS_ETC);
    assert(sys != NULL);
    assert(sys->n == 3);
    assert(sys->m == 2);
    assert(sys->paradigm == EBC_CONTINUOUS_ETC);
    PASS();

    TEST("system_create_invalid_n");
    EBC_System* s2 = ebc_system_create(0, 1, EBC_CONTINUOUS_ETC);
    assert(s2 == NULL);
    PASS();

    TEST("system_create_invalid_m");
    EBC_System* s3 = ebc_system_create(3, 0, EBC_CONTINUOUS_ETC);
    assert(s3 == NULL);
    PASS();

    TEST("system_set_state");
    double x0[] = {1.0, 2.0, 3.0};
    ebc_system_set_state(sys, x0);
    assert(fabs(sys->x[0] - 1.0) < 1e-10);
    assert(fabs(sys->x[1] - 2.0) < 1e-10);
    assert(fabs(sys->x[2] - 3.0) < 1e-10);
    PASS();

    TEST("system_reset");
    sys->x[0] = 5.0;
    ebc_system_reset(sys);
    assert(fabs(sys->x_last[0] - 5.0) < 1e-10);
    assert(sys->event_count == 0);
    PASS();

    TEST("system_free");
    ebc_system_free(sys);
    PASS();
}

/* ================================================================
 * Test 2: Event detection (L2)
 * ================================================================ */
static void test_event_detection(void) {
    TEST("check_event_below_threshold");
    EBC_System* sys = ebc_system_create(2, 1, EBC_CONTINUOUS_ETC);
    double x0[] = {1.0, 0.0};
    ebc_system_set_state(sys, x0);
    EBC_TriggerParams tp = ebc_trigger_mixed(0.5, 0.1);
    sys->x[0] = 1.01;
    bool evt = ebc_check_event(sys, &tp);
    assert(evt == false);
    PASS();

    TEST("check_event_above_threshold");
    sys->x[0] = 3.0;
    evt = ebc_check_event(sys, &tp);
    assert(evt == true);
    PASS();

    TEST("mark_event");
    int cnt = ebc_mark_event(sys, 1.0);
    assert(cnt == 1);
    assert(sys->event_count == 1);
    assert(sys->t_last == 1.0);
    PASS();

    TEST("error_norm_computation");
    sys->x[0] = 2.0;
    sys->x[1] = 2.0;
    sys->x_last[0] = 1.0;
    sys->x_last[1] = 1.0;
    double en = ebc_compute_error_norm(sys);
    assert(fabs(en - sqrt(2.0)) < 1e-8);
    PASS();

    ebc_system_free(sys);
}

/* ================================================================
 * Test 3: Trigger condition types (L5)
 * ================================================================ */
static void test_trigger_types(void) {
    TEST("trigger_absolute");
    EBC_System* sys = ebc_system_create(2, 1, EBC_CONTINUOUS_ETC);
    double x0[] = {1.0, 0.0};
    ebc_system_set_state(sys, x0);
    EBC_TriggerParams tp = ebc_trigger_default();
    tp.type = EBC_ABSOLUTE_ERROR;
    tp.epsilon = 0.5;
    sys->x[0] = 1.6;
    assert(ebc_trigger_absolute(sys, &tp) == true);
    sys->x[0] = 1.2;
    assert(ebc_trigger_absolute(sys, &tp) == false);
    PASS();

    TEST("trigger_relative");
    tp.type = EBC_RELATIVE_ERROR;
    tp.sigma = 0.3;
    sys->x[0] = 1.0; sys->x_last[0] = 2.0;
    sys->x[1] = 0.0; sys->x_last[1] = 0.0;
    assert(ebc_trigger_relative(sys, &tp) == true);
    PASS();

    TEST("trigger_mixed_threshold");
    tp.type = EBC_MIXED_THRESHOLD;
    tp.sigma = 0.1; tp.epsilon = 0.01;
    sys->x[0] = 1.0; sys->x_last[0] = 1.0;
    assert(ebc_trigger_mixed_threshold(sys, &tp) == false);
    sys->x_last[0] = 2.0;
    assert(ebc_trigger_mixed_threshold(sys, &tp) == true);
    PASS();

    TEST("compute_threshold");
    tp.type = EBC_MIXED_THRESHOLD;
    tp.sigma = 0.5; tp.epsilon = 0.1;
    double th = ebc_compute_threshold(sys, &tp);
    assert(fabs(th - 0.6) < 1e-10);
    PASS();

    ebc_system_free(sys);
}

/* ================================================================
 * Test 4: Lyapunov equation solver (L4)
 * ================================================================ */
static void test_lyapunov_solver(void) {
    TEST("lyapunov_solve_2x2");
    double A[] = {-2.0, 1.0, 0.0, -3.0};
    double B[] = {1.0, 0.0};
    double K[] = {-1.0, 0.0};
    double Q[] = {1.0, 0.0, 0.0, 1.0};
    double P[4] = {0};
    int ret = ebc_lyapunov_solve(A, B, K, 2, 1, Q, P);
    assert(ret == 0);
    assert(P[0] > 0.0);
    assert(P[3] > 0.0);
    PASS();

    TEST("lyapunov_solve_null_input");
    ret = ebc_lyapunov_solve(NULL, B, K, 2, 1, Q, P);
    assert(ret == -1);
    PASS();
}

/* ================================================================
 * Test 5: ISS gain and critical sigma (L4)
 * ================================================================ */
static void test_iss_and_sigma(void) {
    TEST("iss_gain_computation");
    double A[] = {-2.0, 0.0, 0.0, -2.0};
    double B[] = {1.0, 0.0};
    double K[] = {-1.0, 0.0};
    double gamma = ebc_iss_gain_linear(A, B, K, 2, 1);
    assert(gamma > 0.0);
    PASS();

    TEST("critical_sigma_positive");
    double sig_crit = ebc_critical_sigma(A, B, K, 2, 1);
    assert(sig_crit > 0.0);
    assert(sig_crit < 1.0);
    PASS();

    TEST("minimum_iet_positive");
    double tau_min = ebc_minimum_iet_linear(A, B, K, 2, 1, 0.1, 0.01);
    assert(tau_min > 0.0);
    PASS();
}

/* ================================================================
 * Test 6: Stability certificate (L4)
 * ================================================================ */
static void test_stability_certificate(void) {
    TEST("stability_certify_linear");
    double A[] = {-2.0, 0.0, 0.0, -3.0};
    double B[] = {1.0, 0.0};
    double K[] = {-1.0, 0.0};
    EBC_StabilityCert cert;
    int ret = ebc_stability_certify_linear(A, B, K, 2, 1, 0.1, 0.01, &cert);
    assert(ret == 0);
    assert(cert.alpha1 > 0.0);
    assert(cert.alpha2 > cert.alpha1);
    assert(cert.sigma_critical > 0.0);
    ebc_stability_cert_free(&cert);
    PASS();

    TEST("stability_cert_free_null");
    ebc_stability_cert_free(NULL);
    PASS();
}

/* ================================================================
 * Test 7: Matrix exponential (L5)
 * ================================================================ */
static void test_matrix_exponential(void) {
    TEST("matrix_exp_identity");
    double A[] = {1.0, 0.0, 0.0, 1.0};
    double E[4];
    ebc_matrix_exponential(A, 2, 0.0, E);
    assert(fabs(E[0] - 1.0) < 0.01);
    assert(fabs(E[1] - 0.0) < 0.01);
    assert(fabs(E[3] - 1.0) < 0.01);
    PASS();

    TEST("matrix_exp_null");
    ebc_matrix_exponential(NULL, 2, 1.0, E);
    PASS();
}

/* ================================================================
 * Test 8: Performance metrics (L6)
 * ================================================================ */
static void test_performance_metrics(void) {
    TEST("ise_computation");
    double traj[] = {1.0, 0.0, 0.5, 0.0, 0.25, 0.0, 0.0, 0.0};
    double ise = ebc_ise(traj, 4, 2, 0.1);
    assert(ise > 0.0);
    assert(ise < 1.0);
    PASS();

    TEST("iae_computation");
    double iae = ebc_iae(traj, 4, 2, 0.1);
    assert(iae > 0.0);
    PASS();

    TEST("itae_computation");
    double itae = ebc_itae(traj, 4, 2, 0.1);
    assert(itae > 0.0);
    PASS();

    TEST("average_iet");
    double events[] = {0.0, 0.5, 1.0, 1.5};
    double avg = ebc_average_iet(events, 4);
    assert(fabs(avg - 0.5) < 1e-10);
    PASS();

    TEST("settling_time");
    double st = ebc_settling_time(traj, 4, 2, 0.1, 0.02);
    assert(st >= 0.0);
    PASS();

    TEST("overshoot_non_negative");
    double os = ebc_overshoot(traj, 4, 2, 0.1);
    assert(os >= -0.01);
    PASS();
}

/* ================================================================
 * Test 9: PETC (L5)
 * ================================================================ */
static void test_petc(void) {
    TEST("petc_config_default");
    EBC_PETC_Config cfg = ebc_petc_config_default();
    assert(cfg.h == 0.01);
    assert(cfg.sigma == 0.1);
    assert(cfg.max_skip == 100);
    PASS();

    TEST("petc_create");
    EBC_PETC_System petsc = ebc_petc_create(cfg, 1.0, 0.1);
    assert(petsc.state == PETC_WAIT);
    assert(petsc.sample_count == 0);
    PASS();
}

/* ================================================================
 * Test 10: Edge cases (robustness)
 * ================================================================ */
static void test_edge_cases(void) {
    TEST("vector_norm_zero");
    double v[] = {0.0, 0.0, 0.0};
    double n = ebc_vector_norm(v, 3);
    assert(fabs(n) < 1e-10);
    PASS();

    TEST("vector_norm_null");
    n = ebc_vector_norm(NULL, 3);
    assert(fabs(n) < 1e-10);
    PASS();

    TEST("matrix_norm_zero");
    double M[] = {0.0, 0.0, 0.0, 0.0};
    double mn = ebc_matrix_norm(M, 2);
    assert(fabs(mn) < 1e-10);
    PASS();

    TEST("trigger_free_null");
    ebc_trigger_free(NULL);
    PASS();

    TEST("zeno_detect_normal");
    double evts[] = {0.0, 0.5, 1.0, 1.5, 2.0};
    bool zeno = ebc_detect_zeno(evts, 5, 0.001);
    assert(zeno == false);
    PASS();

    TEST("classify_iet_positive");
    EBC_IET_Class cls = ebc_classify_iet(evts, 5);
    assert(cls == EBC_IET_POSITIVE);
    PASS();

    TEST("integral_square");
    double xs[] = {3.0, 4.0};
    double isq = ebc_integral_square(xs, 2);
    assert(fabs(isq - 25.0) < 1e-10);
    PASS();
}

/* ================================================================
 * Test 11: Trigger parameters (L2)
 * ================================================================ */
static void test_trigger_params(void) {
    TEST("trigger_default_values");
    EBC_TriggerParams tp = ebc_trigger_default();
    assert(tp.type == EBC_MIXED_THRESHOLD);
    assert(fabs(tp.sigma - 0.1) < 1e-10);
    assert(fabs(tp.epsilon - 0.01) < 1e-10);
    PASS();

    TEST("trigger_relative_constructor");
    tp = ebc_trigger_make_relative(0.5);
    assert(tp.type == EBC_RELATIVE_ERROR);
    assert(fabs(tp.sigma - 0.5) < 1e-10);
    PASS();

    TEST("trigger_mixed_constructor");
    tp = ebc_trigger_mixed(0.3, 0.05);
    assert(tp.type == EBC_MIXED_THRESHOLD);
    assert(fabs(tp.sigma - 0.3) < 1e-10);
    assert(fabs(tp.epsilon - 0.05) < 1e-10);
    PASS();

    TEST("trigger_clamp_sigma");
    tp = ebc_trigger_make_relative(5.0);
    assert(fabs(tp.sigma - 0.1) < 1e-10);
    PASS();
}

/* ================================================================
 * Test 12: Self-triggered (L5)
 * ================================================================ */
static void test_self_triggered(void) {
    TEST("self_next_time_linear");
    double A[] = {-1.0, 0.0, 0.0, -2.0};
    double B[] = {1.0, 0.0};
    double K[] = {0.0, 0.0};
    double xk[] = {1.0, 0.0};
    double tau = ebc_self_next_time_linear(A, B, K, xk, 2, 1, 0.1, 0.01, 1.0, 1e-4);
    assert(tau > 0.0);
    assert(tau <= 1.0);
    PASS();

    TEST("matrix_exp_integral");
    double G[4];
    ebc_matrix_exp_integral(A, 2, 0.1, G);
    assert(G[0] > 0.0);
    PASS();
}

/* ================================================================
 * Test 13: Dynamic trigger (L8)
 * ================================================================ */
static void test_dynamic_trigger(void) {
    TEST("dynamic_trigger_create");
    EBC_DynamicTrigger dt = ebc_dynamic_trigger_create(1.0, 0.5, 0.3, 0.1);
    assert(fabs(dt.eta - 1.0) < 1e-10);
    assert(fabs(dt.beta - 0.5) < 1e-10);
    PASS();

    TEST("dynamic_trigger_update");
    double x[] = {1.0, 0.0};
    double e[] = {0.2, 0.0};
    ebc_dynamic_trigger_update(&dt, x, e, 2, 0.01);
    assert(dt.eta > 0.0);
    PASS();

    TEST("dynamic_trigger_evaluate");
    bool fire = ebc_dynamic_trigger_evaluate(&dt, x, e, 2);
    assert(fire == false);
    PASS();
}

/* ================================================================
 * Test 14: Simulation (L5)
 * ================================================================ */
static void test_simulation(void) {
    TEST("system_step_euler");
    EBC_System* sys = ebc_system_create(2, 1, EBC_CONTINUOUS_ETC);
    double x0[] = {1.0, 0.0};
    ebc_system_set_state(sys, x0);
    EBC_TriggerParams tp = ebc_trigger_mixed(10.0, 100.0);
    EBC_Controller ctrl;
    memset(&ctrl, 0, sizeof(ctrl));
    int stat = ebc_step_euler(sys, &ctrl, 0.01, &tp);
    assert(stat >= -1);
    PASS();

    TEST("system_step_rk4");
    stat = ebc_step_rk4(sys, &ctrl, 0.01, &tp);
    assert(stat >= -1);
    PASS();

    TEST("system_simulate");
    double *traj = NULL, *events = NULL;
    int traj_len = 0, evt_len = 0;
    ebc_system_set_state(sys, x0);
    tp = ebc_trigger_mixed(0.5, 0.1);
    int ret = ebc_simulate(sys, &ctrl, 0.1, 0.01, &tp,
                            &traj, &traj_len, &events, &evt_len);
    assert(ret == 0);
    assert(traj_len > 0);
    assert(evt_len > 0);
    free(traj);
    free(events);
    PASS();

    ebc_system_free(sys);
}

/* ================================================================
 * Main
 * ================================================================ */
int main(void) {
    printf("=== mini-event-based-control Test Suite ===\n\n");

    test_system_lifecycle();
    test_event_detection();
    test_trigger_types();
    test_lyapunov_solver();
    test_iss_and_sigma();
    test_stability_certificate();
    test_matrix_exponential();
    test_performance_metrics();
    test_petc();
    test_edge_cases();
    test_trigger_params();
    test_self_triggered();
    test_dynamic_trigger();
    test_simulation();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
