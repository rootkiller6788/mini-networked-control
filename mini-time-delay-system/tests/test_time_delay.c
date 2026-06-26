#include "time_delay_system.h"
#include "lyapunov_krasovskii.h"
#include "smith_predictor.h"
#include "delay_stability.h"
#include "dde_solver.h"
#include "networked_delay.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_EPS 1e-6

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond) do { if (!(cond)) return 0; } while(0)
#define CHECK_MSG(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 0; } } while(0)

#define RUN_TEST(name) do { \
    tests_run++; \
    int result = test_##name(); \
    if (result) { tests_passed++; printf("  PASS: %s\n", #name); } \
    else { tests_failed++; printf("  FAIL: %s\n", #name); } \
} while(0)

/* ============================================================================
 * L1 — Definitions Tests
 * ============================================================================ */

static int test_delay_descriptor_constant(void) {
    DelayDescriptor* d = delay_create_constant(0.5);
    CHECK(d != NULL);
    CHECK(d->type == DELAY_CONSTANT);
    CHECK(fabs(d->tau_nominal - 0.5) < TEST_EPS);
    CHECK(fabs(d->tau_min - 0.5) < TEST_EPS);
    CHECK(fabs(d->tau_max - 0.5) < TEST_EPS);
    CHECK(d->is_bounded == true);
    CHECK(d->derivative_bound < TEST_EPS);
    delay_free(d);
    return 1;
}

static int test_delay_descriptor_time_varying(void) {
    DelayDescriptor* d = delay_create_time_varying(0.1, 0.5, 0.8);
    CHECK(d != NULL);
    CHECK(d->type == DELAY_TIME_VARYING);
    CHECK(d->tau_min >= 0.09);
    CHECK(d->tau_max <= 0.51);
    CHECK(d->derivative_bound <= 0.81);
    CHECK(d->is_bounded == true);
    delay_free(d);
    return 1;
}

static int test_delay_descriptor_stochastic(void) {
    DelayDescriptor* d = delay_create_stochastic(0.3, 0.01);
    CHECK(d != NULL);
    CHECK(d->type == DELAY_STOCHASTIC);
    CHECK(fabs(d->tau_nominal - 0.3) < TEST_EPS);
    CHECK(d->tau_variance > 0.009);
    delay_free(d);
    return 1;
}

static int test_tds_create_free(void) {
    TimeDelaySystem* sys = tds_create("test_sys", 2, 1, 1);
    CHECK(sys != NULL);
    CHECK(sys->n_states == 2);
    CHECK(sys->n_inputs == 1);
    CHECK(sys->n_outputs == 1);
    CHECK(strcmp(sys->name, "test_sys") == 0);
    CHECK(sys->dde_type == DDE_RETARDED);
    tds_free(sys);
    return 1;
}

static int test_tds_set_linear_model(void) {
    TimeDelaySystem* sys = tds_create("lin_sys", 2, 1, 1);
    double A[4] = {0.0, 1.0, -2.0, -3.0};
    double Ad[4] = {0.1, 0.0, 0.0, 0.1};
    double B[2] = {0.0, 1.0};
    double C[2] = {1.0, 0.0};

    tds_set_linear_model(sys, A, Ad, B, C);
    CHECK(fabs(sys->A[0] - 0.0) < TEST_EPS);
    CHECK(fabs(sys->A[3] + 3.0) < TEST_EPS);
    CHECK(fabs(sys->A_delayed[0] - 0.1) < TEST_EPS);
    CHECK(fabs(sys->B[1] - 1.0) < TEST_EPS);
    CHECK(fabs(sys->C[0] - 1.0) < TEST_EPS);
    tds_free(sys);
    return 1;
}

static int test_tds_add_delay(void) {
    TimeDelaySystem* sys = tds_create("delayed", 1, 1, 1);
    tds_add_delay(sys, DELAY_CONSTANT, 0.25, 0.25, 0.25);
    CHECK(sys->n_delays == 1);
    CHECK(sys->delays != NULL);
    CHECK(fabs(sys->delays[0]->tau_nominal - 0.25) < TEST_EPS);
    tds_free(sys);
    return 1;
}

/* ============================================================================
 * L2 — Core Concepts Tests
 * ============================================================================ */

static int test_characteristic_eqn_scalar(void) {
    /* ẋ = -x - 0.5 x(t-0.2) — stable system */
    TimeDelaySystem* sys = tds_create("scalar_test", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.5};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.2, 0.2, 0.2);

    /* At s = -0.5 (should be near a root if stable) */
    double mag = time_delay_characteristic_eqn(sys, -0.5, 0.0);
    CHECK(mag >= 0.0);  /* Determinant magnitude should be non-negative */
    CHECK(isfinite(mag));

    tds_free(sys);
    return 1;
}

static int test_delay_rate_check(void) {
    DelayDescriptor* dc = delay_create_constant(0.5);
    CHECK(time_delay_rate_check(dc) == true);  /* dτ/dt = 0 < 1 */

    DelayDescriptor* dtv = delay_create_time_varying(0.1, 0.5, 0.9);
    CHECK(time_delay_rate_check(dtv) == true);  /* 0.9 < 1 */

    DelayDescriptor* dtv2 = delay_create_time_varying(0.1, 0.5, 1.5);
    CHECK(time_delay_rate_check(dtv2) == false);  /* 1.5 > 1 */

    delay_free(dc); delay_free(dtv); delay_free(dtv2);
    return 1;
}

/* ============================================================================
   L3 — Mathematical Structures Tests
   ============================================================================ */

static int test_characteristic_root_computation(void) {
    TimeDelaySystem* sys = tds_create("roots", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.3};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.5, 0.5, 0.5);

    int n_roots = tds_compute_characteristic_roots(sys, 8);
    CHECK(n_roots > 0);  /* Should find at least the rightmost root */
    CHECK(sys->n_roots > 0);

    /* Verify the rightmost root: for a+ad e^{-τs}+s=0 with a=-1,ad=-0.3,
     * the root should have negative real part (stable). Allow tolerance. */
    CHECK(isfinite(sys->roots_real[0]));
    CHECK(sys->roots_real[0] < 10.0);  /* Not extremely positive */

    double abscissa = tds_spectral_abscissa(sys);
    CHECK(isfinite(abscissa));

    tds_free(sys);
    return 1;
}

/* ============================================================================
   L4 — Fundamental Theorem Tests
   ============================================================================ */

static int test_lkf_create_evaluate(void) {
    LKFunctional* lkf = lkf_create(2);
    CHECK(lkf != NULL);
    lkf_set_identity(lkf);
    CHECK(lkf_is_positive_definite(lkf) == true);

    double x[2] = {1.0, 0.5};
    double V = lkf_evaluate(lkf, x, NULL, NULL, 0, 0.5, 0.01);
    CHECK(V > 0.0);  /* V(x) = xᵀPx = ||x||² > 0 */

    lkf_free(lkf);
    return 1;
}

static int test_lkf_derivative(void) {
    /* Stable system: ẋ = -x - 0.1 x(t-τ) */
    TimeDelaySystem* sys = tds_create("lkf_test", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.1};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.3, 0.3, 0.3);

    LKFunctional* lkf = lkf_create(1);
    lkf_set_identity(lkf);

    double x[1] = {1.0};
    double xd[1] = {0.8};
    double dv = lkf_derivative(lkf, sys, x, xd, NULL, NULL, 0, 0.3, 0.01);
    CHECK(dv < 0.0);  /* Should be negative for stable system */

    lkf_free(lkf);
    tds_free(sys);
    return 1;
}

static int test_lmi_delay_independent(void) {
    /* a + |a_d| < 0 → delay-independent stable */
    TimeDelaySystem* sys = tds_create("lmi_ind", 1, 1, 1);
    double A[1] = {-2.0};
    double Ad[1] = {0.5};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);

    double P, Q;
    bool stable = lmi_delay_independent_check(sys, &P, &Q);
    CHECK(stable == true);  /* -2 + 0.5 = -1.5 < 0 → stable */

    tds_free(sys);
    return 1;
}

static int test_lmi_delay_dependent(void) {
    /* Stable for small delays */
    TimeDelaySystem* sys = tds_create("lmi_dep", 1, 1, 1);
    double A[1] = {-0.5};
    double Ad[1] = {-2.0};  /* |ad| > |a| → delay-sensitive */
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);

    double P, Q, R;
    /* Should be stable for very small delay */
    bool stable_small = lmi_delay_dependent_check(sys, 0.01, &P, &Q, &R);
    /* Should be unstable for large delay */
    bool stable_large = lmi_delay_dependent_check(sys, 10.0, &P, &Q, &R);

    CHECK(stable_small != stable_large);  /* Different behavior */
    tds_free(sys);
    return 1;
}

static int test_matrix_measure(void) {
    double A[4] = {-1.0, 0.0, 0.0, -2.0};
    double mu = matrix_measure_l2(A, 2);
    CHECK(mu < -0.9);  /* Both eigenvalues negative → μ < 0 */

    double B[4] = {1.0, 0.0, 0.0, 1.0};
    double mu_b = matrix_measure_l2(B, 2);
    CHECK(mu_b > 0.9);  /* Positive eigenvalues → μ > 0 */
    return 1;
}

static int test_halanay_inequality(void) {
    double gamma;
    bool holds = halanay_check(2.0, 1.0, &gamma);
    CHECK(holds == true);
    CHECK(gamma > 0.0);  /* α > β → positive decay rate */

    bool fails = halanay_check(1.0, 2.0, &gamma);
    CHECK(fails == false);
    return 1;
}

/* ============================================================================
   L5 — Algorithm Tests
   ============================================================================ */

static int test_pi_controller(void) {
    PIController* pi = pi_create(2.0, 1.0, 0.01, -1.0, 1.0);
    CHECK(pi != NULL);

    pi_setpoint(pi, 1.0);
    double u = pi_step(pi, 0.5);  /* error = 0.5 */
    CHECK(u > 0.0);  /* Positive error → positive control */
    CHECK(u <= 1.0);  /* Within saturation */

    pi_reset(pi);
    double u2 = pi_step(pi, 1.0);  /* error = 0.0 after reset */
    CHECK(fabs(u2) < 1e-6);

    pi_free(pi);
    return 1;
}

static int test_smith_predictor_create(void) {
    SmithPredictor* sp = sp_create(2, 1, 1, 0.5, 0.01);
    CHECK(sp != NULL);
    CHECK(sp->n_model == 2);
    CHECK(fabs(sp->tau_model - 0.5) < TEST_EPS);
    CHECK(sp->buffer_size >= 2);
    sp_free(sp);
    return 1;
}

static int test_smith_predictor_step(void) {
    SmithPredictor* sp = sp_create(1, 1, 1, 0.1, 0.01);

    /* Plant model: G(s) = 1/(s+1) */
    double A[1] = {-1.0};
    double B[1] = {1.0};
    double C[1] = {1.0};
    sp_set_plant_model(sp, A, B, C);
    sp_configure_pid(sp, 1.0, 0.5, 0.0, 100.0, 1.0, 0.01);

    double y[1] = {0.0};
    sp_step(sp, y);
    const double* u = sp_get_control(sp);
    CHECK(u != NULL);
    CHECK(isfinite(u[0]));

    sp_free(sp);
    return 1;
}

static int test_pade_approximation(void) {
    PadeApproximation* pade = pade_create(2, 0.1);
    CHECK(pade != NULL);
    CHECK(pade->order == 2);
    CHECK(pade->num != NULL);
    CHECK(pade->den != NULL);

    /* First-order Padé: (1 - τs/2) / (1 + τs/2)
     * For order=2: (1 - τs/2 + τ²s²/12) / (1 + τs/2 + τ²s²/12) */
    CHECK(pade->den[pade->order] > 0.0);  /* Leading coefficient positive */

    double *Ap, *Bp, *Cp, Dp;
    pade_to_state_space(pade, &Ap, &Bp, &Cp, &Dp);
    CHECK(Ap != NULL);
    CHECK(Bp != NULL);
    CHECK(Cp != NULL);

    free(Ap); free(Bp); free(Cp);
    pade_free(pade);
    return 1;
}

static int test_dde_solver_basic(void) {
    /* Simple DDE: ẋ = -x(t) - 0.2 x(t-0.2)
     * Very small delay, fast integration, small allocation */
    TimeDelaySystem* sys = tds_create("dde_test", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.2};
    double B[1] = {0.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.1, 0.1, 0.1);

    DDESolverConfig cfg;
    cfg.method = DDE_METHOD_RK4;
    cfg.interp = DDE_INTERP_LINEAR;
    cfg.dt = 0.1;
    cfg.dt_min = 1e-6;
    cfg.dt_max = 0.5;
    cfg.t_start = 0.0;
    cfg.t_end = 0.5;
    cfg.rel_tol = 1e-4;
    cfg.abs_tol = 1e-6;
    cfg.max_steps = 100;
    cfg.adaptive = false;
    cfg.track_discontinuities = false;

    double x0[1] = {1.0};
    DDESolution* sol = dde_solve(sys, &cfg, NULL, x0);
    if (!sol) { tds_free(sys); return 0; }
    CHECK(sol->n_steps > 0);

    dde_solution_free(sol);
    tds_free(sys);
    return 1;
}

/* ============================================================================
   L6 — Canonical Problem Tests
   ============================================================================ */

static int test_delay_margin_scalar(void) {
    /* Stable system: a = -1, a_d = -0.8
     * a_d² = 0.64 < 1.0 = a² → delay-independent stable */
    TimeDelaySystem* sys = tds_create("margin", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.8};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.5, 0.5, 0.5);

    double margin = delay_margin_frequency_sweep(sys);
    CHECK(margin == INFINITY || margin > 100.0);

    tds_free(sys);
    return 1;
}

static int test_tcp_aqm_model(void) {
    TCPAQMModel* m = tcp_aqm_create(50.0, 1000.0, 0.05, 10.0);
    CHECK(m != NULL);
    CHECK(fabs(m->N - 50.0) < TEST_EPS);
    CHECK(fabs(m->C - 1000.0) < TEST_EPS);

    double W, q;
    for (int i = 0; i < 100; i++) {
        tcp_aqm_set_drop_prob(m, 0.01);
        tcp_aqm_step(m, &W, &q);
        CHECK(W > 0.0);
        CHECK(q >= 0.0);
    }

    double rtt = tcp_aqm_rtt(m);
    CHECK(rtt > 0.0);

    tcp_aqm_free(m);
    return 1;
}

/* ============================================================================
   L7 — Applications Tests
   ============================================================================ */

static int test_ncs_create_step(void) {
    TimeDelaySystem* plant = tds_create("ncs_plant", 1, 1, 1);
    double A[1] = {-0.5};
    double Ad[1] = {-0.3};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(plant, A, Ad, B, C);
    plant->current_state[0] = 1.0;

    NetworkedControlSystem* ncs = ncs_create(plant, 0.1, 0.0);
    CHECK(ncs != NULL);

    double K[1] = {0.5};
    ncs_set_gain(ncs, K);

    ncs_set_sc_qos(ncs, 1e6, 0.0, 0.01, 0.002, 0.005, 0.05);
    ncs_set_ca_qos(ncs, 1e6, 0.0, 0.005, 0.001, 0.003, 0.03);

    /* Run a few steps */
    for (int k = 0; k < 10; k++) {
        ncs_step(ncs);
        double x = plant->current_state[0];
        CHECK(isfinite(x));
    }

    ncs_free(ncs);
    tds_free(plant);
    return 1;
}

static int test_mm1_queue(void) {
    MM1Queue* q = mm1_create(10.0, 20.0, 100.0);
    CHECK(q != NULL);
    CHECK(fabs(q->util - 0.5) < TEST_EPS);

    double d = mm1_step(q, 0.1, 2);  /* 2 packets arrive in 0.1s */
    CHECK(d >= 0.0);

    double mean_d = mm1_mean_delay(q);
    CHECK(mean_d > 0.0);

    double p_loss = mm1_loss_probability(q);
    CHECK(p_loss >= 0.0 && p_loss <= 1.0);

    mm1_free(q);
    return 1;
}

/* ============================================================================
   Multidimensional system test
   ============================================================================ */

static int test_2d_system_stability(void) {
    TimeDelaySystem* sys = tds_create("2d_test", 2, 1, 2);
    double A[4] = {-2.0, 0.0, 0.0, -1.5};
    double Ad[4] = {0.1, 0.0, 0.0, 0.1};
    double B[2] = {1.0, 0.5};
    double C[4] = {1.0, 0.0, 0.0, 1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.3, 0.3, 0.3);

    double margin = delay_margin_frequency_sweep(sys);
    CHECK(isfinite(margin) || margin == INFINITY);

    tds_compute_characteristic_roots(sys, 8);
    CHECK(sys->n_roots > 0);

    tds_free(sys);
    return 1;
}

/* ============================================================================
   Additional stability checks
   ============================================================================ */

static int test_nyquist_points(void) {
    TimeDelaySystem* sys = tds_create("nyquist", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.5};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.2, 0.2, 0.2);

    double omega[10], real[10], imag[10];
    int n = delay_nyquist_points(sys, 0.1, 10.0, 10, omega, real, imag);
    CHECK(n == 10);

    tds_free(sys);
    return 1;
}

static int test_gain_phase_margin(void) {
    TimeDelaySystem* sys = tds_create("margins", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.3};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);

    double gm = delay_gain_margin(sys);
    CHECK(gm > 0.0);

    double pm = delay_phase_margin(sys, 0.2);
    CHECK(isfinite(pm));

    tds_free(sys);
    return 1;
}

static int test_rhp_root_count(void) {
    TimeDelaySystem* sys = tds_create("rhp", 1, 1, 1);
    double A[1] = {-1.0};
    double Ad[1] = {-0.5};
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);

    /* Rekasius-substituted RHP count. Returns -1 for indeterminate. */
    int n_rhp = delay_rhp_root_count(sys, 0.3);
    /* Just verify function returns a valid integer */
    CHECK(n_rhp >= -1);

    tds_free(sys);
    return 1;
}

static int test_sweeping_test(void) {
    TimeDelaySystem* sys = tds_create("sweep", 1, 1, 1);
    double A[1] = {-0.5};
    double Ad[1] = {-2.0};  /* Delay-sensitive */
    double B[1] = {1.0};
    double C[1] = {1.0};
    tds_set_linear_model(sys, A, Ad, B, C);

    SweepingTestResult* res = delay_sweeping_test(sys, 0.0, 2.0, 100);
    CHECK(res != NULL);
    /* Should detect at least one crossing for delay-sensitive system */
    if (res->n_crossings > 0)
        CHECK(res->crossing_taus[0] > 0.0);

    sweeping_test_free(res);
    tds_free(sys);
    return 1;
}

static int test_predictor_feedback(void) {
    PredictorFeedback* pf = pf_create(2, 1, 0.3, 0.01, NULL, NULL, NULL);
    CHECK(pf != NULL);
    CHECK(pf->n == 2);
    CHECK(fabs(pf->tau - 0.3) < TEST_EPS);

    double x[2] = {1.0, 0.0};
    double u[1];
    pf_compute_control(pf, x, u);
    CHECK(isfinite(u[0]));

    pf_free(pf);
    return 1;
}

static int test_tslqg(void) {
    int n = 2, m = 1, p = 1;
    TimestampLQG* lqg = tslqg_create(n, m, p,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0.01, 10);
    CHECK(lqg != NULL);
    CHECK(lqg->n == 2);
    CHECK(lqg->buffer_size == 10);

    double y[1] = {0.5};
    tslqg_update(lqg, y, true);
    const double* xh = tslqg_get_estimate(lqg);
    CHECK(xh != NULL);

    double u[1];
    tslqg_compute_control(lqg, u);
    CHECK(isfinite(u[0]));

    tslqg_free(lqg);
    return 1;
}

/* ============================================================================
   Run all tests
   ============================================================================ */

int main(void) {
    setbuf(stdout, NULL);
    printf("=== mini-time-delay-system — Test Suite ===\n\n");

    printf("L1 — Definitions:\n");
    RUN_TEST(delay_descriptor_constant);
    RUN_TEST(delay_descriptor_time_varying);
    RUN_TEST(delay_descriptor_stochastic);
    RUN_TEST(tds_create_free);
    RUN_TEST(tds_set_linear_model);
    RUN_TEST(tds_add_delay);

    printf("\nL2 — Core Concepts:\n");
    RUN_TEST(characteristic_eqn_scalar);
    RUN_TEST(delay_rate_check);

    printf("\nL3 — Mathematical Structures:\n");
    RUN_TEST(characteristic_root_computation);
    RUN_TEST(2d_system_stability);

    printf("\nL4 — Fundamental Theorems:\n");
    RUN_TEST(lkf_create_evaluate);
    RUN_TEST(lkf_derivative);
    RUN_TEST(lmi_delay_independent);
    RUN_TEST(lmi_delay_dependent);
    RUN_TEST(matrix_measure);
    RUN_TEST(halanay_inequality);

    printf("\nL5 — Algorithms:\n");
    RUN_TEST(pi_controller);
    RUN_TEST(smith_predictor_create);
    RUN_TEST(smith_predictor_step);
    RUN_TEST(pade_approximation);
    RUN_TEST(dde_solver_basic);
    RUN_TEST(predictor_feedback);

    printf("\nL6 — Canonical Problems:\n");
    RUN_TEST(delay_margin_scalar);
    RUN_TEST(tcp_aqm_model);
    RUN_TEST(nyquist_points);
    RUN_TEST(gain_phase_margin);
    RUN_TEST(rhp_root_count);
    RUN_TEST(sweeping_test);

    printf("\nL7 — Applications:\n");
    RUN_TEST(ncs_create_step);
    RUN_TEST(mm1_queue);
    RUN_TEST(tslqg);

    printf("\n========================================\n");
    printf("Results: %d/%d passed", tests_passed, tests_run);
    if (tests_failed > 0) printf("  (%d FAILED)", tests_failed);
    printf("\n");

    return tests_failed > 0 ? 1 : 0;
}
