#include "time_delay_system.h"
#include "delay_stability.h"
#include "dde_solver.h"
#include "lyapunov_krasovskii.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ============================================================================
 * Example 3: DDE Stability Analysis — Delay Margin Sweep
 *
 * End-to-end demonstration of:
 *   1. Delay margin computation for a scalar DDE
 *   2. Numerical DDE integration via method of steps
 *   3. Lyapunov-Krasovskii functional evaluation
 *   4. Stability classification and visualization
 *
 * System: ẋ(t) = -a x(t) - a_d x(t-τ)
 * Classic problem: find τ* such that system is stable iff τ < τ*
 * ============================================================================ */

/* Example nonlinear history: step function */
static void step_history(double t, int n, double* x) {
    (void)n;
    if (t < -0.5) {
        x[0] = 0.5;  /* Old value */
    } else {
        x[0] = 1.0;  /* Recent value */
    }
}

int main(void) {
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  DDE Stability — Delay Margin Analysis       ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* --- Part 1: Delay Margin for Multiple Systems --- */
    printf("=== Part 1: Delay Margin Computation ===\n\n");

    struct { double a; double ad; const char* desc; } cases[] = {
        {-1.0, -0.3, "Stable, |ad| < |a| → delay-independent"},
        {-0.5, -2.0, "Delay-sensitive, |ad| > |a|"},
        {-2.0,  0.5, "Unstable without delay"},
        {-0.1, -3.0, "Very delay-sensitive"}
    };
    int n_cases = 4;

    printf("%-45s %12s %12s %20s\n",
           "System", "a", "a_d", "τ* (s)");
    printf("────────────────────────────────────────────────────────────────\n");

    for (int c = 0; c < n_cases; c++) {
        TimeDelaySystem* sys = tds_create("case", 1, 1, 1);
        double A[1] = {cases[c].a};
        double Ad[1] = {cases[c].ad};
        double B[1] = {1.0};
        double C[1] = {1.0};
        tds_set_linear_model(sys, A, Ad, B, C);
        tds_add_delay(sys, DELAY_CONSTANT, 0.5, 0.0, 0.5);

        double tau_star = delay_margin_frequency_sweep(sys);

        printf("%-45s %12.3f %12.3f ",
               cases[c].desc, cases[c].a, cases[c].ad);
        if (tau_star == INFINITY)
            printf("%20s\n", "∞ (always stable)");
        else if (tau_star <= 0)
            printf("%20s\n", "0 (always unstable)");
        else
            printf("%18.6f\n", tau_star);

        tds_free(sys);
    }
    printf("\n");

    /* --- Part 2: Numerical DDE Integration --- */
    printf("=== Part 2: DDE Integration (Method of Steps + RK4) ===\n\n");

    TimeDelaySystem* sys_dde = tds_create("integrate", 1, 1, 1);
    double A_dde[1] = {-1.0};
    double Ad_dde[1] = {-0.5};
    double B_dde[1] = {1.0};
    double C_dde[1] = {1.0};
    tds_set_linear_model(sys_dde, A_dde, Ad_dde, B_dde, C_dde);
    tds_add_delay(sys_dde, DELAY_CONSTANT, 0.5, 0.5, 0.5);
    tds_set_history(sys_dde, step_history);

    /* Delay margin for this system */
    double tau_star = delay_margin_frequency_sweep(sys_dde);
    printf("System: ẋ = -x - 0.5 x(t-τ), τ=0.5\n");
    printf("Delay margin τ* = %.4f s\n", tau_star);
    printf("At τ=0.5: %s\n\n",
           (0.5 < tau_star) ? "STABLE ✓" : "UNSTABLE ⚠");

    /* Integrate using RK4 */
    DDESolverConfig cfg = dde_config_rk4(0.01, 3.0);
    cfg.track_discontinuities = true;
    double x0[1] = {1.0};

    DDESolution* sol = dde_solve(sys_dde, &cfg, step_history, x0);
    if (sol && sol->success) {
        printf("Integration: %d steps, final t=%.2f\n",
               sol->n_steps, sol->final_time);

        /* Print trajectory at selected points */
        printf("\nTrajectory samples:\n");
        printf("  t=0.0:  x=%.6f\n", sol->x[0]);
        for (int i = 20; i < sol->n_steps; i += sol->n_steps / 6) {
            if (i < sol->n_steps) {
                printf("  t=%.2f: x=%.6f\n",
                       sol->t[i], sol->x[(size_t)i * sol->n_states]);
            }
        }
        size_t last = (size_t)(sol->n_steps - 1);
        printf("  t=%.2f: x=%.6f\n",
               sol->t[last], sol->x[last * sol->n_states]);
        printf("\n");
    } else {
        printf("Integration failed!\n");
    }

    /* --- Part 3: LK Functional Evaluation --- */
    printf("=== Part 3: Lyapunov-Krasovskii Functional ===\n\n");

    LKFunctional* lkf = lkf_create(1);
    double P[1] = {1.0}, Q[1] = {2.0}, R[1] = {0.5};
    lkf_set_P(lkf, P);
    lkf_set_Q(lkf, Q);
    lkf_set_R(lkf, R);

    bool is_pd = lkf_is_positive_definite(lkf);
    printf("P, Q, R positive definite: %s\n", is_pd ? "YES ✓" : "NO");

    /* Evaluate V and dV/dt at t=0 */
    if (sol && sol->n_steps > 1) {
        const double* x0_state = dde_solution_state_at(sol, 0);
        const double* x1_state = dde_solution_state_at(sol, 1);
        /* Approximate xdot */
        double xdot_0 = (x1_state[0] - x0_state[0]) / (sol->t[1] - sol->t[0]);
        double xdot_arr[1] = {xdot_0};

        /* Approximate delayed state at t=0: use history */
        double xd0[1] = {0.5};  /* Step history value at t=-τ=-0.5 */

        double V0 = lkf_evaluate(lkf, x0_state, NULL, NULL, 0, 0.5, 0.01);
        double dV0 = lkf_derivative(lkf, sys_dde, x0_state, xd0,
                                     xdot_arr, NULL, 0, 0.5, 0.01);

        printf("V(x_0) = %.6f\n", V0);
        printf("dV/dt(x_0) = %.6f\n", dV0);
        printf("Decay condition dV/dt < 0: %s\n", (dV0 < 0) ? "YES ✓" : "NO");

        double eps = lkf_decay_rate(lkf, sys_dde, x0_state, xd0, 0.5);
        printf("Decay rate ε = %.6f\n\n", eps);
    }

    /* --- Part 4: Characteristic Roots --- */
    printf("=== Part 4: Characteristic Roots ===\n\n");

    tds_compute_characteristic_roots(sys_dde, 12);
    printf("Rightmost characteristic roots (n=%d):\n", sys_dde->n_roots);
    int show = sys_dde->n_roots < 8 ? sys_dde->n_roots : 8;
    for (int i = 0; i < show; i++) {
        printf("  λ_%d = %.6f + j %.6f  ", i,
               sys_dde->roots_real[i], sys_dde->roots_imag[i]);
        if (sys_dde->roots_real[i] < -1e-6) printf("(stable)");
        else if (fabs(sys_dde->roots_real[i]) < 1e-6) printf("(marginally stable)");
        else printf("(UNSTABLE)");
        printf("\n");
    }

    double abscissa = tds_spectral_abscissa(sys_dde);
    printf("\nSpectral abscissa α = %.6f\n", abscissa);
    printf("Exponential decay rate = %.6f\n", delay_exponential_decay_rate(sys_dde));
    printf("Stability class: %d\n\n", sys_dde->stability_class);

    /* --- Cleanup --- */
    lkf_free(lkf);
    dde_solution_free(sol);
    tds_free(sys_dde);

    printf("Done.\n");
    return 0;
}
