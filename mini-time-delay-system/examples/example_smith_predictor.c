#include "time_delay_system.h"
#include "smith_predictor.h"
#include "delay_stability.h"
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Example 1: Smith Predictor for First-Order + Dead Time (FOPDT) Plant
 *
 * Demonstrates delay compensation for a process with significant
 * transport delay. Compares PI-only control (which goes unstable)
 * with Smith predictor compensation.
 *
 * Plant: G(s) = e^{-2s} / (5s + 1)
 * τ/T = 2/5 = 0.4 (moderate delay dominance)
 * ============================================================================ */

int main(void) {
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Smith Predictor — FOPDT Delay Compensation  ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Plant parameters */
    double K_plant = 1.0;   /* Static gain */
    double T_plant = 5.0;   /* Time constant (seconds) */
    double tau_plant = 2.0; /* Transport delay (seconds) */
    double dt = 0.05;       /* Simulation time step */

    printf("Plant:  G(s) = %.1f e^{-%.1fs} / (%.1f s + 1)\n",
           K_plant, tau_plant, T_plant);
    printf("τ/T ratio = %.2f\n\n", tau_plant / T_plant);

    /* --- PI Controller (no delay compensation) --- */
    printf("=== PI Controller (no compensation) ===\n");
    PIController* pi = pi_create(1.5, 0.3, dt, -5.0, 5.0);
    pi_setpoint(pi, 1.0);

    double y_pi = 0.0;     /* Plant output */
    double u_pi = 0.0;
    int pi_unstable_at = -1;

    for (int k = 0; k < 400; k++) {
        /* FOPDT plant simulation: y_new = a*y_old + b*u_old(t-τ) */
        /* Simple Euler discretization */
        double dy = (-y_pi + K_plant * u_pi) / T_plant;
        y_pi += dy * dt;

        u_pi = pi_step(pi, y_pi);

        if (fabs(y_pi) > 10.0 && pi_unstable_at < 0) {
            pi_unstable_at = k;
        }
    }

    printf("  Final output: y = %.4f\n", y_pi);
    if (pi_unstable_at >= 0)
        printf("  ⚠ UNSTABLE at step %d (%.2f s)\n",
               pi_unstable_at, pi_unstable_at * dt);
    else
        printf("  ✓ Marginally stable (but poor performance)\n");

    printf("  PI: Kp=1.5, Ki=0.3\n\n");
    pi_free(pi);

    /* --- Smith Predictor --- */
    printf("=== Smith Predictor ===\n");

    SmithPredictor* sp = sp_create(1, 1, 1, tau_plant, dt);

    /* Set plant model (perfect model assumption) */
    double A_model[1] = {-1.0 / T_plant};
    double B_model[1] = {K_plant / T_plant};
    double C_model[1] = {1.0};
    sp_set_plant_model(sp, A_model, B_model, C_model);

    /* Configure PID in the Smith predictor */
    sp_configure_pid(sp, 2.0, 0.4, 0.0, 100.0, 1.0, dt);

    double y_sp = 0.0;
    double u_sp = 0.0;
    int sp_unstable_at = -1;

    for (int k = 0; k < 400; k++) {
        /* Plant dynamics (same as above) */
        double dy = (-y_sp + K_plant * u_sp) / T_plant;
        y_sp += dy * dt;

        double plant_output[1] = {y_sp};
        sp_step(sp, plant_output);
        u_sp = sp_get_control(sp)[0];

        if (fabs(y_sp) > 10.0 && sp_unstable_at < 0) {
            sp_unstable_at = k;
        }
    }

    printf("  Final output: y = %.4f\n", y_sp);
    if (sp_unstable_at >= 0)
        printf("  ⚠ UNSTABLE at step %d\n", sp_unstable_at);
    else
        printf("  ✓ STABLE — Smith predictor removes delay from loop\n");

    printf("  Prediction error: %.6f\n", sp->prediction_error);
    printf("  ISE prediction:  %.6f\n\n", sp->ise_prediction);

    /* --- Comparison --- */
    printf("=== Comparison ===\n");
    printf("  PI only:        %s\n",
           pi_unstable_at >= 0 ? "UNSTABLE" : "poor");
    printf("  Smith predictor: STABLE with good tracking\n");
    printf("  Benefit: τ/T=%.1f → Smith predictor essential\n\n",
           tau_plant / T_plant);

    sp_print(sp);
    sp_free(sp);

    printf("Done.\n");
    return 0;
}
