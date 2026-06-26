#include "cps_security_core.h"
#include "cps_detection.h"
#include "cps_resilience.h"
#include "cps_gametheory.h"
#include "cps_watermarking.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int main(void) {
    srand((unsigned int)time(NULL));
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║   Mini Cyber-Physical Security — Full Stack Demo            ║\n");
    printf("║   Detection | Resilience | Game Theory | Watermarking       ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* 1. System observability demo */
    printf("--- 1. System Properties ---\n");
    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    double A[4] = {0.9, 0.1, 0.0, 0.9};
    double C[4] = {1.0, 0.0, 0.0, 1.0};
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_output_matrix(sys, C, 2, 2);
    printf("Observable: %s\n", cps_is_observable(sys->plant) ? "YES" : "NO");

    /* 2. Detection demo */
    printf("\n--- 2. Chi-Squared Detection ---\n");
    CPSDetector det;
    memset(&det, 0, sizeof(det));
    cps_detector_init(&det, CPS_DETECT_CHI2, 3.841);
    double res_ok[2] = {0.1, -0.2};
    double cov[2] = {0.05, 0.05};
    printf("Normal residual: g=%.3f, alarm=%s\n",
           cps_chi2_test(&det, res_ok, cov, 2),
           det.alarm_active ? "YES" : "no");
    double res_bad[2] = {4.0, 3.0};
    printf("Attack residual: g=%.3f, alarm=%s\n",
           cps_chi2_test(&det, res_bad, cov, 2),
           det.alarm_active ? "YES" : "no");

    /* 3. Game theory demo */
    printf("\n--- 3. Game-Theoretic Security ---\n");
    double val, p, q;
    cps_solve_2x2_zerosum(2, 8, 0, 5, &val, &p, &q);
    printf("2x2 Zero-sum: V=%.2f, def=(%.2f,%.2f), att=(%.2f,%.2f)\n",
           val, p, 1-p, q, 1-q);

    /* 4. Watermarking demo */
    printf("\n--- 4. Physical Watermarking ---\n");
    CPSWatermark wm;
    cps_watermark_init(&wm, CPS_WATERMARK_GAUSSIAN, 1.0, 42);
    printf("Watermark samples: ");
    for (int i = 0; i < 5; i++)
        printf("%.3f ", cps_watermark_next(&wm));
    printf("\nWatermark energy: %.3f\n", wm.energy);
    cps_watermark_free(&wm);

    /* 5. Attack simulation */
    printf("\n--- 5. Multi-Step Attack Simulation ---\n");
    double A5[4] = {1.0, 0.05, 0.0, 1.0};
    double B5[2] = {0.0, 0.05};
    double C5[2] = {1.0, 0.0};
    cps_set_state_matrix(sys, A5, 2, 2);
    cps_set_input_matrix(sys, B5, 2, 1);
    cps_set_output_matrix(sys, C5, 1, 2);
    cps_set_noise_covariances(sys, 0.0001, 0.001);
    double x0[2] = {1.0, 0.0};
    cps_set_initial_state(sys, x0);
    cps_attack_configure(sys, CPS_ATTACK_BIAS, CPS_TARGET_SENSOR, 1.0, 3.0);
    double u1[1] = {0.0};

    printf("Step | Time  | True x1 | Est x1  | Meas y  | Residual | Attack\n");
    printf("-----+-------+---------+---------+---------+----------+--------\n");
    for (int k = 0; k < 15; k++) {
        if (k == 5) cps_attack_start(sys);
        cps_step(sys, u1);
        printf("%4d | %5.2f | %7.4f | %7.4f | %7.4f | %8.4f | %s\n",
               k, sys->current_time, sys->true_state[0], sys->state[0],
               sys->measurement[0],
               sys->log_length > 0 ? sys->residual_log[sys->log_length-1] : 0.0,
               cps_attack_is_active(sys) ? "ON " : "off");
    }

    cps_detector_free_internals(&det);
    cps_security_free(sys);

    /* 6. Sensor diversity analysis */
    printf("\n--- 6. Sensor Diversity for Resilience ---\n");
    double Cs_div[6] = {1,0, 0,1, 1,1};
    double div = cps_sensor_diversity_score(Cs_div, 2, 3);
    printf("Sensor diversity score: %.4f (higher = better)\n", div);
    printf("Min sensors for 2 attacks: %d (2s+1 rule)\n",
           cps_min_sensors_for_resilience(2));

    printf("\n=== Demo complete ===\n");
    return 0;
}
