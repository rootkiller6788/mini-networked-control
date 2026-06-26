#include "cps_security_core.h"
#include "cps_detection.h"
#include "cps_resilience.h"
#include "cps_gametheory.h"
#include "cps_watermarking.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(void) {
    setvbuf(stdout, NULL, _IONBF, 0);
    printf("=== Mini Cyber-Physical Security Demo ===\n\n");

    /* Observability */
    printf("1. System Properties\n");
    CPSSecuritySystem* sys = cps_security_create(2, 1, 2);
    double A[4] = {0.9, 0.1, 0.0, 0.9};
    double C[4] = {1.0, 0.0, 0.0, 1.0};
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_output_matrix(sys, C, 2, 2);
    printf("   Observable: %s\n", cps_is_observable(sys->plant) ? "YES" : "NO");
    cps_security_free(sys);

    /* Detection */
    printf("\n2. Chi-Squared Detection\n");
    CPSDetector det;
    memset(&det, 0, sizeof(det));
    cps_detector_init(&det, CPS_DETECT_CHI2, 3.841);
    double res_ok[2] = {0.1, -0.2}, cov[2] = {0.05, 0.05};
    double g_ok = cps_chi2_test(&det, res_ok, cov, 2);
    printf("   Normal residual: g=%.3f, alarm=%s\n", g_ok, det.alarm_active ? "YES" : "no");
    double res_bad[2] = {4.0, 3.0};
    double g_bad = cps_chi2_test(&det, res_bad, cov, 2);
    printf("   Attack residual: g=%.3f, alarm=%s\n", g_bad, det.alarm_active ? "YES" : "no");
    cps_detector_free_internals(&det);

    /* Game Theory */
    printf("\n3. Game-Theoretic Security\n");
    double val, p, q;
    cps_solve_2x2_zerosum(2, 8, 0, 5, &val, &p, &q);
    printf("   2x2 Zero-sum: V=%.2f, def=(%.2f,%.2f), att=(%.2f,%.2f)\n", val, p, 1-p, q, 1-q);

    /* Watermarking */
    printf("\n4. Physical Watermarking\n");
    CPSWatermark wm;
    cps_watermark_init(&wm, CPS_WATERMARK_GAUSSIAN, 1.0, 42);
    printf("   Watermark samples: ");
    for (int i = 0; i < 5; i++) printf("%.3f ", cps_watermark_next(&wm));
    printf("\n   Watermark energy: %.3f\n", wm.energy);
    cps_watermark_free(&wm);

    /* Attack Simulation */
    printf("\n5. FDI Attack Simulation\n");
    sys = cps_security_create(2, 1, 1);
    double A5[4] = {1.0, 0.05, 0.0, 1.0};
    double B5[2] = {0.0, 0.05};
    double C5[2] = {1.0, 0.0};
    double x0[2] = {1.0, 0.0};
    cps_set_state_matrix(sys, A5, 2, 2);
    cps_set_input_matrix(sys, B5, 2, 1);
    cps_set_output_matrix(sys, C5, 1, 2);
    cps_set_noise_covariances(sys, 0.0001, 0.001);
    cps_set_initial_state(sys, x0);
    cps_attack_configure(sys, CPS_ATTACK_BIAS, CPS_TARGET_SENSOR, 1.0, 3.0);
    cps_attack_generate_signal(sys);
    double u1[1] = {0.0};
    printf("   Step | Time  | True x1 | Est x1  | Meas y  | Residual | Attack\n");
    printf("   -----+-------+---------+---------+---------+----------+-------\n");
    for (int k = 0; k < 10; k++) {
        if (k == 4) cps_attack_start(sys);
        cps_step(sys, u1);
        printf("   %4d | %5.2f | %7.4f | %7.4f | %7.4f | %8.4f | %s\n",
               k, sys->current_time, sys->true_state[0], sys->state[0],
               sys->measurement[0],
               sys->log_length > 0 ? sys->residual_log[sys->log_length-1] : 0.0,
               cps_attack_is_active(sys) ? "ON " : "off");
    }
    cps_security_free(sys);

    /* Sensor Diversity */
    printf("\n6. Sensor Diversity for Resilience\n");
    double Cs_div[6] = {1,0, 0,1, 1,1};
    double div = cps_sensor_diversity_score(Cs_div, 2, 3);
    printf("   Sensor diversity score: %.4f\n", div);
    printf("   Min sensors for s=2 attacks: %d (2s+1 rule)\n", cps_min_sensors_for_resilience(2));

    printf("\n=== Demo complete ===\n");
    return 0;
}
