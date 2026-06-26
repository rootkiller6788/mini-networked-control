#include "cps_security_core.h"
#include "cps_detection.h"
#include "cps_resilience.h"
#include "cps_gametheory.h"
#include "cps_watermarking.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int passed = 0, failed = 0;
#define CHECK(cond, msg) do { \
    if (cond) { passed++; printf("  PASS: %s\n", msg); } \
    else { failed++; printf("  FAIL: %s\n", msg); } \
    fflush(stdout); \
} while(0)

int main(void) {
    printf("=== CPS Security Test v2 ===\n"); fflush(stdout);
    
    /* L1 */
    printf("\n-- L1 Definitions --\n"); fflush(stdout);
    CHECK(strcmp(cps_attack_type_name(CPS_ATTACK_FDI), "False Data Injection") == 0, "FDI name");
    CHECK(strcmp(cps_security_state_name(CPS_SECURE_NORMAL), "Normal") == 0, "Normal state");
    CPSSecuritySystem* s = cps_security_create(2, 1, 1);
    CHECK(s != NULL, "System create");
    CHECK(s->plant != NULL, "Plant exists");
    CHECK(s->attack != NULL, "Attack model exists");
    CHECK(s->detector != NULL, "Detector exists");
    cps_security_free(s);
    
    /* L3 - Matrix ops */
    printf("\n-- L3 Math Structures --\n"); fflush(stdout);
    double A[4] = {1,2,3,4}, B[4] = {5,6,7,8}, C[4];
    cps_matrix_multiply(C, A, B, 2, 2, 2);
    CHECK(fabs(C[0]-19) < 1e-6 && fabs(C[3]-50) < 1e-6, "Matrix multiply");
    double I9[9] = {1,0,0,0,1,0,0,0,1};
    CHECK(cps_matrix_rank(I9, 3, 3, 1e-6) == 3, "Rank of identity");
    double v[3] = {3,4,0};
    CHECK(fabs(cps_vector_norm(v,3)-5) < 1e-6, "Vector norm");
    
    /* L4 */
    printf("\n-- L4 Fundamental Laws --\n"); fflush(stdout);
    s = cps_security_create(2, 1, 2);
    double A2[4] = {0.5, 1.0, 0.0, 0.5};
    double C2[4] = {1.0, 0.0, 0.0, 1.0};
    cps_set_state_matrix(s, A2, 2, 2);
    cps_set_output_matrix(s, C2, 2, 2);
    CHECK(cps_is_observable(s->plant), "Observability (C=I)");
    double g_obs[4];
    double logd = cps_observability_gramian(s->plant, g_obs, 5);
    CHECK(logd > -1e10, "Obs Gramian");
    cps_security_free(s);
    
    s = cps_security_create(2, 1, 1);
    double Actrl[4] = {0, 1, -1, 0};
    double Bctrl[2] = {0, 1};
    double Cdum[2] = {1, 0};
    cps_set_state_matrix(s, Actrl, 2, 2);
    cps_set_input_matrix(s, Bctrl, 2, 1);
    cps_set_output_matrix(s, Cdum, 1, 2);
    CHECK(cps_is_controllable(s->plant), "Controllability");
    cps_security_free(s);
    
    /* L5 */
    printf("\n-- L5 Algorithms --\n"); fflush(stdout);
    CPSDetector det;
    memset(&det, 0, sizeof(det));
    cps_detector_init(&det, CPS_DETECT_CHI2, 3.841);
    double res[2] = {0.1, -0.1}, cov[2] = {1, 1};
    double g = cps_chi2_test(&det, res, cov, 2);
    CHECK(g >= 0 && !det.alarm_active, "Chi2 normal");
    double big_res[2] = {5, 5};
    cps_chi2_test(&det, big_res, cov, 2);
    CHECK(det.alarm_active, "Chi2 alarm");
    cps_detector_free_internals(&det);
    
    cps_detector_init(&det, CPS_DETECT_CUSUM, 5);
    det.cusum_drift = 0.5;
    for (int i = 0; i < 20; i++) cps_cusum_update(&det, 1.0);
    CHECK(det.statistic > 0, "CUSUM accumulation");
    cps_detector_free_internals(&det);
    
    /* Game theory */
    double val, p, q;
    cps_solve_2x2_zerosum(1, 5, 0, 3, &val, &p, &q);
    CHECK(p >= 0 && p <= 1, "2x2 game solved");
    
    /* Watermark */
    CPSWatermark wm;
    cps_watermark_init(&wm, CPS_WATERMARK_GAUSSIAN, 1, 42);
    double wv = cps_watermark_next(&wm);
    CHECK(fabs(wv) < 10, "Watermark generation");
    cps_watermark_free(&wm);
    
    /* Attack config */
    s = cps_security_create(2, 1, 1);
    cps_attack_configure(s, CPS_ATTACK_FDI, CPS_TARGET_SENSOR, 5, 0);
    CHECK(!cps_attack_is_active(s), "Attack inactive");
    cps_attack_start(s);
    CHECK(cps_attack_is_active(s), "Attack active");
    cps_security_free(s);
    
    /* Sensor diversity */
    double Cs_div[6] = {1,0, 0,1, 1,1};
    double div = cps_sensor_diversity_score(Cs_div, 2, 3);
    CHECK(div >= 0 && div <= 1, "Sensor diversity");
    int min_s = cps_min_sensors_for_resilience(2);
    CHECK(min_s == 5, "2s+1 rule");
    
    printf("\n=== Results: %d/%d passed ===\n", passed, passed+failed);
    fflush(stdout);
    return (failed > 0) ? 1 : 0;
}
