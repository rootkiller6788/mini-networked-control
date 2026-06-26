#include "cps_security_core.h"
#include "cps_watermarking.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 2: Replay Attack Defense via Physical Watermarking ===\n\n");

    CPSWatermark wm;
    cps_watermark_init(&wm, CPS_WATERMARK_GAUSSIAN, 0.5, 12345);
    cps_watermark_generate_sequence(&wm, 50);

    double A[4] = {1.0, 0.1, 0.0, 1.0};
    double B[2] = {0.0, 0.1};
    double C[2] = {1.0, 0.0};

    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_input_matrix(sys, B, 2, 1);
    cps_set_output_matrix(sys, C, 1, 2);
    double x0[2] = {0.5, 0.0};
    cps_set_initial_state(sys, x0);

    printf("Step | Level   | Watermark | Meas(no atk) | Meas(replay) | Detected\n");
    printf("-----+---------+-----------+--------------+---------------+----------\n");

    double normal_meas[30];
    double u[1] = {0.2};

    for (int k = 0; k < 30; k++) {
        double wm_val = cps_watermark_next(&wm);
        double u_wm[1] = {u[0] + wm_val * 0.1};
        cps_step(sys, u_wm);
        normal_meas[k] = sys->measurement[0];
    }

    for (int k = 0; k < 30; k++) {
        double replay_meas = normal_meas[(k + 15) % 30];
        double wm_resp = cps_watermark_next(&wm);
        double diff = replay_meas - (sys->state[0] + wm_resp * 0.1);
        int detected = (fabs(diff) > 0.5) ? 1 : 0;
        printf("%4d | %7.4f | %9.4f | %12.4f | %13.4f | %s\n",
               k, sys->state[0], wm_resp, normal_meas[k], replay_meas,
               detected ? "YES" : "no");
    }

    printf("\nPhysical watermarking enables replay attack detection.\n");
    printf("Without watermarking, replay would be indistinguishable from normal data.\n");

    cps_security_free(sys);
    cps_watermark_free(&wm);
    return 0;
}
