#include "cps_security_core.h"
#include "cps_detection.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 1: False Data Injection (FDI) Attack on CPS ===\n\n");

    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    double A[4] = {0.98, 0.10, 0.00, 0.95};
    double B[2] = {0.00, 0.10};
    double C[2] = {1.00, 0.00};
    double x0[2] = {1.0, 0.0};
    cps_set_state_matrix(sys, A, 2, 2);
    cps_set_input_matrix(sys, B, 2, 1);
    cps_set_output_matrix(sys, C, 1, 2);
    cps_set_initial_state(sys, x0);
    cps_set_noise_covariances(sys, 0.001, 0.01);

    cps_attack_configure(sys, CPS_ATTACK_FDI, CPS_TARGET_SENSOR, 2.0, 5.0);
    cps_attack_generate_signal(sys);
    cps_detector_init(sys->detector, CPS_DETECT_CHI2, 3.841);

    double u[1] = {0.0};
    printf("Time  | State[0] | State[1] | Measurement | Residual | Alarm\n");
    printf("------+----------+----------+-------------+----------+------\n");

    for (int k = 0; k < 30; k++) {
        if (k == 10) cps_attack_start(sys);
        cps_step(sys, u);
        double r = sys->residual_log[sys->log_length - 1];
        printf("%5.1f | %8.4f | %8.4f | %11.4f | %8.4f | %s\n",
               sys->current_time, sys->state[0], sys->state[1],
               sys->measurement[0], r,
               cps_attack_is_active(sys) ? "ATTACK" : "normal");
    }

    printf("\nAttack was: %s\n", cps_attack_is_active(sys) ? "Active" : "Inactive");
    printf("Security state: %s\n", cps_security_state_name(sys->security_state));
    cps_security_free(sys);
    return 0;
}
