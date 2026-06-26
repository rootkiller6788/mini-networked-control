#include "cps_security_core.h"
#include "cps_detection.h"
#include "cps_resilience.h"
#include "cps_watermarking.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CPS_EPS 1e-12

/* ============================================================================
 * CPS Security Applications (L7: Applications)
 *
 * Real-world CPS security scenarios demonstrating the integration
 * of detection, resilience, and watermarking subsystems.
 *
 * Applications covered:
 *   1. Smart Grid FDI Attack Detection (Mo & Sinopoli, 2010)
 *   2. Autonomous Vehicle GPS Spoofing Detection
 *   3. Industrial Control System (ICS) Replay Attack Defense
 *   4. Water Distribution SCADA Security
 * ============================================================================ */

/* ============================================================================
 * Application 1: Smart Grid State Estimation Security
 *
 * Scenario: An attacker injects false data into power flow measurements
 * to deceive the state estimator. The chi-squared detector monitors
 * the measurement residual for anomalies.
 *
 * System: DC power flow model (linear)
 *   x = voltage angles (n buses)
 *   y = power flow measurements (p measurements)
 *   C = measurement matrix (power flow Jacobian)
 *
 * Reference: Liu, Ning & Reiter (2011), "False Data Injection Attacks
 *            against State Estimation in Electric Power Grids"
 *            IEEE Transactions on Information Forensics and Security
 * ============================================================================ */

typedef struct {
    int n_buses;
    int n_measurements;
    double* bus_angles;         /* True voltage angles */
    double* bus_angles_est;     /* Estimated angles */
    double* measurements;       /* Power measurements */
    double* measurement_matrix; /* p x n Jacobian */
    double* attack_vector;      /* FDI attack vector */
    int attacked_measurement;
    double alarm_threshold;
    int alarm_raised;
    double detection_time;
} SmartGridSecurity;

SmartGridSecurity* smart_grid_security_create(int n_buses, int n_meas) {
    SmartGridSecurity* sg = (SmartGridSecurity*)calloc(1,
        sizeof(SmartGridSecurity));
    sg->n_buses = n_buses;
    sg->n_measurements = n_meas;
    sg->bus_angles = (double*)calloc((size_t)n_buses, sizeof(double));
    sg->bus_angles_est = (double*)calloc((size_t)n_buses, sizeof(double));
    sg->measurements = (double*)calloc((size_t)n_meas, sizeof(double));
    sg->measurement_matrix = (double*)calloc(
        (size_t)(n_meas * n_buses), sizeof(double));
    sg->attack_vector = (double*)calloc((size_t)n_meas, sizeof(double));
    sg->alarm_threshold = 3.841; /* chi2_0.95(1) */
    return sg;
}

void smart_grid_security_free(SmartGridSecurity* sg) {
    if (!sg) return;
    free(sg->bus_angles); free(sg->bus_angles_est);
    free(sg->measurements); free(sg->measurement_matrix);
    free(sg->attack_vector);
    free(sg);
}

/* Set up IEEE 5-bus system measurement matrix */
void smart_grid_setup_ieee5(SmartGridSecurity* sg) {
    if (!sg || sg->n_buses < 5) return;
    /* Simple 5-bus topology:
     * Bus 1 (reference, angle=0) -> Bus 2, Bus 3
     * Bus 2 -> Bus 4, Bus 5
     * Measurements: P12, P13, P24, P25, P1(ref) */
    int n = sg->n_buses;
    int p = sg->n_measurements;
    /* Set up DC power flow measurement matrix (simplified) */
    double C_example[] = {
        /* P12 */  1.0, -1.0,  0.0,  0.0,  0.0,
        /* P13 */  1.0,  0.0, -1.0,  0.0,  0.0,
        /* P24 */  0.0,  1.0,  0.0, -1.0,  0.0,
        /* P25 */  0.0,  1.0,  0.0,  0.0, -1.0,
        /* Pref*/  1.0,  0.0,  0.0,  0.0,  0.0,
    };
    for (int i = 0; i < p * n && i < 25; i++)
        sg->measurement_matrix[i] = C_example[i];
}

/* Configure FDI attack targeting measurement idx */
void smart_grid_configure_fdi_attack(SmartGridSecurity* sg,
                                      int measurement_idx,
                                      double attack_magnitude) {
    if (!sg || measurement_idx < 0
        || measurement_idx >= sg->n_measurements) return;
    sg->attacked_measurement = measurement_idx;
    for (int i = 0; i < sg->n_measurements; i++)
        sg->attack_vector[i] = (i == measurement_idx)
            ? attack_magnitude : 0.0;
}

/* Run state estimation and check for FDI attack */
int smart_grid_detect_fdi(SmartGridSecurity* sg) {
    if (!sg) return 0;
    int n = sg->n_buses;
    int p = sg->n_measurements;

    /* True measurements (with possible attack) */
    for (int i = 0; i < p; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++)
            s += sg->measurement_matrix[i * n + j]
                 * sg->bus_angles[j];
        sg->measurements[i] = s + sg->attack_vector[i]
            + 0.01 * ((double)rand() / RAND_MAX - 0.5);
    }

    /* State estimation: min ||y - C*x||^2
     * x_est = (C'C)^{-1} C' y */
    double* CtC = (double*)calloc((size_t)(n * n), sizeof(double));
    double* CtY = (double*)calloc((size_t)n, sizeof(double));

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++)
                CtC[j * n + k] += sg->measurement_matrix[i * n + j]
                    * sg->measurement_matrix[i * n + k];
            CtY[j] += sg->measurement_matrix[i * n + j]
                * sg->measurements[i];
        }
    }

    /* Solve via diagonal approximation */
    for (int i = 0; i < n; i++)
        sg->bus_angles_est[i] = (fabs(CtC[i*n+i]) > CPS_EPS)
            ? CtY[i] / CtC[i*n+i] : 0.0;

    /* Compute residual: r = y - C*x_est */
    double residual_norm = 0.0;
    for (int i = 0; i < p; i++) {
        double y_hat = 0.0;
        for (int j = 0; j < n; j++)
            y_hat += sg->measurement_matrix[i * n + j]
                     * sg->bus_angles_est[j];
        double r = sg->measurements[i] - y_hat;
        residual_norm += r * r;
    }

    /* Chi-squared test */
    sg->alarm_raised = (residual_norm > sg->alarm_threshold) ? 1 : 0;

    free(CtC); free(CtY);
    return sg->alarm_raised;
}

void smart_grid_print_report(SmartGridSecurity* sg) {
    if (!sg) return;
    printf("=== Smart Grid Security Report ===\n");
    printf("Buses: %d, Measurements: %d\n",
           sg->n_buses, sg->n_measurements);
    printf("Attacked measurement: %d, Magnitude: %.2f\n",
           sg->attacked_measurement,
           sg->attack_vector[sg->attacked_measurement]);
    printf("Alarm raised: %s\n",
           sg->alarm_raised ? "YES (FDI detected!)"
                            : "NO (Attack undetected)");
    printf("Estimated bus angles:\n");
    for (int i = 0; i < sg->n_buses && i < 5; i++)
        printf("  Bus %d: true=%.4f, est=%.4f\n", i + 1,
               sg->bus_angles[i], sg->bus_angles_est[i]);
}

/* ============================================================================
 * Application 2: Autonomous Vehicle GPS Spoofing Detection
 *
 * Scenario: An attacker spoofs GPS signals to mislead an autonomous
 * vehicle's navigation system. The vehicle uses a Kalman filter
 * fusing GPS, IMU, and wheel odometry. GPS spoofing creates a
 * discrepancy between GPS and dead-reckoning position estimates.
 *
 * Detection: The innovation (residual) of the Kalman filter grows
 * when GPS is spoofed because the GPS measurement no longer matches
 * the motion model prediction from IMU + odometry.
 * ============================================================================ */

typedef struct {
    double pos_x;         /* True position */
    double pos_y;
    double vel_x;         /* True velocity */
    double vel_y;
    double pos_x_est;     /* Estimated position */
    double pos_y_est;
    double vel_x_est;     /* Estimated velocity */
    double vel_y_est;
    double gps_x;         /* GPS measurement (possibly spoofed) */
    double gps_y;
    double imu_accel_x;   /* IMU acceleration */
    double imu_accel_y;
    double wheel_speed;   /* Wheel odometry */
    double spoof_offset_x;/* GPS spoofing offset */
    double spoof_offset_y;
    int spoof_active;
    double residual_magnitude;
    double detection_threshold;
    int alarm_raised;
} AVSecurity;

AVSecurity* av_security_create(void) {
    AVSecurity* av = (AVSecurity*)calloc(1, sizeof(AVSecurity));
    av->detection_threshold = 5.0; /* 5 meter residual */
    return av;
}

void av_security_free(AVSecurity* av) { free(av); }

void av_security_set_spoof(AVSecurity* av, double offset_x,
                            double offset_y) {
    if (!av) return;
    av->spoof_offset_x = offset_x;
    av->spoof_offset_y = offset_y;
    av->spoof_active = 1;
}

void av_security_step(AVSecurity* av, double dt,
                       double accel_x, double accel_y,
                       double steering) {
    if (!av) return;

    /* True dynamics (bicycle model simplified) */
    av->vel_x += accel_x * dt;
    av->vel_y += accel_y * dt;
    av->pos_x += av->vel_x * dt;
    av->pos_y += av->vel_y * dt;

    /* GPS measurement (with possible spoofing) */
    av->gps_x = av->pos_x + (av->spoof_active ? av->spoof_offset_x : 0.0)
                + 0.5 * ((double)rand() / RAND_MAX - 0.5);
    av->gps_y = av->pos_y + (av->spoof_active ? av->spoof_offset_y : 0.0)
                + 0.5 * ((double)rand() / RAND_MAX - 0.5);

    /* Dead-reckoning prediction (IMU + odometry) */
    double pred_x = av->pos_x_est + av->vel_x_est * dt;
    double pred_y = av->pos_y_est + av->vel_y_est * dt;
    av->vel_x_est += accel_x * dt;
    av->vel_y_est += accel_y * dt;

    /* Kalman update: fuse GPS with prediction */
    double K = 0.3; /* Kalman gain */
    av->pos_x_est = pred_x + K * (av->gps_x - pred_x);
    av->pos_y_est = pred_y + K * (av->gps_y - pred_y);

    /* Residual = GPS - prediction */
    double rx = av->gps_x - pred_x;
    double ry = av->gps_y - pred_y;
    av->residual_magnitude = sqrt(rx * rx + ry * ry);

    /* Detection */
    av->alarm_raised = (av->residual_magnitude
                        > av->detection_threshold) ? 1 : 0;
}

void av_security_print_status(AVSecurity* av) {
    if (!av) return;
    printf("=== Autonomous Vehicle Security ===\n");
    printf("True position:     (%.2f, %.2f)\n",
           av->pos_x, av->pos_y);
    printf("Estimated position:(%.2f, %.2f)\n",
           av->pos_x_est, av->pos_y_est);
    printf("GPS measurement:   (%.2f, %.2f)\n",
           av->gps_x, av->gps_y);
    printf("Residual magnitude: %.4f m\n",
           av->residual_magnitude);
    printf("Spoof active: %s, Alarm: %s\n",
           av->spoof_active ? "YES" : "NO",
           av->alarm_raised ? "SPOOFING DETECTED!" : "Normal");
}

/* ============================================================================
 * Application 3: ICS Replay Attack Defense with Watermarking
 *
 * Scenario: A Stuxnet-class attacker records sensor measurements
 * during normal operation and replays them during the attack to
 * hide the physical damage being caused. Physical watermarking
 * defeats this by injecting a known signal that the replay cannot
 * replicate.
 * ============================================================================ */

typedef struct {
    double tank_level;         /* True physical process value */
    double tank_level_meas;    /* Measured value */
    double replay_buffer[100]; /* Attacker's recording */
    int replay_pos;
    int replay_active;
    double watermark_signal;   /* Physical watermark value */
    double expected_response;  /* Expected measurement response */
    double detection_stat;
    double threshold;
    int alarm_raised;
    double valve_position;     /* Control input */
} ICSReplayDefense;

ICSReplayDefense* ics_replay_defense_create(void) {
    ICSReplayDefense* ics = (ICSReplayDefense*)calloc(1,
        sizeof(ICSReplayDefense));
    ics->tank_level = 10.0;
    ics->threshold = 2.0;
    ics->valve_position = 0.5;
    return ics;
}

void ics_replay_defense_free(ICSReplayDefense* ics) { free(ics); }

void ics_replay_record(ICSReplayDefense* ics, double measurement) {
    if (!ics) return;
    ics->replay_buffer[ics->replay_pos % 100] = measurement;
    ics->replay_pos++;
}

void ics_replay_attack_start(ICSReplayDefense* ics) {
    if (!ics) return;
    ics->replay_active = 1;
    ics->replay_pos = 0;
}

void ics_step_with_watermark(ICSReplayDefense* ics, double dt,
                              double inflow, double setpoint) {
    if (!ics) return;
    double A = 0.02; /* Tank cross-section */

    /* Generate physical watermark */
    double wm = 0.1 * sin(2.0 * 3.14159 * 0.5 * dt * 100.0);
    ics->watermark_signal = wm;

    /* Control with watermark: valve = PI(setpoint - level) + watermark */
    double error = setpoint - ics->tank_level;
    ics->valve_position = 0.3 * error + wm;
    if (ics->valve_position < 0.0) ics->valve_position = 0.0;
    if (ics->valve_position > 1.0) ics->valve_position = 1.0;

    /* True process: d(level)/dt = inflow - outflow*sqrt(level) * valve */
    double outflow = ics->valve_position * sqrt(ics->tank_level);
    ics->tank_level += (inflow - outflow) * dt / A;
    if (ics->tank_level < 0.0) ics->tank_level = 0.0;

    /* Measurement */
    double true_meas = ics->tank_level
        + 0.02 * ((double)rand() / RAND_MAX - 0.5);
    ics->expected_response = ics->tank_level + wm * 0.1;

    if (ics->replay_active) {
        ics->tank_level_meas = ics->replay_buffer[
            ics->replay_pos % 100];
        ics->replay_pos++;
    } else {
        ics->tank_level_meas = true_meas;
        ics_replay_record(ics, true_meas);
    }

    /* Watermark detection: compare measurement to expected */
    double diff = ics->tank_level_meas - ics->expected_response;
    ics->detection_stat = fabs(diff);
    ics->alarm_raised = (ics->detection_stat > ics->threshold) ? 1 : 0;
}

void ics_replay_print_status(ICSReplayDefense* ics) {
    if (!ics) return;
    printf("=== ICS Replay Attack Defense ===\n");
    printf("Tank level:    true=%.4f, meas=%.4f\n",
           ics->tank_level, ics->tank_level_meas);
    printf("Watermark:     signal=%.4f, expected_resp=%.4f\n",
           ics->watermark_signal, ics->expected_response);
    printf("Detection stat:%.4f (threshold=%.2f)\n",
           ics->detection_stat, ics->threshold);
    printf("Replay active: %s, Alarm: %s\n",
           ics->replay_active ? "YES" : "NO",
           ics->alarm_raised ? "REPLAY DETECTED!" : "Normal");
}

/* ============================================================================
 * Application 4: Water Distribution SCADA Security
 *
 * Scenario: An attacker compromises the SCADA system of a water
 * distribution network and manipulates chlorine residual sensor
 * readings to hide contamination. The CUSUM detector identifies
 * subtle, sustained deviations in residual patterns.
 * ============================================================================ */

typedef struct {
    double chlorine_level;     /* True chlorine concentration (mg/L) */
    double chlorine_meas;      /* SCADA measurement */
    double attack_bias;        /* Attacker's bias on measurement */
    int attack_active;
    double cusum_accumulator;
    double cusum_threshold;
    double cusum_drift;
    int alarm_raised;
    double flow_rate;
    double chlorine_dosing;    /* Chlorine injection rate */
} WaterSCADASecurity;

WaterSCADASecurity* water_scada_create(void) {
    WaterSCADASecurity* ws = (WaterSCADASecurity*)calloc(1,
        sizeof(WaterSCADASecurity));
    ws->chlorine_level = 2.0; /* mg/L (safe level) */
    ws->cusum_threshold = 5.0;
    ws->cusum_drift = 0.2;
    ws->flow_rate = 100.0; /* L/s */
    ws->chlorine_dosing = 0.5; /* mg/L */
    return ws;
}

void water_scada_free(WaterSCADASecurity* ws) { free(ws); }

void water_scada_attack_start(WaterSCADASecurity* ws,
                               double bias) {
    if (!ws) return;
    ws->attack_active = 1;
    ws->attack_bias = bias;
}

void water_scada_step(WaterSCADASecurity* ws, double dt) {
    if (!ws) return;

    /* Chlorine decay: dC/dt = -k*C + dosing/volume
     * k=0.01 (decay rate for chlorine) */
    double decay = -0.01 * ws->chlorine_level;
    double dosing = ws->chlorine_dosing / (ws->flow_rate * dt + 1.0);
    ws->chlorine_level += (decay + dosing) * dt;
    if (ws->chlorine_level < 0.0) ws->chlorine_level = 0.0;

    /* True measurement with noise */
    double true_meas = ws->chlorine_level
        + 0.03 * ((double)rand() / RAND_MAX - 0.5);

    /* Attacked measurement */
    ws->chlorine_meas = true_meas
        + (ws->attack_active ? ws->attack_bias : 0.0);

    /* Expected level (from model) */
    double expected = ws->chlorine_level; /* Ideal model */

    /* CUSUM update */
    double residual = ws->chlorine_meas - expected;
    ws->cusum_accumulator += residual - ws->cusum_drift;
    if (ws->cusum_accumulator < 0.0) ws->cusum_accumulator = 0.0;

    ws->alarm_raised = (ws->cusum_accumulator
                        > ws->cusum_threshold) ? 1 : 0;
}

void water_scada_print_status(WaterSCADASecurity* ws) {
    if (!ws) return;
    printf("=== Water Distribution SCADA Security ===\n");
    printf("Chlorine: true=%.4f mg/L, meas=%.4f mg/L\n",
           ws->chlorine_level, ws->chlorine_meas);
    printf("CUSUM:    accumulator=%.4f (threshold=%.2f, drift=%.2f)\n",
           ws->cusum_accumulator,
           ws->cusum_threshold, ws->cusum_drift);
    printf("Attack: %s, Alarm: %s\n",
           ws->attack_active ? "ACTIVE" : "None",
           ws->alarm_raised
               ? "CONTAMINATION DETECTED! Shut down."
               : "Normal operation");
}

/* ============================================================================
 * Integrated CPS Security Demo (L7: Applications)
 *
 * Combines all four application scenarios into a single demonstration
 * showing the complete CPS security lifecycle.
 * ============================================================================ */

void cps_applications_demo_all(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║     CPS Security Applications Demonstration          ║\n");
    printf("║   Smart Grid | Autonomous Vehicle | ICS | SCADA      ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n");

    /* 1. Smart Grid FDI */
    printf("\n--- 1. Smart Grid FDI Attack Detection ---\n");
    SmartGridSecurity* sg = smart_grid_security_create(5, 5);
    smart_grid_setup_ieee5(sg);
    sg->bus_angles[0] = 0.0;
    sg->bus_angles[1] = -0.05;
    sg->bus_angles[2] = -0.03;
    sg->bus_angles[3] = -0.08;
    sg->bus_angles[4] = -0.06;

    smart_grid_configure_fdi_attack(sg, 2, 5.0);
    smart_grid_detect_fdi(sg);
    smart_grid_print_report(sg);
    smart_grid_security_free(sg);

    /* 2. Autonomous Vehicle GPS Spoofing */
    printf("\n--- 2. Autonomous Vehicle GPS Spoofing ---\n");
    AVSecurity* av = av_security_create();
    for (int i = 0; i < 20; i++) {
        av_security_step(av, 0.1, 0.5, 0.0, 0.0);
        if (i == 10) {
            av_security_set_spoof(av, 30.0, 0.0);
            printf("[GPS Spoofing injected at step 10]\n");
        }
    }
    av_security_print_status(av);
    av_security_free(av);

    /* 3. ICS Replay Attack */
    printf("\n--- 3. ICS Replay Attack Defense with Watermarking ---\n");
    ICSReplayDefense* ics = ics_replay_defense_create();
    for (int i = 0; i < 30; i++) {
        ics_step_with_watermark(ics, 0.1, 0.2, 10.0);
        if (i == 15) {
            ics_replay_attack_start(ics);
            printf("[Replay attack started at step 15]\n");
        }
    }
    ics_replay_print_status(ics);
    ics_replay_defense_free(ics);

    /* 4. Water SCADA */
    printf("\n--- 4. Water Distribution SCADA Security ---\n");
    WaterSCADASecurity* ws = water_scada_create();
    for (int i = 0; i < 50; i++) {
        water_scada_step(ws, 0.1);
        if (i == 25) {
            water_scada_attack_start(ws, -0.8);
            printf("[Chlorine sensor attack started at step 25]\n");
        }
    }
    water_scada_print_status(ws);
    water_scada_free(ws);

    printf("\n=== All application demos complete ===\n");
}
