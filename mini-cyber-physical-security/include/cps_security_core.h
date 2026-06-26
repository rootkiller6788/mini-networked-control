#ifndef CPS_SECURITY_CORE_H
#define CPS_SECURITY_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Cyber-Physical System (CPS) Security Core Types
 *
 * Based on foundational works:
 *   Mo & Sinopoli (2010) — "False Data Injection Attacks in Control Systems"
 *   Pasqualetti, Dorfler & Bullo (2013) — "Attack Detection and Identification
 *                                          in Cyber-Physical Systems"
 *   Fawzi, Tabuada & Diggavi (2014) — "Secure Estimation and Control for CPS"
 *   Cardenas, Amin & Sastry (2008) — "Secure Control: Towards Survivable CPS"
 *   Smith (2015) — "Covert Misappropriation of Networked Control Systems"
 *   Ding, Han, Xiang, Ge & Zhang (2018) — "A Survey on Security Control and
 *                                          Attack Detection for Industrial CPS"
 * ============================================================================ */

/* --- Attack Type Enumeration (L1: Definitions) --- */

typedef enum {
    CPS_ATTACK_NONE = 0,
    CPS_ATTACK_DOS = 1,
    CPS_ATTACK_FDI = 2,
    CPS_ATTACK_REPLAY = 3,
    CPS_ATTACK_BIAS = 4,
    CPS_ATTACK_COVERT = 5,
    CPS_ATTACK_SURGE = 6,
    CPS_ATTACK_ZERO_DYNAMICS = 7
} CPSAttackType;

typedef enum {
    CPS_TARGET_SENSOR = 0,
    CPS_TARGET_ACTUATOR = 1,
    CPS_TARGET_CONTROLLER = 2,
    CPS_TARGET_NETWORK = 3,
    CPS_TARGET_ESTIMATOR = 4,
    CPS_TARGET_REFERENCE = 5
} CPSAttackTarget;

typedef enum {
    CPS_SECURE_NORMAL = 0,
    CPS_SECURE_SUSPICIOUS = 1,
    CPS_SECURE_ATTACKED = 2,
    CPS_SECURE_DEGRADED = 3,
    CPS_SECURE_RECOVERING = 4,
    CPS_SECURE_COMPROMISED = 5
} CPSSecurityState;

typedef enum {
    CPS_DETECT_CHI2 = 0,
    CPS_DETECT_CUSUM = 1,
    CPS_DETECT_WATERMARK = 2,
    CPS_DETECT_MMD = 3,
    CPS_DETECT_KL_DIVERGENCE = 4,
    CPS_DETECT_INTERVAL = 5,
    CPS_DETECT_SET_MEMBERSHIP = 6
} CPSDetectionMethod;

typedef enum {
    CPS_RESILIENT_HOLD = 0,
    CPS_RESILIENT_FALLBACK = 1,
    CPS_RESILIENT_RECONFIG = 2,
    CPS_RESILIENT_GAME = 3,
    CPS_RESILIENT_MPC = 4,
    CPS_RESILIENT_REDUNDANT = 5
} CPSResilientStrategy;

/* --- Matrix and Vector Structures (L3: Mathematical Structures) --- */

typedef struct {
    double* data;
    int rows;
    int cols;
    int owns_data;
} CPSMatrix;

typedef struct {
    double* components;
    int dimension;
} CPSVector;

/* --- Linear Dynamic System Model --- */

typedef struct {
    CPSMatrix A;        /* State: x[k+1] = A x[k] + B u[k] */
    CPSMatrix B;        /* Input matrix */
    CPSMatrix C;        /* Output: y[k] = C x[k] */
    CPSMatrix D;        /* Feedthrough */
    int n_states;
    int n_inputs;
    int n_outputs;
    double Q_scale;     /* Process noise */
    double R_scale;     /* Measurement noise */
    int* observable_indices;
    int n_observable;
    int* controllable_indices;
    int n_controllable;
} CPSDynamicalSystem;

/* --- Attack Model --- */

typedef struct {
    CPSAttackType type;
    CPSAttackTarget target;
    double start_time;
    double end_time;
    double magnitude;
    double frequency;
    double stealthiness;
    int target_index;
    double* attack_signal;
    int signal_length;
    int signal_capacity;
    int is_active;
    double* Gamma_a;
    double* Gamma_y;
    int n_attack_inputs;
    int n_attack_outputs;
} CPSAttackModel;

/* --- Detector State --- */

typedef struct {
    CPSDetectionMethod method;
    double threshold;
    double statistic;
    double alarm_time;
    int alarm_count;
    double false_positive_rate;
    double detection_rate;
    double* residual_history;
    int history_length;
    int history_capacity;
    int alarm_active;
    double cusum_pos;
    double cusum_neg;
    double cusum_drift;
    double cusum_reset;
    int chi2_df;
    double chi2_pvalue;
    double watermark_energy;
    double watermark_expected;
} CPSDetector;

/* --- Resilient Controller --- */

typedef struct {
    CPSResilientStrategy strategy;
    double* safe_input;
    double* fallback_input;
    double* current_input;
    int input_dim;
    double* measurement_buffer;
    int buffer_length;
    int buffer_capacity;
    double defender_cost;
    double attacker_cost;
    double saddle_value;
    int prediction_horizon;
    double constraint_margin;
    double* tightened_lb;
    double* tightened_ub;
    int* active_sensors;
    int* active_actuators;
    int n_active_sensors;
    int n_active_actuators;
} CPSResilientController;

/* --- Full CPS Security Model --- */

typedef struct {
    CPSDynamicalSystem* plant;
    CPSAttackModel* attack;
    CPSDetector* detector;
    CPSResilientController* resilient;
    double* state;
    double* true_state;
    double* measurement;
    double* control_input;
    int n_states;
    int n_measurements;
    int n_inputs;
    double current_time;
    double time_step;
    int step_count;
    CPSSecurityState security_state;
    double time_under_attack;
    double time_in_degraded;
    double accumulated_cost;
    double detection_delay;
    double recovery_time;
    double* state_log;
    double* residual_log;
    int* alarm_log;
    int log_length;
    int log_capacity;
} CPSSecuritySystem;

/* --- Core API --- */

CPSSecuritySystem* cps_security_create(int n_states, int n_inputs, int n_outputs);
void cps_security_free(CPSSecuritySystem* sys);
void cps_set_state_matrix(CPSSecuritySystem* sys, const double* A_data, int rows, int cols);
void cps_set_input_matrix(CPSSecuritySystem* sys, const double* B_data, int rows, int cols);
void cps_set_output_matrix(CPSSecuritySystem* sys, const double* C_data, int rows, int cols);
void cps_set_noise_covariances(CPSSecuritySystem* sys, double Q_scale, double R_scale);
void cps_set_initial_state(CPSSecuritySystem* sys, const double* x0);
void cps_measure(CPSSecuritySystem* sys);
void cps_step(CPSSecuritySystem* sys, const double* u);
void cps_evolve_true(CPSSecuritySystem* sys);
void cps_attack_configure(CPSSecuritySystem* sys, CPSAttackType type, CPSAttackTarget target, double magnitude, double start_time);
void cps_attack_start(CPSSecuritySystem* sys);
void cps_attack_stop(CPSSecuritySystem* sys);
void cps_attack_set_target_index(CPSSecuritySystem* sys, int index);
void cps_attack_generate_signal(CPSSecuritySystem* sys);
bool cps_attack_is_active(CPSSecuritySystem* sys);
void cps_attack_set_stealthiness(CPSSecuritySystem* sys, double stealth);
CPSSecurityState cps_security_get_state(CPSSecuritySystem* sys);
const char* cps_attack_type_name(CPSAttackType type);
const char* cps_security_state_name(CPSSecurityState state);
void cps_matrix_multiply(double* C, const double* A, const double* B, int m, int n, int p);
void cps_matrix_transpose(double* AT, const double* A, int rows, int cols);
double cps_matrix_det_2x2(double a, double b, double c, double d);
double cps_matrix_det_3x3(const double* A);
int cps_matrix_rank(const double* A, int rows, int cols, double tol);
void cps_matrix_inv_2x2(double* inv, double a, double b, double c, double d);
void cps_matrix_inv_3x3(double* inv, const double* A);
double cps_vector_norm(const double* v, int n);
double cps_vector_dot(const double* a, const double* b, int n);
void cps_gaussian_sample(double* out, int n, double sigma);
bool cps_is_observable(const CPSDynamicalSystem* sys);
bool cps_is_controllable(const CPSDynamicalSystem* sys);
double cps_observability_gramian(const CPSDynamicalSystem* sys, double* gramian, int horizon);
double cps_controllability_gramian(const CPSDynamicalSystem* sys, double* gramian, int horizon);

#endif /* CPS_SECURITY_CORE_H */
