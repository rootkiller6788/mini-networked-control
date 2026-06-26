#ifndef SMITH_PREDICTOR_H
#define SMITH_PREDICTOR_H

#include "time_delay_system.h"

/* ============================================================================
 * Smith Predictor — Delay Compensation for Control Systems
 *
 * Reference:
 *   O. J. M. Smith, "A controller to overcome dead time," ISA Journal, 1959
 *     — Original Smith predictor paper
 *   Z. J. Palmor, "Time-delay compensation — Smith predictor and its
 *     modifications" in The Control Handbook (1996)
 *   K. J. Åström & T. Hägglund, "PID Controllers" (1995) — Ch. 6 on
 *     dead-time compensation
 *   J. E. Normey-Rico & E. F. Camacho, "Control of Dead-time Processes"
 *     (2007) — comprehensive treatment
 *
 * Core idea:
 *   For plant G(s) = G₀(s) e^{-τs}, use internal model G₀(s)(1 - e^{-τs})
 *   to predict the delay-free output. The controller acts on the predicted
 *   output, removing the delay from the characteristic equation.
 *
 * Level 5 — Algorithm: Smith Predictor Implementation
 * ============================================================================ */

/* ============================================================================
 * Smith Predictor Configuration
 * ============================================================================ */

/* --- Controller type for the primary (delay-free) loop --- */
typedef enum {
    SP_CONTROLLER_PID = 0,          /* Standard PID: Kp, Ki, Kd */
    SP_CONTROLLER_PI = 1,           /* PI only */
    SP_CONTROLLER_LEAD_LAG = 2,     /* Lead-lag compensator */
    SP_CONTROLLER_STATE_FEEDBACK = 3/* State feedback with observer */
} SPControllerType;

/* --- PID parameters --- */
typedef struct {
    double Kp;          /* Proportional gain */
    double Ki;          /* Integral gain */
    double Kd;          /* Derivative gain */
    double N;           /* Derivative filter coefficient (for proper D-term) */
    double b;           /* Setpoint weighting factor */
    double setpoint;    /* Reference signal */
    double integral;    /* Accumulated integral error (internal state) */
    double prev_error;  /* Previous error for derivative computation */
    double dt;          /* Time step */
} SP_PID_Params;

/* --- Lead-Lag compensator parameters --- */
typedef struct {
    double K;           /* DC gain */
    double T_lead;      /* Lead time constant */
    double T_lag;       /* Lag time constant */
    double setpoint;
    double prev_output;
    double prev_error;
    double dt;
} SP_LeadLag_Params;

/* --- Smith Predictor main structure --- */
typedef struct {
    /* Primary controller (acts on predicted output) */
    SPControllerType ctrl_type;
    void* ctrl_params;          /* SP_PID_Params* or SP_LeadLag_Params* */

    /* Plant model: Ĝ₀(s) — the delay-free part */
    /* Represented as state-space: ẋ_m = A_m x_m + B_m u */
    double* A_model;            /* n×n */
    double* B_model;            /* n×m */
    double* C_model;            /* p×n */
    double* x_model;            /* n×1 model state */
    int n_model;                /* Model state dimension */
    int m_inputs;
    int p_outputs;

    /* Delay model */
    double tau_model;           /* Modeled delay (seconds) */
    double* delay_buffer;       /* Ring buffer for delayed signals */
    int buffer_size;            /* Number of samples in delay buffer */
    int buffer_index;           /* Write index */
    double dt;                  /* Sample time */

    /* Internal signals */
    double* predicted_output;   /* ŷ(t) = C x_m(t)  (delay-free prediction) */
    double* delayed_output;     /* ŷ(t-τ)  (delayed prediction) */
    double* corrected_error;    /* e*(t) = r - ŷ(t) - (y(t) - ŷ(t-τ)) */
    double* control_signal;     /* u(t) */
    double* plant_output;       /* Latest plant measurement y(t) */

    /* Performance tracking */
    double ise_prediction;      /* ISE between prediction and actual */
    double prediction_error;    /* Current prediction error ||y - ŷ(t-τ)|| */

} SmithPredictor;

/* ============================================================================
 * Predictor-based Feedback Control (Krstic, 2009)
 * ============================================================================ */

/* Predictor feedback for linear systems with input delay:
 *   ẋ(t) = A x(t) + B u(t-τ)
 *   u(t) = K P(t)  where P(t) = x(t) + ∫_{t-τ}^{t} e^{A(t-τ-s)} B u(s) ds
 *
 * This is the continuous-time version, solved via PDE backstepping. */
typedef struct {
    double* A;            /* n×n */
    double* B;            /* n×m */
    double* K;            /* m×n gain matrix */
    double tau;           /* Input delay */
    double* P;            /* Predicted state P(t) */
    double* x_current;    /* Current state */
    double* u_history;    /* Control history buffer */
    int n_hist;           /* History length */
    int n;                /* State dimension */
    int m;                /* Input dimension */
    double dt;
} PredictorFeedback;

/* ============================================================================
 * Smith Predictor API
 * ============================================================================ */

/* Create Smith Predictor */
SmithPredictor* sp_create(int n_model, int m_inputs, int p_outputs,
                           double tau_model, double dt);

/* Set the plant model */
void sp_set_plant_model(SmithPredictor* sp,
                        const double* A, const double* B, const double* C);

/* Configure PID controller */
void sp_configure_pid(SmithPredictor* sp,
                      double Kp, double Ki, double Kd, double N,
                      double setpoint, double dt);

/* Configure Lead-Lag controller */
void sp_configure_leadlag(SmithPredictor* sp,
                          double K, double T_lead, double T_lag,
                          double setpoint, double dt);

/* Execute one time step of the Smith predictor:
 *   Input:  plant_output = y(t)  (actual plant output)
 *   Output: control_signal = u(t) (computed control)
 * Also updates internal model state and delayed prediction. */
void sp_step(SmithPredictor* sp, const double* plant_output);

/* Get the current control signal */
const double* sp_get_control(const SmithPredictor* sp);

/* Get the predicted (delay-free) output — useful for monitoring */
const double* sp_get_predicted_output(const SmithPredictor* sp);

/* Reset internal states */
void sp_reset(SmithPredictor* sp);

/* Free Smith predictor */
void sp_free(SmithPredictor* sp);

/* Print Smith predictor status */
void sp_print(const SmithPredictor* sp);

/* ============================================================================
 * PI Controller (standalone, for comparison with Smith predictor)
 * ============================================================================ */

/* Standard PI with anti-windup */
typedef struct {
    double Kp, Ki;
    double integral;
    double integral_max;    /* Anti-windup limit */
    double integral_min;
    double setpoint;
    double dt;
    double u_min;           /* Control saturation limits */
    double u_max;
} PIController;

PIController* pi_create(double Kp, double Ki, double dt,
                         double u_min, double u_max);
void pi_setpoint(PIController* pi, double r);
double pi_step(PIController* pi, double y);
void pi_reset(PIController* pi);
void pi_free(PIController* pi);

/* ============================================================================
 * Predictor Feedback (Krstic) API
 * ============================================================================ */

/* Create predictor feedback controller */
PredictorFeedback* pf_create(int n, int m, double tau, double dt,
                              const double* A, const double* B,
                              const double* K);

/* Compute control u(t) = K P(t) given current state x(t) */
void pf_compute_control(PredictorFeedback* pf,
                        const double* x, double* u);

/* Update predictor state (call each time step) */
void pf_update_predictor(PredictorFeedback* pf, const double* x);

/* Free predictor feedback */
void pf_free(PredictorFeedback* pf);

#endif /* SMITH_PREDICTOR_H */
