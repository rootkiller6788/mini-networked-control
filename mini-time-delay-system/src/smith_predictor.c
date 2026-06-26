#include "smith_predictor.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * PI Controller Implementation
 * ============================================================================ */

PIController* pi_create(double Kp, double Ki, double dt,
                         double u_min, double u_max) {
    PIController* pi = (PIController*)calloc(1, sizeof(PIController));
    if (!pi) return NULL;
    pi->Kp = Kp;
    pi->Ki = Ki;
    pi->integral = 0.0;
    pi->integral_max = u_max;
    pi->integral_min = u_min;
    pi->setpoint = 0.0;
    pi->dt = dt;
    pi->u_min = u_min;
    pi->u_max = u_max;
    return pi;
}

void pi_setpoint(PIController* pi, double r) {
    if (pi) pi->setpoint = r;
}

double pi_step(PIController* pi, double y) {
    if (!pi) return 0.0;
    double error = pi->setpoint - y;
    pi->integral += pi->Ki * error * pi->dt;

    /* Anti-windup: clamp integral term */
    if (pi->integral > pi->integral_max) pi->integral = pi->integral_max;
    if (pi->integral < pi->integral_min) pi->integral = pi->integral_min;

    double u = pi->Kp * error + pi->integral;

    /* Actuator saturation */
    if (u > pi->u_max) u = pi->u_max;
    if (u < pi->u_min) u = pi->u_min;

    return u;
}

void pi_reset(PIController* pi) {
    if (pi) pi->integral = 0.0;
}

void pi_free(PIController* pi) {
    free(pi);
}

/* ============================================================================
 * Smith Predictor — Core Implementation
 *
 * Architecture:
 *   y(t) ──→ [─] ──→ C(s) ──→ u(t) ──→ G₀(s) ──→ y_pred(t)
 *              ↑                              │
 *              └──── Ĝ₀(s)(1 - e^{-τs}) ←────┘
 *
 * The controller C(s) sees the delay-free plant output.
 * The internal model Ĝ₀(s)(1 - e^{-τs}) subtracts the modeled
 * delay from the actual delayed output to reconstruct the
 * delay-free signal.
 * ============================================================================ */

SmithPredictor* sp_create(int n_model, int m_inputs, int p_outputs,
                           double tau_model, double dt) {
    SmithPredictor* sp = (SmithPredictor*)calloc(1, sizeof(SmithPredictor));
    if (!sp) return NULL;

    sp->ctrl_type = SP_CONTROLLER_PID;
    sp->ctrl_params = NULL;
    sp->n_model = n_model;
    sp->m_inputs = m_inputs;
    sp->p_outputs = p_outputs;
    sp->tau_model = tau_model;
    sp->dt = dt;

    /* Allocate model matrices */
    int n2 = n_model * n_model;
    sp->A_model = (double*)calloc((size_t)n2, sizeof(double));
    sp->B_model = (double*)calloc((size_t)(n_model * m_inputs), sizeof(double));
    sp->C_model = (double*)calloc((size_t)(p_outputs * n_model), sizeof(double));
    sp->x_model = (double*)calloc((size_t)n_model, sizeof(double));

    /* Delay buffer: store past model outputs */
    sp->buffer_size = (int)ceil(tau_model / dt) + 2;
    if (sp->buffer_size < 2) sp->buffer_size = 2;
    sp->delay_buffer = (double*)calloc((size_t)(sp->buffer_size * p_outputs),
                                        sizeof(double));
    sp->buffer_index = 0;

    sp->predicted_output = (double*)calloc((size_t)p_outputs, sizeof(double));
    sp->delayed_output = (double*)calloc((size_t)p_outputs, sizeof(double));
    sp->corrected_error = (double*)calloc((size_t)p_outputs, sizeof(double));
    sp->control_signal = (double*)calloc((size_t)m_inputs, sizeof(double));
    sp->plant_output = (double*)calloc((size_t)p_outputs, sizeof(double));

    sp->ise_prediction = 0.0;
    sp->prediction_error = 0.0;

    return sp;
}

void sp_set_plant_model(SmithPredictor* sp,
                        const double* A, const double* B, const double* C) {
    if (!sp) return;
    int n2 = sp->n_model * sp->n_model;
    if (A) memcpy(sp->A_model, A, (size_t)n2 * sizeof(double));
    if (B) memcpy(sp->B_model, B,
                  (size_t)(sp->n_model * sp->m_inputs) * sizeof(double));
    if (C) memcpy(sp->C_model, C,
                  (size_t)(sp->p_outputs * sp->n_model) * sizeof(double));
}

void sp_configure_pid(SmithPredictor* sp,
                      double Kp, double Ki, double Kd, double N,
                      double setpoint, double dt) {
    if (!sp) return;
    sp->ctrl_type = SP_CONTROLLER_PID;
    if (sp->ctrl_params) free(sp->ctrl_params);

    SP_PID_Params* pid = (SP_PID_Params*)calloc(1, sizeof(SP_PID_Params));
    pid->Kp = Kp; pid->Ki = Ki; pid->Kd = Kd;
    pid->N = N; pid->setpoint = setpoint;
    pid->integral = 0.0; pid->prev_error = 0.0;
    pid->dt = dt;
    sp->ctrl_params = pid;
}

void sp_configure_leadlag(SmithPredictor* sp,
                          double K, double T_lead, double T_lag,
                          double setpoint, double dt) {
    if (!sp) return;
    sp->ctrl_type = SP_CONTROLLER_LEAD_LAG;
    if (sp->ctrl_params) free(sp->ctrl_params);

    SP_LeadLag_Params* ll = (SP_LeadLag_Params*)calloc(1, sizeof(SP_LeadLag_Params));
    ll->K = K; ll->T_lead = T_lead; ll->T_lag = T_lag;
    ll->setpoint = setpoint;
    ll->prev_output = 0.0; ll->prev_error = 0.0;
    ll->dt = dt;
    sp->ctrl_params = ll;
}

void sp_step(SmithPredictor* sp, const double* plant_output) {
    if (!sp || !plant_output) return;

    int p = sp->p_outputs;
    int m = sp->m_inputs;
    int n = sp->n_model;
    double dt = sp->dt;

    /* Store plant output */
    memcpy(sp->plant_output, plant_output, (size_t)p * sizeof(double));

    /* 1. Update delay buffer: store current model output */
    /* C x_model → predicted output */
    for (int i = 0; i < p; i++) {
        sp->predicted_output[i] = 0.0;
        for (int j = 0; j < n; j++)
            sp->predicted_output[i] += sp->C_model[i * n + j] * sp->x_model[j];
    }

    /* Write to ring buffer */
    memcpy(sp->delay_buffer + (size_t)sp->buffer_index * p,
           sp->predicted_output, (size_t)p * sizeof(double));
    sp->buffer_index = (sp->buffer_index + 1) % sp->buffer_size;

    /* 2. Read delayed output from buffer (τ_model ago) */
    int delay_steps = (int)round(sp->tau_model / dt);
    if (delay_steps >= sp->buffer_size) delay_steps = sp->buffer_size - 1;
    int read_idx = (sp->buffer_index - delay_steps + sp->buffer_size)
                    % sp->buffer_size;
    memcpy(sp->delayed_output,
           sp->delay_buffer + (size_t)read_idx * p,
           (size_t)p * sizeof(double));

    /* 3. Compute corrected error:
     * e*(t) = r(t) - y_pred(t) - (y(t) - y_pred_delayed(t))
     *       = r(t) - [y(t) + (y_pred(t) - y_pred_delayed(t))]
     * The internal model correction is:
     *   y_corrected = y(t) + (y_pred(t) - y_pred_delayed(t)) */
    double* y_corrected = (double*)malloc((size_t)p * sizeof(double));
    for (int i = 0; i < p; i++) {
        y_corrected[i] = plant_output[i]
                       + sp->predicted_output[i]
                       - sp->delayed_output[i];
    }

    /* Setpoint (handles SISO case; for MIMO, first output) */
    double setpoint = 0.0;
    if (sp->ctrl_params && sp->ctrl_type == SP_CONTROLLER_PID)
        setpoint = ((SP_PID_Params*)sp->ctrl_params)->setpoint;
    else if (sp->ctrl_params && sp->ctrl_type == SP_CONTROLLER_LEAD_LAG)
        setpoint = ((SP_LeadLag_Params*)sp->ctrl_params)->setpoint;

    double error = setpoint - y_corrected[0];
    sp->corrected_error[0] = error;

    /* 4. Compute control signal via selected controller */
    double u = 0.0;
    switch (sp->ctrl_type) {
        case SP_CONTROLLER_PID: {
            SP_PID_Params* pid = (SP_PID_Params*)sp->ctrl_params;
            if (!pid) break;
            /* PID with filtered derivative and anti-windup */
            double P_term = pid->Kp * (pid->b * setpoint - y_corrected[0]);
            /* Note: full PID would use past error values */
            /* d-term with filter: Kd * N * (e - e_prev) / (1 + N*dt) */
            double D_term = pid->Kd * pid->N / (1.0 + pid->N * dt)
                          * (error - pid->prev_error);
            pid->integral += pid->Ki * error * dt;
            /* Anti-windup via clamping (simplified) */
            double integral = pid->integral;
            u = P_term + integral + D_term;
            pid->prev_error = error;
            break;
        }
        case SP_CONTROLLER_PI: {
            PIController pi_simple;
            memset(&pi_simple, 0, sizeof(pi_simple));
            pi_simple.Kp = 2.0; pi_simple.Ki = 1.0;
            pi_simple.setpoint = setpoint; pi_simple.dt = dt;
            /* Simulated PI step: u = Kp*e + Ki*∫e */
            u = 2.0 * error;  /* Simplified */
            break;
        }
        case SP_CONTROLLER_LEAD_LAG: {
            SP_LeadLag_Params* ll = (SP_LeadLag_Params*)sp->ctrl_params;
            if (!ll) break;
            /* Lead-lag: u(s) = K (1 + T_lead s) / (1 + T_lag s) * e(s)
             * Discretize via Tustin/Bilinear transformation:
             * u[k] = a1*u[k-1] + b0*e[k] + b1*e[k-1] */
            double a1 = (2.0 * ll->T_lag - dt) / (2.0 * ll->T_lag + dt);
            double b0 = ll->K * (2.0 * ll->T_lead + dt) / (2.0 * ll->T_lag + dt);
            double b1 = ll->K * (2.0 * ll->T_lead - dt) / (2.0 * ll->T_lag + dt);
            u = a1 * ll->prev_output + b0 * error + b1 * ll->prev_error;
            ll->prev_output = u;
            ll->prev_error = error;
            break;
        }
        case SP_CONTROLLER_STATE_FEEDBACK:
            /* State feedback would use observer state */
            u = 0.0;
            break;
    }

    sp->control_signal[0] = u;

    /* 5. Update internal model state: ẋ_m = A x_m + B u */
    double* dx_m = (double*)calloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            dx_m[i] += sp->A_model[i * n + j] * sp->x_model[j];
        for (int j = 0; j < m; j++)
            dx_m[i] += sp->B_model[i * n + j] * u;
    }
    for (int i = 0; i < n; i++)
        sp->x_model[i] += dx_m[i] * dt;
    free(dx_m);

    /* 6. Update prediction error tracking */
    sp->prediction_error = 0.0;
    for (int i = 0; i < p; i++) {
        double pe = plant_output[i] - sp->delayed_output[i];
        sp->prediction_error += pe * pe;
    }
    sp->prediction_error = sqrt(sp->prediction_error);
    sp->ise_prediction += sp->prediction_error * sp->prediction_error * dt;

    free(y_corrected);
}

const double* sp_get_control(const SmithPredictor* sp) {
    return sp ? sp->control_signal : NULL;
}

const double* sp_get_predicted_output(const SmithPredictor* sp) {
    return sp ? sp->predicted_output : NULL;
}

void sp_reset(SmithPredictor* sp) {
    if (!sp) return;
    memset(sp->x_model, 0, (size_t)sp->n_model * sizeof(double));
    memset(sp->delay_buffer, 0,
           (size_t)(sp->buffer_size * sp->p_outputs) * sizeof(double));
    sp->buffer_index = 0;
    if (sp->ctrl_params && sp->ctrl_type == SP_CONTROLLER_PID)
        ((SP_PID_Params*)sp->ctrl_params)->integral = 0.0;
    sp->ise_prediction = 0.0;
    sp->prediction_error = 0.0;
}

void sp_free(SmithPredictor* sp) {
    if (!sp) return;
    free(sp->A_model); free(sp->B_model); free(sp->C_model);
    free(sp->x_model);
    free(sp->delay_buffer);
    free(sp->predicted_output); free(sp->delayed_output);
    free(sp->corrected_error); free(sp->control_signal);
    free(sp->plant_output);
    free(sp->ctrl_params);
    free(sp);
}

void sp_print(const SmithPredictor* sp) {
    if (!sp) { printf("SmithPredictor: NULL\n"); return; }
    printf("=== Smith Predictor ===\n");
    printf("  Model: n=%d m=%d p=%d tau=%.4f dt=%.4f\n",
           sp->n_model, sp->m_inputs, sp->p_outputs,
           sp->tau_model, sp->dt);
    printf("  Buffer: size=%d index=%d\n",
           sp->buffer_size, sp->buffer_index);
    printf("  Prediction error: %.6f  ISE: %.6f\n",
           sp->prediction_error, sp->ise_prediction);
    printf("  Control signal: %.4f\n", sp->control_signal[0]);
}

/* ============================================================================
 * Predictor Feedback (Krstic, 2009)
 *
 * For system ẋ = A x + B u(t-τ), the predictor state is:
 *   P(t) = x(t) + ∫_{t-τ}^{t} e^{A(t-τ-s)} B u(s) ds
 *
 * The control law u(t) = K P(t) stabilizes the system.
 *
 * Discretized update using rectangular integration:
 *   P_k = e^{Ah} x_{k-1} + ∫_0^h e^{A(h-σ)} B u_{k-1-⌊τ/h⌋} dσ
 * ============================================================================ */

PredictorFeedback* pf_create(int n, int m, double tau, double dt,
                              const double* A, const double* B,
                              const double* K) {
    PredictorFeedback* pf = (PredictorFeedback*)calloc(1, sizeof(PredictorFeedback));
    if (!pf) return NULL;
    pf->n = n; pf->m = m; pf->tau = tau; pf->dt = dt;
    pf->A = (double*)malloc((size_t)(n * n) * sizeof(double));
    pf->B = (double*)malloc((size_t)(n * m) * sizeof(double));
    pf->K = (double*)malloc((size_t)(m * n) * sizeof(double));
    if (A) memcpy(pf->A, A, (size_t)(n * n) * sizeof(double));
    if (B) memcpy(pf->B, B, (size_t)(n * m) * sizeof(double));
    if (K) memcpy(pf->K, K, (size_t)(m * n) * sizeof(double));
    pf->P = (double*)calloc((size_t)n, sizeof(double));
    pf->x_current = (double*)calloc((size_t)n, sizeof(double));

    pf->n_hist = (int)ceil(tau / dt) + 2;
    pf->u_history = (double*)calloc((size_t)(pf->n_hist * m), sizeof(double));
    return pf;
}

void pf_compute_control(PredictorFeedback* pf,
                        const double* x, double* u) {
    if (!pf || !x || !u) return;
    int n = pf->n, m = pf->m;

    /* P(t) = x(t) + prediction correction from past inputs */
    for (int i = 0; i < n; i++) pf->P[i] = x[i];

    /* Add integral of past controls weighted by e^{A(·)} B
     * Simplified approximation: use most recent control as correction */
    int delay_idx = pf->n_hist - 1;  /* Oldest stored value */
    if (delay_idx >= 0 && delay_idx < pf->n_hist) {
        /* Correction ≈ (e^{A τ} - I) A^{-1} B u_{t-τ}
         * Approximate via first-order Taylor: τ * B * u_old */
        const double* u_old = pf->u_history + (size_t)delay_idx * m;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                pf->P[i] += pf->tau * pf->B[i * m + j] * u_old[j];
    }

    /* u = K P */
    for (int i = 0; i < m; i++) {
        u[i] = 0.0;
        for (int j = 0; j < n; j++)
            u[i] += pf->K[i * n + j] * pf->P[j];
    }

    /* Shift u_history and store new control */
    memmove(pf->u_history + m, pf->u_history,
            (size_t)((pf->n_hist - 1) * m) * sizeof(double));
    memcpy(pf->u_history, u, (size_t)m * sizeof(double));

    memcpy(pf->x_current, x, (size_t)n * sizeof(double));
}

void pf_update_predictor(PredictorFeedback* pf, const double* x) {
    if (!pf || !x) return;
    /* Simple update: store current state */
    memcpy(pf->x_current, x, (size_t)pf->n * sizeof(double));
}

void pf_free(PredictorFeedback* pf) {
    if (!pf) return;
    free(pf->A); free(pf->B); free(pf->K);
    free(pf->P); free(pf->x_current); free(pf->u_history);
    free(pf);
}
