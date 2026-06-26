/**
 * @file    qc_dynamic_quantizer.c
 * @brief   Dynamic (zooming) quantizer implementation (Liberzon 2003)
 *
 * Implements the zoom-in/zoom-out quantization strategy for achieving
 * asymptotic stability with finite data rate. The quantizer range
 * parameter mu(t) adapts based on the system state.
 *
 * Algorithm (Liberzon 2003, Section 2.3):
 *   - Initialize with mu(0) > |x(0)| / M
 *   - At each step, quantize: q(x) = mu * round(x/(mu * Delta))
 *   - If |x| > M*mu: zoom_out (mu <- rho_out * mu)
 *   - If |x| < epsilon*mu for dwell_time: zoom_in (mu <- rho_in * mu)
 *
 * This hybrid feedback strategy guarantees:
 *   1. State remains bounded for all time
 *   2. State converges to the origin asymptotically
 *   3. Only a finite number of zoom-out events occur
 *
 * Key reference:
 *   - Liberzon, D. (2003). Switching in Systems and Control. Birkhauser.
 *   - Brockett & Liberzon (2000). IEEE TAC 45(7): 1271-1285.
 */

#include "quantized_control.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * Dynamic Quantizer Initialization and Operations
 * ================================================================ */

void qc_dyn_quantizer_init(QCDynamicQuantizer *dq, int bits, double M) {
    if (!dq) return;
    memset(dq, 0, sizeof(QCDynamicQuantizer));
    dq->bits = (bits > 0 && bits <= QC_MAX_BITS) ? bits : 8;
    dq->M = (M > 0.0) ? M : 1.0;
    dq->mu = 1.0;
    dq->mu_min = 1e-6;
    dq->mu_max = 1e6;
    dq->rho_in = 0.5;
    dq->rho_out = 2.0;
    dq->switching_threshold = 0.1;
    dq->zoom_level = 0;
    dq->max_zoom_level = 100;
    dq->zoom_in_count = 0;
    dq->zoom_out_count = 0;
    dq->dwell_time = 5;
    dq->time_since_last_zoom = 0.0;
}

double qc_dyn_quantize(QCDynamicQuantizer *dq, double x) {
    if (!dq) return 0.0;

    /* The quantizer has range [-M*mu, M*mu] with 2^bits levels.
     * Step size: Delta_mu = M*mu / (2^(bits-1))
     * q(x) = Delta_mu * round(x / Delta_mu)
     * Saturation: clamp to [-M*mu, M*mu]
     */
    int n_levels = 1 << dq->bits;
    double range = dq->M * dq->mu;
    double step = 2.0 * range / (double)(n_levels - 1);

    /* Clamp to range */
    double x_clamped = x;
    if (x_clamped > range) x_clamped = range;
    if (x_clamped < -range) x_clamped = -range;

    double xq = step * round(x_clamped / step);
    return xq;
}

void qc_dyn_zoom_in(QCDynamicQuantizer *dq) {
    if (!dq) return;
    double new_mu = dq->mu * dq->rho_in;
    if (new_mu >= dq->mu_min) {
        dq->mu = new_mu;
        dq->zoom_level--;
        dq->zoom_in_count++;
        dq->time_since_last_zoom = 0.0;
    }
}

void qc_dyn_zoom_out(QCDynamicQuantizer *dq) {
    if (!dq) return;
    double new_mu = dq->mu * dq->rho_out;
    if (new_mu <= dq->mu_max) {
        dq->mu = new_mu;
        dq->zoom_level++;
        dq->zoom_out_count++;
        dq->time_since_last_zoom = 0.0;
    }
}

int qc_dyn_should_zoom(const QCDynamicQuantizer *dq, double x) {
    if (!dq) return 0;
    double range = dq->M * dq->mu;

    /* Zoom out if state exceeds range (saturation imminent) */
    if (fabs(x) > range * 0.95) return 1;  /* zoom out */

    /* Zoom in if state is well within range and dwell time satisfied */
    if (fabs(x) < dq->switching_threshold * range &&
        dq->time_since_last_zoom >= dq->dwell_time) return -1; /* zoom in */

    return 0; /* no zoom needed */
}

double qc_dyn_get_range(const QCDynamicQuantizer *dq) {
    if (!dq) return 1.0;
    return dq->M * dq->mu;
}

void qc_dyn_reset(QCDynamicQuantizer *dq) {
    if (!dq) return;
    dq->mu = 1.0;
    dq->zoom_level = 0;
    dq->zoom_in_count = 0;
    dq->zoom_out_count = 0;
    dq->time_since_last_zoom = 0.0;
}

double qc_dyn_get_resolution(const QCDynamicQuantizer *dq) {
    if (!dq) return 0.0;
    int n_levels = 1 << dq->bits;
    double range = dq->M * dq->mu;
    return 2.0 * range / (double)(n_levels - 1);
}

/* ================================================================
 * Zoom Verification (Lyapunov-based)
 * ================================================================ */

void qc_zoom_verify_init(QCZoomVerification *zv) {
    if (!zv) return;
    memset(zv, 0, sizeof(QCZoomVerification));
    zv->zoom_sequence_len = 0;
    zv->bounded = 1;
    zv->converges_to_zero = 0;
}

int qc_zoom_verify_step(QCZoomVerification *zv, double V_now, double V_prev,
                         double state_norm) {
    if (!zv) return -1;

    /* Record state norm */
    zv->zoom_sequence_len++;
    zv->zoom_directions = realloc(zv->zoom_directions, zv->zoom_sequence_len * sizeof(int));
    zv->lyapunov_values = realloc(zv->lyapunov_values, zv->zoom_sequence_len * sizeof(double));
    zv->state_norms = realloc(zv->state_norms, zv->zoom_sequence_len * sizeof(double));
    if (!zv->zoom_directions || !zv->lyapunov_values || !zv->state_norms) return -1;

    int idx = zv->zoom_sequence_len - 1;
    zv->state_norms[idx] = state_norm;
    zv->lyapunov_values[idx] = V_now;

    if (state_norm > zv->max_state_norm) zv->max_state_norm = state_norm;

    /* Determine zoom direction from Lyapunov change */
    if (V_now < V_prev) {
        zv->zoom_directions[idx] = -1; /* zoom in (decreasing) */
        zv->num_zoom_ins++;
    } else if (V_now > V_prev * 1.1) {
        zv->zoom_directions[idx] = 1; /* zoom out (increasing) */
        zv->num_zoom_outs++;
    } else {
        zv->zoom_directions[idx] = 0; /* no change */
    }

    return 0;
}

int qc_zoom_verify_conclusion(QCZoomVerification *zv) {
    if (!zv || zv->zoom_sequence_len == 0) return 0;

    /* Check bounded: all state norms finite */
    zv->bounded = (zv->max_state_norm < INFINITY) ? 1 : 0;

    /* Check convergence: last V value much smaller than initial */
    if (zv->zoom_sequence_len >= 2) {
        double V_init = zv->lyapunov_values[0];
        double V_final = zv->lyapunov_values[zv->zoom_sequence_len - 1];
        zv->converges_to_zero = (V_final < V_init * 1e-3) ? 1 : 0;
        zv->convergence_rate = (V_init > 0) ? V_final / V_init : 0.0;
    }

    /* Successful zoom strategy:
     * - bounded (no escape)
     * - converges to zero
     * - finite number of zoom-outs (switching does not accumulate)
     */
    return (zv->bounded && zv->converges_to_zero) ? 1 : 0;
}

void qc_zoom_verify_free(QCZoomVerification *zv) {
    if (!zv) return;
    free(zv->zoom_directions);
    free(zv->lyapunov_values);
    free(zv->state_norms);
    zv->zoom_directions = NULL;
    zv->lyapunov_values = NULL;
    zv->state_norms = NULL;
    zv->zoom_sequence_len = 0;
}

/* ================================================================
 * Dynamic Quantizer with Full State Feedback
 * ================================================================ */

/**
 * Simulate one step of the zoom-based quantized control loop.
 * This implements the full Liberzon zoom strategy.
 *
 * Given:
 *   - System x_{k+1} = A x_k + B u_k
 *   - Quantized control u_k = q_mu(-K x_k)
 *   - Dynamic quantizer with parameter mu
 *
 * Update rule:
 *   1. Compute ideal control: u_ideal = -K x
 *   2. Quantize: u_q = q_mu(u_ideal)
 *   3. Evolve state: x_next = A x + B u_q
 *   4. Update mu based on |x|
 *
 * @return 1 if zoom-out occurred, -1 if zoom-in, 0 otherwise
 */
int qc_dyn_quantizer_step(QCDynamicQuantizer *dq, const double *A,
                           const double *B, const double *K,
                           int nx, int nu, double *x,
                           double *u, double dt) {
    if (!dq || !A || !B || !K || !x || !u || nx <= 0 || nu <= 0) return -1;

    /* Compute unconstrained LQR control */
    for (int i = 0; i < nu; i++) {
        u[i] = 0.0;
        for (int j = 0; j < nx; j++) {
            u[i] -= K[i * nx + j] * x[j];
        }
    }

    /* Quantize each control channel */
    for (int i = 0; i < nu; i++) {
        u[i] = qc_dyn_quantize(dq, u[i]);
    }

    /* Euler integration */
    double *x_new = calloc(nx, sizeof(double));
    if (!x_new) return -1;
    for (int i = 0; i < nx; i++) {
        x_new[i] = x[i];
        for (int j = 0; j < nx; j++) {
            x_new[i] += dt * A[i * nx + j] * x[j];
        }
        for (int j = 0; j < nu; j++) {
            x_new[i] += dt * B[i * nu + j] * u[j];
        }
    }

    /* Compute state norm and decide zoom */
    double norm_x = 0.0;
    for (int i = 0; i < nx; i++) norm_x += x_new[i] * x_new[i];
    norm_x = sqrt(norm_x);

    int zoom = qc_dyn_should_zoom(dq, norm_x);
    if (zoom > 0) qc_dyn_zoom_out(dq);
    else if (zoom < 0) qc_dyn_zoom_in(dq);

    dq->time_since_last_zoom += dt;

    /* Update state */
    memcpy(x, x_new, nx * sizeof(double));
    free(x_new);

    return zoom;
}
