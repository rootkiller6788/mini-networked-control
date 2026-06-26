/**
 * @file    quantized_control.c
 * @brief   Core implementation of quantized control system operations
 *
 * Implements system initialization, configuration, validation,
 * simulation, and the fundamental quantized control loop with
 * both input and output quantization.
 *
 * Mathematical foundation:
 *   dx/dt = A x + B q_u(u),  y_q = q_y(C x)
 *   where q_u, q_y are quantizers on actuator and sensor channels.
 *
 * Key references:
 *   - Elia & Mitter (2001). IEEE TAC.
 *   - Liberzon (2003). Switching in Systems and Control.
 *   - Nair & Evans (2004). SIAM J. Control Optim.
 *
 * Applications (L7):
 *   - Automotive: CAN bus quantized control in Tesla/ISO 11898 systems
 *   - Embedded: microcontroller ADC/DAC quantization in IoT sensor networks
 *   - Industrial: PROFINET/Fieldbus with limited feedback data rates
 */

#include "quantized_control.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * System Initialization and Configuration
 * ================================================================ */

void qc_system_init(QCSystem *sys) {
    if (!sys) return;
    memset(sys, 0, sizeof(QCSystem));
    sys->state_dim = 0;
    sys->input_dim = 0;
    sys->output_dim = 0;
    sys->sampling_period = 0.01;
    sys->sector_delta = 0.0;
    sys->sector_upper = 1.0;
    sys->stability = QC_STABLE_MARGINAL;
    qc_quantizer_init(&sys->input_quantizer, QC_QTYPE_UNIFORM, 8);
    qc_quantizer_init(&sys->output_quantizer, QC_QTYPE_UNIFORM, 8);
    qc_encoder_init(&sys->input_encoder, QC_ENC_FIXED_LENGTH, 8);
    qc_decoder_init(&sys->input_decoder, QC_ENC_FIXED_LENGTH, 8);
    qc_encoder_init(&sys->output_encoder, QC_ENC_FIXED_LENGTH, 8);
    qc_decoder_init(&sys->output_decoder, QC_ENC_FIXED_LENGTH, 8);
    qc_data_rate_init(&sys->data_rate, 1000.0, 1);
}

int qc_system_configure(QCSystem *sys, int nx, int nu, int ny) {
    if (!sys) return -1;
    if (nx <= 0 || nx > QC_MAX_STATE_DIM) return -1;
    if (nu < 0 || nu > QC_MAX_INPUT_DIM) return -1;
    if (ny < 0 || ny > QC_MAX_OUTPUT_DIM) return -1;
    sys->state_dim = nx;
    sys->input_dim = nu;
    sys->output_dim = ny;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            sys->A[i * nx + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < nx * nu; i++) sys->B[i] = 0.0;
    for (int i = 0; i < ny * nx; i++) sys->C[i] = 0.0;
    for (int i = 0; i < ny * nu; i++) sys->D[i] = 0.0;
    return 0;
}

int qc_system_set_A(QCSystem *sys, const double *A) {
    if (!sys || !A || sys->state_dim <= 0) return -1;
    memcpy(sys->A, A, sys->state_dim * sys->state_dim * sizeof(double));
    return 0;
}

int qc_system_set_B(QCSystem *sys, const double *B) {
    if (!sys || !B || sys->state_dim <= 0 || sys->input_dim <= 0) return -1;
    memcpy(sys->B, B, sys->state_dim * sys->input_dim * sizeof(double));
    return 0;
}

int qc_system_set_C(QCSystem *sys, const double *C) {
    if (!sys || !C || sys->output_dim <= 0 || sys->state_dim <= 0) return -1;
    memcpy(sys->C, C, sys->output_dim * sys->state_dim * sizeof(double));
    return 0;
}

int qc_system_set_D(QCSystem *sys, const double *D) {
    if (!sys || !D) return -1;
    int p = sys->output_dim, m = sys->input_dim;
    if (p <= 0 || m <= 0) return 0;
    memcpy(sys->D, D, p * m * sizeof(double));
    return 0;
}

int qc_system_validate(const QCSystem *sys) {
    if (!sys) return -1;
    if (sys->state_dim <= 0 || sys->state_dim > QC_MAX_STATE_DIM) return -2;
    if (sys->input_dim > QC_MAX_INPUT_DIM) return -3;
    if (sys->output_dim > QC_MAX_OUTPUT_DIM) return -4;
    if (sys->sampling_period <= 0.0) return -5;
    for (int i = 0; i < sys->state_dim * sys->state_dim; i++) {
        if (isnan(sys->A[i])) return -6;
    }
    return 0;
}

void qc_system_print(const QCSystem *sys) {
    if (!sys) { printf("QCSystem: NULL\n"); return; }
    printf("========================================\n");
    printf(" Quantized Control System\n");
    printf(" States: %d | Inputs: %d | Outputs: %d\n",
           sys->state_dim, sys->input_dim, sys->output_dim);
    printf(" Sampling: %.6f s | Sector: [%.4f, %.4f]\n",
           sys->sampling_period, sys->sector_delta, sys->sector_upper);
    printf(" Input Q: %s %d-bit | Output Q: %s %d-bit\n",
           qc_quantizer_type_name(sys->input_quantizer.type), sys->input_quantizer.bits,
           qc_quantizer_type_name(sys->output_quantizer.type), sys->output_quantizer.bits);
    printf(" Data rate: %.2f bps (min: %.2f, margin: %.2f)\n",
           sys->data_rate.total_rate_bps, sys->data_rate.min_rate, sys->data_rate.rate_margin);
    printf(" Stability: %s\n", qc_stability_status_name(sys->stability));
    printf(" A matrix (%dx%d):\n", sys->state_dim, sys->state_dim);
    for (int i = 0; i < sys->state_dim; i++) {
        printf("  ");
        for (int j = 0; j < sys->state_dim; j++) {
            printf("%8.4f ", sys->A[i * sys->state_dim + j]);
        }
        printf("\n");
    }
    printf("========================================\n");
}

size_t qc_system_memory_estimate(const QCSystem *sys) {
    if (!sys) return 0;
    size_t s = sizeof(QCSystem);
    if (sys->input_quantizer.levels) s += sys->input_quantizer.levels_len * sizeof(double);
    if (sys->output_quantizer.levels) s += sys->output_quantizer.levels_len * sizeof(double);
    if (sys->input_encoder.buffer) s += sys->input_encoder.buffer_size;
    if (sys->output_encoder.buffer) s += sys->output_encoder.buffer_size;
    if (sys->data_rate.channel_bits) s += sys->data_rate.num_channels * sizeof(int);
    if (sys->data_rate.channel_rate) s += sys->data_rate.num_channels * sizeof(double);
    return s;
}

/* ================================================================
 * Simulation Result Management
 * ================================================================ */

void qc_sim_result_init(QCSimulationResult *res) {
    if (!res) return;
    memset(res, 0, sizeof(QCSimulationResult));
    res->final_status = QC_STABLE_MARGINAL;
}

void qc_sim_result_free(QCSimulationResult *res) {
    if (!res) return;
    free(res->state_trajectory);
    free(res->input_trajectory);
    free(res->output_trajectory);
    free(res->time_points);
    res->state_trajectory = NULL;
    res->input_trajectory = NULL;
    res->output_trajectory = NULL;
    res->time_points = NULL;
    res->trajectory_len = 0;
    res->trajectory_capacity = 0;
}

static int sim_ensure_capacity(QCSimulationResult *res, int steps) {
    if (res->trajectory_capacity >= steps) return 0;
    int new_cap = steps * 2;
    int chunk_s = new_cap * QC_MAX_STATE_DIM;
    int chunk_i = new_cap * QC_MAX_INPUT_DIM;
    int chunk_o = new_cap * QC_MAX_OUTPUT_DIM;
    double *ns = realloc(res->state_trajectory, chunk_s * sizeof(double));
    double *ni = realloc(res->input_trajectory, chunk_i * sizeof(double));
    double *no = realloc(res->output_trajectory, chunk_o * sizeof(double));
    double *nt = realloc(res->time_points, new_cap * sizeof(double));
    if (!ns || !nt) {
        free(ns); free(ni); free(no); free(nt); return -1;
    }
    res->state_trajectory = ns;
    res->input_trajectory = ni;
    res->output_trajectory = no;
    res->time_points = nt;
    res->trajectory_capacity = new_cap;
    return 0;
}

static void sim_store_step(QCSimulationResult *res, const double *x, int nx,
                            const double *u, int nu, const double *y, int ny,
                            double t, int step) {
    int off_s = step * QC_MAX_STATE_DIM;
    int off_i = step * QC_MAX_INPUT_DIM;
    int off_o = step * QC_MAX_OUTPUT_DIM;
    for (int i = 0; i < nx; i++) res->state_trajectory[off_s + i] = x[i];
    for (int i = 0; i < nu; i++) res->input_trajectory[off_i + i] = u[i];
    for (int i = 0; i < ny; i++) res->output_trajectory[off_o + i] = y[i];
    res->time_points[step] = t;
}

/* ================================================================
 * Closed-Loop Simulation with Quantized Channels
 * ================================================================ */

int qc_simulate_closed_loop(QCSystem *sys, const double *x0,
                             double t0, double tf, double dt,
                             double (*controller)(const double*, int, double, double*),
                             QCSimulationResult *res) {
    if (!sys || !x0 || !res || tf <= t0 || dt <= 0) return -1;
    int nx = sys->state_dim, nu = sys->input_dim, ny = sys->output_dim;
    if (nx <= 0) return -1;
    int steps = (int)((tf - t0) / dt) + 1;
    if (sim_ensure_capacity(res, steps) != 0) return -1;

    double *x = calloc(nx, sizeof(double));
    double *dx = calloc(nx, sizeof(double));
    double *u = calloc(nu > 0 ? nu : 1, sizeof(double));
    double *u_q = calloc(nu > 0 ? nu : 1, sizeof(double));
    double *y = calloc(ny > 0 ? ny : 1, sizeof(double));
    double *y_q = calloc(ny > 0 ? ny : 1, sizeof(double));
    if (!x || !dx || !u || !u_q || !y || !y_q) {
        free(x); free(dx); free(u); free(u_q); free(y); free(y_q);
        return -1;
    }
    memcpy(x, x0, nx * sizeof(double));

    double max_q_err = 0.0, total_bits = 0.0;
    for (int k = 0; k <= steps; k++) {
        double t = t0 + k * dt;
        /* Output quantization */
        for (int i = 0; i < ny; i++) {
            y[i] = 0.0;
            for (int j = 0; j < nx; j++) y[i] += sys->C[i * nx + j] * x[j];
        }
        qc_quantize_vector(&sys->output_quantizer, y, ny, y_q);

        /* Controller acts on quantized measurement */
        double u_ideal[QC_MAX_INPUT_DIM];
        for (int i = 0; i < nu; i++) u_ideal[i] = 0.0;
        if (controller) controller(y_q, ny, t, u_ideal);

        /* Input quantization */
        qc_quantize_vector(&sys->input_quantizer, u_ideal, nu, u_q);

        /* Euler integration */
        for (int i = 0; i < nx; i++) {
            dx[i] = 0.0;
            for (int j = 0; j < nx; j++) dx[i] += sys->A[i * nx + j] * x[j];
            for (int j = 0; j < nu; j++) dx[i] += sys->B[i * nu + j] * u_q[j];
        }
        for (int i = 0; i < nu; i++) {
            double e = fabs(u_q[i] - u_ideal[i]);
            if (e > max_q_err) max_q_err = e;
        }
        qc_encoder_encode_vector(&sys->input_encoder, &sys->input_quantizer, u_ideal, nu);
        total_bits += sys->input_encoder.bits_per_symbol * nu;
        sim_store_step(res, x, nx, u_q, nu, y_q, ny, t, k);
        for (int i = 0; i < nx; i++) x[i] += dt * dx[i];
    }

    double ferr = 0.0;
    for (int i = 0; i < nx; i++) ferr += x[i] * x[i];
    ferr = sqrt(ferr);
    res->converged = (ferr < QC_EPSILON);
    res->steps = steps;
    res->final_error = ferr;
    res->max_quantization_error = max_q_err;
    res->avg_bit_rate = total_bits / (tf - t0);
    res->trajectory_len = steps + 1;
    if (ferr < QC_EPSILON) res->final_status = QC_STABLE_ASYMPTOTIC;
    else if (ferr < 1.0) res->final_status = QC_STABLE_PRACTICAL;
    else res->final_status = QC_STABLE_MARGINAL;

    free(x); free(dx); free(u); free(u_q); free(y); free(y_q);
    return 0;
}

int qc_simulate_quantized_lqr(QCSystem *sys, const double *x0,
                               const double *K, double t0, double tf,
                               double dt, QCSimulationResult *res) {
    if (!sys || !x0 || !K || !res || tf <= t0 || dt <= 0) return -1;
    int nx = sys->state_dim, nu = sys->input_dim;
    if (nx <= 0 || nu <= 0) return -1;
    int steps = (int)((tf - t0) / dt) + 1;
    if (sim_ensure_capacity(res, steps) != 0) return -1;

    double *x = calloc(nx, sizeof(double));
    double *dx = calloc(nx, sizeof(double));
    double *u_ideal = calloc(nu, sizeof(double));
    double *u_q = calloc(nu, sizeof(double));
    double *y = calloc(sys->output_dim > 0 ? sys->output_dim : 1, sizeof(double));
    if (!x || !dx || !u_ideal || !u_q || !y) {
        free(x); free(dx); free(u_ideal); free(u_q); free(y); return -1;
    }
    memcpy(x, x0, nx * sizeof(double));

    double max_q_err = 0.0;
    for (int k = 0; k <= steps; k++) {
        double t = t0 + k * dt;
        for (int i = 0; i < nu; i++) {
            u_ideal[i] = 0.0;
            for (int j = 0; j < nx; j++) u_ideal[i] -= K[i * nx + j] * x[j];
        }
        qc_quantize_vector(&sys->input_quantizer, u_ideal, nu, u_q);
        for (int i = 0; i < nx; i++) {
            dx[i] = 0.0;
            for (int j = 0; j < nx; j++) dx[i] += sys->A[i * nx + j] * x[j];
            for (int j = 0; j < nu; j++) dx[i] += sys->B[i * nu + j] * u_q[j];
        }
        for (int i = 0; i < nu; i++) {
            double e = fabs(u_q[i] - u_ideal[i]);
            if (e > max_q_err) max_q_err = e;
        }
        int ny = sys->output_dim > 0 ? sys->output_dim : 1;
        if (sys->output_dim > 0) {
            for (int i = 0; i < ny; i++) {
                y[i] = 0.0;
                for (int j = 0; j < nx; j++) y[i] += sys->C[i * nx + j] * x[j];
            }
        } else {
            y[0] = 0.0;
            for (int j = 0; j < nx; j++) y[0] += x[j] * x[j];
            y[0] = sqrt(y[0]);
        }
        sim_store_step(res, x, nx, u_q, nu, y, ny, t, k);
        for (int i = 0; i < nx; i++) x[i] += dt * dx[i];
    }

    double ferr = 0.0;
    for (int i = 0; i < nx; i++) ferr += x[i] * x[i];
    ferr = sqrt(ferr);
    res->converged = (ferr < QC_EPSILON);
    res->steps = steps; res->final_error = ferr;
    res->max_quantization_error = max_q_err;
    res->trajectory_len = steps + 1;
    if (ferr < QC_EPSILON) res->final_status = QC_STABLE_ASYMPTOTIC;
    else if (ferr < 1.0) res->final_status = QC_STABLE_PRACTICAL;
    else res->final_status = QC_STABLE_MARGINAL;

    free(x); free(dx); free(u_ideal); free(u_q); free(y);
    return 0;
}

int qc_simulate_output_feedback(QCSystem *sys, const double *x0,
                                 const double *K, const double *L,
                                 double t0, double tf, double dt,
                                 QCSimulationResult *res) {
    if (!sys || !x0 || !K || !L || !res) return -1;
    int nx = sys->state_dim, nu = sys->input_dim, ny = sys->output_dim;
    if (nx <= 0 || nu <= 0 || ny <= 0) return -1;
    int steps = (int)((tf - t0) / dt) + 1;
    if (sim_ensure_capacity(res, steps) != 0) return -1;

    double *x = calloc(nx, sizeof(double));
    double *xh = calloc(nx, sizeof(double));
    double *dx = calloc(nx, sizeof(double));
    double *dxh = calloc(nx, sizeof(double));
    double *u_ideal = calloc(nu, sizeof(double));
    double *u_q = calloc(nu, sizeof(double));
    double *y = calloc(ny, sizeof(double));
    double *y_q = calloc(ny, sizeof(double));
    if (!x || !xh || !dx || !dxh || !u_ideal || !u_q || !y || !y_q) {
        free(x); free(xh); free(dx); free(dxh);
        free(u_ideal); free(u_q); free(y); free(y_q); return -1;
    }
    memcpy(x, x0, nx * sizeof(double));
    memset(xh, 0, nx * sizeof(double));

    double max_q_err = 0.0;
    for (int k = 0; k <= steps; k++) {
        double t = t0 + k * dt;
        for (int i = 0; i < ny; i++) {
            y[i] = 0.0;
            for (int j = 0; j < nx; j++) y[i] += sys->C[i * nx + j] * x[j];
        }
        qc_quantize_vector(&sys->output_quantizer, y, ny, y_q);

        /* Observer: xh_dot = A*xh + B*u_q + L*(y_q - C*xh) */
        for (int i = 0; i < nx; i++) {
            double yh_i = 0.0;
            for (int j = 0; j < nx; j++) yh_i += sys->C[i % ny * nx + j] * xh[j];
            dxh[i] = 0.0;
            for (int j = 0; j < nx; j++) dxh[i] += sys->A[i * nx + j] * xh[j];
            for (int j = 0; j < nu; j++) dxh[i] += sys->B[i * nu + j] * u_q[j];
            for (int j = 0; j < ny; j++) dxh[i] += L[i * ny + j] * (y_q[j] - yh_i);
        }

        for (int i = 0; i < nu; i++) {
            u_ideal[i] = 0.0;
            for (int j = 0; j < nx; j++) u_ideal[i] -= K[i * nx + j] * xh[j];
        }
        qc_quantize_vector(&sys->input_quantizer, u_ideal, nu, u_q);

        for (int i = 0; i < nx; i++) {
            dx[i] = 0.0;
            for (int j = 0; j < nx; j++) dx[i] += sys->A[i * nx + j] * x[j];
            for (int j = 0; j < nu; j++) dx[i] += sys->B[i * nu + j] * u_q[j];
        }

        double qe = 0.0;
        for (int i = 0; i < nu; i++) qe += (u_q[i] - u_ideal[i]) * (u_q[i] - u_ideal[i]);
        qe = sqrt(qe);
        if (qe > max_q_err) max_q_err = qe;

        sim_store_step(res, x, nx, u_q, nu, y_q, ny, t, k);
        for (int i = 0; i < nx; i++) { x[i] += dt * dx[i]; xh[i] += dt * dxh[i]; }
    }

    double ferr = 0.0;
    for (int i = 0; i < nx; i++) ferr += x[i] * x[i];
    ferr = sqrt(ferr);
    res->converged = (ferr < 1e-3); res->steps = steps;
    res->final_error = ferr;
    res->max_quantization_error = max_q_err;
    res->trajectory_len = steps + 1;
    if (ferr < QC_EPSILON) res->final_status = QC_STABLE_ASYMPTOTIC;
    else if (ferr < 1.0) res->final_status = QC_STABLE_PRACTICAL;
    else res->final_status = QC_STABLE_MARGINAL;

    free(x); free(xh); free(dx); free(dxh);
    free(u_ideal); free(u_q); free(y); free(y_q);
    return 0;
}

void qc_sim_result_print(const QCSimulationResult *res) {
    if (!res) { printf("QCSimulationResult: NULL\n"); return; }
    printf("========================================\n");
    printf(" Simulation Result\n");
    printf(" Converged: %s | Steps: %d\n",
           res->converged ? "YES" : "NO", res->steps);
    printf(" Final error: %.6e | Max Q error: %.6e\n",
           res->final_error, res->max_quantization_error);
    printf(" Avg bit rate: %.2f bps | Status: %s\n",
           res->avg_bit_rate, qc_stability_status_name(res->final_status));
    printf("========================================\n");
}

/* ================================================================
 * Matrix Operations and Stability Analysis
 * ================================================================ */

int qc_matrix_eigenvalues(const double *A, int n,
                           double *eig_real, double *eig_imag) {
    if (!A || !eig_real || !eig_imag || n <= 0) return -1;
    /* Simple power iteration for dominant eigenvalue, with deflation.
     * For full eigenvalue decomposition, use QR algorithm. */
    double *B = malloc(n * n * sizeof(double));
    if (!B) return -1;
    memcpy(B, A, n * n * sizeof(double));

    for (int i = 0; i < n; i++) {
        /* Power iteration on remaining (n-i)x(n-i) submatrix */
        int m = n - i;
        double *v = calloc(m, sizeof(double));
        double *Av = calloc(m, sizeof(double));
        if (!v || !Av) { free(v); free(Av); free(B); return -1; }
        for (int j = 0; j < m; j++) v[j] = 1.0 / sqrt((double)m);

        for (int iter = 0; iter < 100; iter++) {
            double norm = 0.0;
            for (int r = 0; r < m; r++) {
                Av[r] = 0.0;
                for (int c = 0; c < m; c++) {
                    Av[r] += B[(i+r) * n + (i+c)] * v[c];
                }
                norm += Av[r] * Av[r];
            }
            norm = sqrt(norm);
            if (norm < 1e-15) break;
            for (int r = 0; r < m; r++) v[r] = Av[r] / norm;
        }

        /* Rayleigh quotient */
        double lambda = 0.0;
        for (int r = 0; r < m; r++) {
            double Avr = 0.0;
            for (int c = 0; c < m; c++) Avr += B[(i+r)*n + (i+c)] * v[c];
            lambda += v[r] * Avr;
        }
        eig_real[i] = lambda;
        eig_imag[i] = 0.0;

        /* Deflate: remove this eigenvalue */
        for (int r = 0; r < m; r++) {
            for (int c = 0; c < m; c++) {
                B[(i+r)*n + (i+c)] -= lambda * v[r] * v[c];
            }
        }
        free(v); free(Av);
    }
    free(B);
    return 0;
}

int qc_is_schur_stable(const double *A, int n) {
    if (!A || n <= 0) return 0;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1; }
    qc_matrix_eigenvalues(A, n, er, ei);
    int stable = 1;
    for (int i = 0; i < n; i++) {
        if (sqrt(er[i]*er[i] + ei[i]*ei[i]) >= 1.0 - QC_EPSILON) { stable = 0; break; }
    }
    free(er); free(ei);
    return stable;
}

int qc_is_hurwitz_stable(const double *A, int n) {
    if (!A || n <= 0) return 0;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1; }
    qc_matrix_eigenvalues(A, n, er, ei);
    int stable = 1;
    for (int i = 0; i < n; i++) {
        if (er[i] >= -QC_EPSILON) { stable = 0; break; }
    }
    free(er); free(ei);
    return stable;
}

double qc_spectral_radius(const double *A, int n) {
    if (!A || n <= 0) return -1.0;
    double *er = calloc(n, sizeof(double));
    double *ei = calloc(n, sizeof(double));
    if (!er || !ei) { free(er); free(ei); return -1.0; }
    qc_matrix_eigenvalues(A, n, er, ei);
    double mr = 0.0;
    for (int i = 0; i < n; i++) {
        double m = sqrt(er[i]*er[i] + ei[i]*ei[i]);
        if (m > mr) mr = m;
    }
    free(er); free(ei);
    return mr;
}

void qc_matrix_multiply(const double *A, const double *B, double *C,
                         int m, int n, int p) {
    if (!A || !B || !C || m <= 0 || n <= 0 || p <= 0) return;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            C[i * p + j] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i * p + j] += A[i * n + k] * B[k * p + j];
            }
        }
    }
}

int qc_solve_lyapunov(const double *A, const double *Q, double *P, int n) {
    if (!A || !Q || !P || n <= 0) return -1;
    if (n > 16) return -1;
    int n2 = n * n;
    double *M = calloc(n2 * n2, sizeof(double));
    double *b = calloc(n2, sizeof(double));
    if (!M || !b) { free(M); free(b); return -1; }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    int row = i * n + k, col = j * n + l;
                    if (i == j) M[row * n2 + col] += A[k * n + l];
                    if (k == l) M[row * n2 + col] += A[i * n + j];
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            b[i * n + j] = -Q[i * n + j];
    /* Gaussian elimination with partial pivoting */
    for (int k = 0; k < n2; k++) {
        int pivot = k;
        double mv = fabs(M[k * n2 + k]);
        for (int i = k + 1; i < n2; i++) {
            if (fabs(M[i * n2 + k]) > mv) { mv = fabs(M[i * n2 + k]); pivot = i; }
        }
        if (mv < 1e-12) { free(M); free(b); return -1; }
        if (pivot != k) {
            for (int j = 0; j < n2; j++) {
                double t = M[k * n2 + j]; M[k * n2 + j] = M[pivot * n2 + j]; M[pivot * n2 + j] = t;
            }
            double t = b[k]; b[k] = b[pivot]; b[pivot] = t;
        }
        for (int i = k + 1; i < n2; i++) {
            double f = M[i * n2 + k] / M[k * n2 + k];
            for (int j = k; j < n2; j++) M[i * n2 + j] -= f * M[k * n2 + j];
            b[i] -= f * b[k];
        }
    }
    for (int i = n2 - 1; i >= 0; i--) {
        double s = b[i];
        for (int j = i + 1; j < n2; j++) s -= M[i * n2 + j] * b[j];
        b[i] = s / M[i * n2 + i];
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P[i * n + j] = b[i * n + j];
    free(M);
    return 0;
}

/* ================================================================
 * Control Utilities
 * ================================================================ */

int qc_control_with_saturation(const QCSystem *sys, const double *x,
                                double *u, const double *u_max) {
    if (!sys || !x || !u || !u_max) return -1;
    int nu = sys->input_dim, nx = sys->state_dim;
    if (nu <= 0) return -1;
    for (int i = 0; i < nu; i++) {
        double ui = 0.0;
        for (int j = 0; j < nx; j++) ui -= sys->B[j * nu + i] * x[j];
        if (ui > u_max[i]) u[i] = u_max[i];
        else if (ui < -u_max[i]) u[i] = -u_max[i];
        else u[i] = ui;
    }
    return 0;
}

double qc_ultimate_bound_estimate(const QCSystem *sys) {
    if (!sys || sys->state_dim <= 0) return INFINITY;
    double rho = qc_spectral_radius(sys->A, sys->state_dim);
    if (rho >= 1.0) return INFINITY;
    double step = sys->input_quantizer.step;
    if (step <= 0) step = 0.01;
    double normB = 0.0;
    int sz = sys->state_dim * sys->input_dim;
    for (int i = 0; i < sz; i++) normB += sys->B[i] * sys->B[i];
    normB = sqrt(normB);
    return step * normB / (1.0 - rho);
}

/* ================================================================
 * String Utilities
 * ================================================================ */

const char* qc_quantizer_type_name(QCQuantizerType t) {
    switch (t) {
        case QC_QTYPE_UNIFORM: return "Uniform";
        case QC_QTYPE_LOGARITHMIC: return "Logarithmic";
        case QC_QTYPE_DYNAMIC: return "Dynamic";
        case QC_QTYPE_VECTOR: return "Vector";
        case QC_QTYPE_RANDOM_DITHER: return "Random-Dither";
        case QC_QTYPE_FIXED_RATE: return "Fixed-Rate";
        case QC_QTYPE_VARIABLE_RATE: return "Variable-Rate";
        default: return "Unknown";
    }
}

const char* qc_overload_strategy_name(QCOverloadStrategy s) {
    switch (s) {
        case QC_OVERLOAD_SATURATE: return "Saturate";
        case QC_OVERLOAD_ZOOM_OUT: return "Zoom-Out";
        case QC_OVERLOAD_MODULO: return "Modulo";
        case QC_OVERLOAD_EXTEND: return "Extend";
        default: return "Unknown";
    }
}

const char* qc_stability_status_name(QCStabilityStatus s) {
    switch (s) {
        case QC_STABLE_ASYMPTOTIC: return "Asymptotically-Stable";
        case QC_STABLE_PRACTICAL: return "Practically-Stable";
        case QC_STABLE_MARGINAL: return "Marginally-Stable";
        case QC_STABLE_UNSTABLE: return "Unstable";
        case QC_STABLE_ULTIMATE_BOUND: return "Ultimately-Bounded";
        default: return "Unknown";
    }
}

const char* qc_encoding_scheme_name(QCEncodingScheme e) {
    switch (e) {
        case QC_ENC_FIXED_LENGTH: return "Fixed-Length";
        case QC_ENC_HUFFMAN: return "Huffman";
        case QC_ENC_DIFFERENTIAL: return "Differential";
        case QC_ENC_ADAPTIVE: return "Adaptive";
        case QC_ENC_ARITHMETIC: return "Arithmetic";
        default: return "Unknown";
    }
}

double qc_bit_error_probability(double snr_db) {
    double snr_lin = pow(10.0, snr_db / 10.0);
    return 0.5 * erfc(sqrt(snr_lin));
}
