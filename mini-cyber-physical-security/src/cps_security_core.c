#include "cps_security_core.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CPS_INF 1e308
#define PI 3.14159265358979323846
#define CPS_DEFAULT_LOG_CAPACITY 1000
#define CPS_EPS 1e-12

/* ============================================================================
 * Internal helpers
 * ============================================================================ */

static double urand(void) {
    return (double)rand() / (double)RAND_MAX;
}

static double box_muller(double mu, double sigma) {
    double u1 = urand();
    double u2 = urand();
    if (u1 < 1e-15) u1 = 1e-15;
    return mu + sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

static double* vec_alloc(int n) {
    return (double*)calloc((size_t)n, sizeof(double));
}

static void mat_data_copy(double* dst, const double* src, int rows, int cols) {
    memcpy(dst, src, (size_t)(rows * cols) * sizeof(double));
}

/* ============================================================================
 * QR Decomposition via Householder Reflections (L3: Mathematical Structures)
 *
 * Decomposes A (m x n) into Q (m x m orthogonal) and R (m x n upper triangular)
 * such that A = Q * R. Used for numerical rank computation.
 * Complexity: O(m * n^2). Reference: Golub & Van Loan (2013)
 * ============================================================================ */

static void householder_qr(const double* A, int m, int n,
                            double* Q, double* R) {
    double* A_copy = (double*)malloc((size_t)(m * n) * sizeof(double));
    memcpy(A_copy, A, (size_t)(m * n) * sizeof(double));
    double* v = (double*)malloc((size_t)m * sizeof(double));

    for (int i = 0; i < m * m; i++) Q[i] = 0.0;
    for (int i = 0; i < m; i++) Q[i * m + i] = 1.0;
    for (int i = 0; i < m * n; i++) R[i] = A_copy[i];

    int k_max = (m < n) ? m : n;
    for (int k = 0; k < k_max; k++) {
        double norm_x = 0.0;
        for (int i = k; i < m; i++)
            norm_x += R[i * n + k] * R[i * n + k];
        norm_x = sqrt(norm_x);
        if (norm_x < CPS_EPS) continue;

        double alpha = (R[k * n + k] > 0) ? -norm_x : norm_x;
        for (int i = k; i < m; i++) v[i] = R[i * n + k];
        v[k] -= alpha;

        double v_norm_sq = 0.0;
        for (int i = k; i < m; i++) v_norm_sq += v[i] * v[i];
        if (v_norm_sq < CPS_EPS) continue;
        double beta = 2.0 / v_norm_sq;

        for (int j = k; j < n; j++) {
            double dot = 0.0;
            for (int i = k; i < m; i++) dot += v[i] * R[i * n + j];
            double tau = beta * dot;
            for (int i = k; i < m; i++) R[i * n + j] -= tau * v[i];
        }

        for (int j = 0; j < m; j++) {
            double dot = 0.0;
            for (int i = k; i < m; i++) dot += v[i] * Q[j * m + i];
            double tau = beta * dot;
            for (int i = k; i < m; i++) Q[j * m + i] -= tau * v[i];
        }

        R[k * n + k] = alpha;
        for (int i = k + 1; i < m; i++) R[i * n + k] = 0.0;
    }

    free(A_copy);
    free(v);
}

/* ============================================================================
 * Matrix Operations (L3: Mathematical Structures)
 * ============================================================================ */

void cps_matrix_multiply(double* C, const double* A, const double* B,
                          int m, int n, int p) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A[i * n + k] * B[k * p + j];
            C[i * p + j] = sum;
        }
    }
}

void cps_matrix_transpose(double* AT, const double* A, int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            AT[j * rows + i] = A[i * cols + j];
}

double cps_matrix_det_2x2(double a, double b, double c, double d) {
    return a * d - b * c;
}

double cps_matrix_det_3x3(const double* A) {
    return A[0] * (A[4] * A[8] - A[5] * A[7])
         - A[1] * (A[3] * A[8] - A[5] * A[6])
         + A[2] * (A[3] * A[7] - A[4] * A[6]);
}

void cps_matrix_inv_2x2(double* inv, double a, double b, double c, double d) {
    double det = a * d - b * c;
    if (fabs(det) < CPS_EPS) {
        inv[0] = inv[1] = inv[2] = inv[3] = 0.0;
        return;
    }
    double det_inv = 1.0 / det;
    inv[0] =  d * det_inv;  inv[1] = -b * det_inv;
    inv[2] = -c * det_inv;  inv[3] =  a * det_inv;
}

void cps_matrix_inv_3x3(double* inv, const double* A) {
    double det = cps_matrix_det_3x3(A);
    if (fabs(det) < CPS_EPS) {
        for (int i = 0; i < 9; i++) inv[i] = 0.0;
        return;
    }
    double d = 1.0 / det;
    inv[0] = (A[4]*A[8] - A[5]*A[7]) * d;
    inv[1] = (A[2]*A[7] - A[1]*A[8]) * d;
    inv[2] = (A[1]*A[5] - A[2]*A[4]) * d;
    inv[3] = (A[5]*A[6] - A[3]*A[8]) * d;
    inv[4] = (A[0]*A[8] - A[2]*A[6]) * d;
    inv[5] = (A[2]*A[3] - A[0]*A[5]) * d;
    inv[6] = (A[3]*A[7] - A[4]*A[6]) * d;
    inv[7] = (A[1]*A[6] - A[0]*A[7]) * d;
    inv[8] = (A[0]*A[4] - A[1]*A[3]) * d;
}

int cps_matrix_rank(const double* A, int rows, int cols, double tol) {
    if (rows <= 0 || cols <= 0) return 0;
    int m = rows, n = cols;
    double* Q = (double*)malloc((size_t)(m * m) * sizeof(double));
    double* R = (double*)calloc((size_t)(m * n), sizeof(double));
    householder_qr(A, m, n, Q, R);
    int rank = 0;
    int k_max = (m < n) ? m : n;
    for (int i = 0; i < k_max; i++) {
        if (fabs(R[i * n + i]) > tol) rank++;
    }
    free(Q); free(R);
    return rank;
}

/* ============================================================================
 * Vector Operations
 * ============================================================================ */

double cps_vector_norm(const double* v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

double cps_vector_dot(const double* a, const double* b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

void cps_gaussian_sample(double* out, int n, double sigma) {
    for (int i = 0; i < n; i++)
        out[i] = box_muller(0.0, sigma);
}

/* ============================================================================
 * Observability Check (L4: Fundamental Laws)
 *
 * Theorem (Kalman, 1960): System (A,C) is observable iff
 *   rank([C; CA; CA^2; ...; CA^{n-1}]) = n
 *
 * CPS Security significance: Observability determines whether an
 * attacker can modify sensor measurements without affecting the
 * residual. Unobservable attacks are inherently undetectable
 * (zero-dynamics attacks exploit this).
 * ============================================================================ */

bool cps_is_observable(const CPSDynamicalSystem* sys) {
    int n = sys->n_states;
    int p = sys->n_outputs;
    if (n <= 0 || p <= 0) return false;
    int m_rows = n * p;
    double* O = (double*)calloc((size_t)(m_rows * n), sizeof(double));
    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp = (double*)malloc((size_t)(n * n) * sizeof(double));

    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            O[i * n + j] = sys->C.data[i * n + j];

    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;

    for (int k = 1; k < n; k++) {
        cps_matrix_multiply(temp, A_pow, sys->A.data, n, n, n);
        memcpy(A_pow, temp, (size_t)(n * n) * sizeof(double));
        double* block = &O[k * p * n];
        cps_matrix_multiply(block, sys->C.data, A_pow, p, n, n);
    }

    int rank = cps_matrix_rank(O, m_rows, n, 1e-8);
    free(O); free(A_pow); free(temp);
    return rank == n;
}

/* ============================================================================
 * Controllability Check (L4: Fundamental Laws)
 *
 * Theorem (Kalman, 1960): System (A,B) is controllable iff
 *   rank([B, AB, A^2B, ..., A^{n-1}B]) = n
 *
 * CPS Security significance: Controllability determines the
 * attacker's ability to drive the system to unsafe states via
 * actuator compromise. The attacker's controllability Gramian
 * (through Gamma_a) bounds the reachable unsafe region.
 * ============================================================================ */

bool cps_is_controllable(const CPSDynamicalSystem* sys) {
    int n = sys->n_states;
    int m = sys->n_inputs;
    if (n <= 0 || m <= 0) return false;

    int m_cols = n * m;
    double* Co = (double*)calloc((size_t)(n * m_cols), sizeof(double));
    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp = (double*)malloc((size_t)(n * n) * sizeof(double));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Co[i * m_cols + j] = sys->B.data[i * m + j];

    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;

    for (int k = 1; k < n; k++) {
        cps_matrix_multiply(temp, A_pow, sys->A.data, n, n, n);
        memcpy(A_pow, temp, (size_t)(n * n) * sizeof(double));
        /* Compute A^k * B and place into block k of Co
         * Co has n rows × (n*m) cols; block k occupies cols k*m..k*m+m-1
         * Element at (i, k*m+j) = flat index i*(n*m) + k*m + j */
        double* AB = (double*)malloc((size_t)(n * m) * sizeof(double));
        cps_matrix_multiply(AB, A_pow, sys->B.data, n, n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                Co[i * m_cols + k * m + j] = AB[i * m + j];
        free(AB);
    }

    int rank = cps_matrix_rank(Co, n, m_cols, 1e-8);
    free(Co); free(A_pow); free(temp);
    return rank == n;
}

/* ============================================================================
 * Observability Gramian (L4: Fundamental Laws)
 *
 * Wo(T) = sum_{k=0}^{T-1} (A')^k * C' * C * A^k
 *
 * Quantifies output energy from initial state. Higher Gramian
 * determinant means better state distinguishability. In CPS
 * security, the Gramian appears in the fundamental limitation
 * of attack detection (Mo & Sinopoli, 2010, Theorem 2).
 * ============================================================================ */

double cps_observability_gramian(const CPSDynamicalSystem* sys,
                                  double* gramian, int horizon) {
    int n = sys->n_states;
    int p = sys->n_outputs;
    if (n <= 0 || p <= 0 || horizon <= 0) return -CPS_INF;

    memset(gramian, 0, (size_t)(n * n) * sizeof(double));
    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* A_pow_T = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* CT = (double*)malloc((size_t)(n * p) * sizeof(double));
    double* C_CA_k = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp = (double*)malloc((size_t)(p * n) * sizeof(double));
    double* temp2 = (double*)malloc((size_t)(n * n) * sizeof(double));

    cps_matrix_transpose(CT, sys->C.data, p, n);
    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;

    for (int k = 0; k < horizon; k++) {
        cps_matrix_multiply(temp, sys->C.data, A_pow, p, n, n);
        cps_matrix_multiply(temp2, CT, temp, n, p, n);
        cps_matrix_transpose(A_pow_T, A_pow, n, n);
        cps_matrix_multiply(C_CA_k, A_pow_T, temp2, n, n, n);

        for (int i = 0; i < n * n; i++) gramian[i] += C_CA_k[i];

        cps_matrix_multiply(temp2, A_pow, sys->A.data, n, n, n);
        memcpy(A_pow, temp2, (size_t)(n * n) * sizeof(double));
    }

    double log_det;
    if (n == 1) {
        log_det = log(gramian[0] + CPS_EPS);
    } else if (n == 2) {
        double det = gramian[0]*gramian[3] - gramian[1]*gramian[2];
        log_det = (det > CPS_EPS) ? log(det) : -CPS_INF;
    } else {
        double trace_log = 0.0;
        for (int i = 0; i < n; i++)
            trace_log += log(fabs(gramian[i * n + i]) + 1.0);
        log_det = trace_log;
    }

    free(A_pow); free(A_pow_T); free(CT);
    free(C_CA_k); free(temp); free(temp2);
    return log_det;
}

/* ============================================================================
 * Controllability Gramian (L4: Fundamental Laws)
 *
 * Wc(T) = sum_{k=0}^{T-1} A^k * B * B' * (A')^k
 *
 * Characterizes minimum control energy to reach a target state.
 * In CPS security, the attacker's controllability Gramian (through
 * Gamma_a) bounds the set of reachable states under attack.
 * ============================================================================ */

double cps_controllability_gramian(const CPSDynamicalSystem* sys,
                                    double* gramian, int horizon) {
    int n = sys->n_states;
    int m = sys->n_inputs;
    if (n <= 0 || m <= 0 || horizon <= 0) return -CPS_INF;

    memset(gramian, 0, (size_t)(n * n) * sizeof(double));
    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* AB = (double*)malloc((size_t)(n * m) * sizeof(double));
    double* BAT = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp_nn = (double*)malloc((size_t)(n * n) * sizeof(double));

    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;

    for (int k = 0; k < horizon; k++) {
        cps_matrix_multiply(AB, A_pow, sys->B.data, n, n, m);
        double* AB_T = (double*)malloc((size_t)(m * n) * sizeof(double));
        cps_matrix_transpose(AB_T, AB, n, m);
        cps_matrix_multiply(BAT, AB, AB_T, n, m, n);
        free(AB_T);

        for (int i = 0; i < n * n; i++) gramian[i] += BAT[i];

        cps_matrix_multiply(temp_nn, A_pow, sys->A.data, n, n, n);
        memcpy(A_pow, temp_nn, (size_t)(n * n) * sizeof(double));
    }

    double log_det;
    if (n == 1) {
        log_det = log(gramian[0] + CPS_EPS);
    } else if (n == 2) {
        double det = gramian[0]*gramian[3] - gramian[1]*gramian[2];
        log_det = (det > CPS_EPS) ? log(det) : -CPS_INF;
    } else {
        double tlog = 0.0;
        for (int i = 0; i < n; i++)
            tlog += log(fabs(gramian[i * n + i]) + 1.0);
        log_det = tlog;
    }

    free(A_pow); free(AB); free(BAT); free(temp_nn);
    return log_det;
}

/* ============================================================================
 * CPS Security System Lifecycle (L1: Definitions)
 * ============================================================================ */

CPSSecuritySystem* cps_security_create(int n_states, int n_inputs,
                                        int n_outputs) {
    CPSSecuritySystem* sys = (CPSSecuritySystem*)calloc(1,
                                    sizeof(CPSSecuritySystem));
    if (!sys) return NULL;

    sys->n_states = n_states;
    sys->n_inputs = n_inputs;
    sys->n_measurements = n_outputs;
    sys->time_step = 0.01;
    sys->security_state = CPS_SECURE_NORMAL;

    sys->plant = (CPSDynamicalSystem*)calloc(1, sizeof(CPSDynamicalSystem));
    if (!sys->plant) { free(sys); return NULL; }
    sys->plant->n_states = n_states;
    sys->plant->n_inputs = n_inputs;
    sys->plant->n_outputs = n_outputs;
    sys->plant->A.data = vec_alloc(n_states * n_states);
    sys->plant->A.rows = n_states; sys->plant->A.cols = n_states;
    sys->plant->A.owns_data = 1;
    sys->plant->B.data = vec_alloc(n_states * n_inputs);
    sys->plant->B.rows = n_states; sys->plant->B.cols = n_inputs;
    sys->plant->B.owns_data = 1;
    sys->plant->C.data = vec_alloc(n_outputs * n_states);
    sys->plant->C.rows = n_outputs; sys->plant->C.cols = n_states;
    sys->plant->C.owns_data = 1;
    sys->plant->D.data = vec_alloc(n_outputs * n_inputs);
    sys->plant->D.rows = n_outputs; sys->plant->D.cols = n_inputs;
    sys->plant->D.owns_data = 1;
    sys->plant->Q_scale = 0.01;
    sys->plant->R_scale = 0.01;

    sys->attack = (CPSAttackModel*)calloc(1, sizeof(CPSAttackModel));
    sys->attack->type = CPS_ATTACK_NONE;
    sys->attack->end_time = CPS_INF;
    sys->attack->signal_capacity = 256;
    sys->attack->attack_signal = vec_alloc(256);

    sys->detector = (CPSDetector*)calloc(1, sizeof(CPSDetector));
    sys->detector->method = CPS_DETECT_CHI2;
    sys->detector->threshold = 3.841;
    sys->detector->history_capacity = 256;
    sys->detector->residual_history = vec_alloc(256);

    sys->resilient = (CPSResilientController*)calloc(1,
                                    sizeof(CPSResilientController));
    sys->resilient->strategy = CPS_RESILIENT_HOLD;
    sys->resilient->input_dim = n_inputs;
    sys->resilient->safe_input = vec_alloc(n_inputs);
    sys->resilient->fallback_input = vec_alloc(n_inputs);
    sys->resilient->current_input = vec_alloc(n_inputs);
    sys->resilient->buffer_capacity = 256;
    sys->resilient->measurement_buffer = vec_alloc(256 * n_outputs);

    sys->state = vec_alloc(n_states);
    sys->true_state = vec_alloc(n_states);
    sys->measurement = vec_alloc(n_outputs);
    sys->control_input = vec_alloc(n_inputs);

    sys->log_capacity = CPS_DEFAULT_LOG_CAPACITY;
    sys->state_log = vec_alloc((size_t)n_states * CPS_DEFAULT_LOG_CAPACITY);
    sys->residual_log = vec_alloc(CPS_DEFAULT_LOG_CAPACITY);
    sys->alarm_log = (int*)calloc((size_t)CPS_DEFAULT_LOG_CAPACITY, sizeof(int));

    return sys;
}

void cps_security_free(CPSSecuritySystem* sys) {
    if (!sys) return;
    if (sys->plant) {
        free(sys->plant->A.data); free(sys->plant->B.data);
        free(sys->plant->C.data); free(sys->plant->D.data);
        free(sys->plant->observable_indices);
        free(sys->plant->controllable_indices);
        free(sys->plant);
    }
    if (sys->attack) {
        free(sys->attack->attack_signal);
        free(sys->attack->Gamma_a); free(sys->attack->Gamma_y);
        free(sys->attack);
    }
    if (sys->detector) {
        free(sys->detector->residual_history);
        free(sys->detector);
    }
    if (sys->resilient) {
        free(sys->resilient->safe_input);
        free(sys->resilient->fallback_input);
        free(sys->resilient->current_input);
        free(sys->resilient->measurement_buffer);
        free(sys->resilient->tightened_lb);
        free(sys->resilient->tightened_ub);
        free(sys->resilient->active_sensors);
        free(sys->resilient->active_actuators);
        free(sys->resilient);
    }
    free(sys->state); free(sys->true_state);
    free(sys->measurement); free(sys->control_input);
    free(sys->state_log); free(sys->residual_log); free(sys->alarm_log);
    free(sys);
}

/* ============================================================================
 * System Configuration (L2: Core Concepts)
 * ============================================================================ */

void cps_set_state_matrix(CPSSecuritySystem* sys, const double* A_data,
                           int rows, int cols) {
    if (!sys || !A_data) return;
    free(sys->plant->A.data);
    sys->plant->A.rows = rows; sys->plant->A.cols = cols;
    sys->plant->n_states = rows;
    sys->plant->A.data = vec_alloc(rows * cols);
    mat_data_copy(sys->plant->A.data, A_data, rows, cols);
    sys->n_states = rows;
    /* Reallocate state vectors for new dimension */
    free(sys->state); sys->state = vec_alloc(rows);
    free(sys->true_state); sys->true_state = vec_alloc(rows);
}

void cps_set_input_matrix(CPSSecuritySystem* sys, const double* B_data,
                           int rows, int cols) {
    if (!sys || !B_data) return;
    free(sys->plant->B.data);
    sys->plant->B.rows = rows; sys->plant->B.cols = cols;
    sys->plant->n_inputs = cols;
    sys->plant->B.data = vec_alloc(rows * cols);
    mat_data_copy(sys->plant->B.data, B_data, rows, cols);
    sys->n_inputs = cols;
    free(sys->control_input); sys->control_input = vec_alloc(cols);
    if (sys->resilient) {
        free(sys->resilient->safe_input);
        free(sys->resilient->fallback_input);
        free(sys->resilient->current_input);
        sys->resilient->safe_input = vec_alloc(cols);
        sys->resilient->fallback_input = vec_alloc(cols);
        sys->resilient->current_input = vec_alloc(cols);
        sys->resilient->input_dim = cols;
    }
}

void cps_set_output_matrix(CPSSecuritySystem* sys, const double* C_data,
                            int rows, int cols) {
    if (!sys || !C_data) return;
    free(sys->plant->C.data);
    sys->plant->C.rows = rows; sys->plant->C.cols = cols;
    sys->plant->n_outputs = rows;
    sys->plant->C.data = vec_alloc(rows * cols);
    mat_data_copy(sys->plant->C.data, C_data, rows, cols);
    sys->n_measurements = rows;
    free(sys->measurement); sys->measurement = vec_alloc(rows);
}

void cps_set_noise_covariances(CPSSecuritySystem* sys, double Q_scale,
                                double R_scale) {
    if (!sys) return;
    sys->plant->Q_scale = Q_scale;
    sys->plant->R_scale = R_scale;
}

void cps_set_initial_state(CPSSecuritySystem* sys, const double* x0) {
    if (!sys || !x0) return;
    memcpy(sys->state, x0, (size_t)sys->n_states * sizeof(double));
    memcpy(sys->true_state, x0, (size_t)sys->n_states * sizeof(double));
}

/* ============================================================================
 * System Dynamics (L3: Mathematical Structures)
 *
 * Discrete-time LTI system with attack:
 *   x[k+1] = A*x[k] + B*u[k] + Gamma_a*a[k] + w[k]
 *   y[k]   = C*x[k] + Gamma_y*a[k] + v[k]
 * where w ~ N(0,Q), v ~ N(0,R)
 * ============================================================================ */

void cps_evolve_true(CPSSecuritySystem* sys) {
    int n = sys->n_states;
    int m = sys->n_inputs;
    if (n <= 0) return;

    double* Ax = vec_alloc(n);
    double* Bu = vec_alloc(n);
    double* noise = vec_alloc(n);

    cps_matrix_multiply(Ax, sys->plant->A.data, sys->true_state, n, n, 1);
    cps_matrix_multiply(Bu, sys->plant->B.data, sys->control_input, n, m, 1);
    cps_gaussian_sample(noise, n, sqrt(sys->plant->Q_scale));

    for (int i = 0; i < n; i++)
        sys->true_state[i] = Ax[i] + Bu[i] + noise[i];

    if (sys->attack->is_active
        && sys->attack->target == CPS_TARGET_ACTUATOR
        && sys->attack->Gamma_a) {
        double att = (sys->attack->attack_signal
            && sys->attack->signal_length > 0)
            ? sys->attack->attack_signal[sys->step_count
                % sys->attack->signal_length]
            : sys->attack->magnitude;
        for (int i = 0; i < n; i++)
            sys->true_state[i] += att
                * sys->attack->Gamma_a[i * sys->attack->n_attack_inputs];
    }

    free(Ax); free(Bu); free(noise);
}

void cps_measure(CPSSecuritySystem* sys) {
    int n = sys->n_states;
    int p = sys->n_measurements;
    if (n <= 0 || p <= 0) return;

    double* noise = vec_alloc(p);
    cps_matrix_multiply(sys->measurement, sys->plant->C.data,
                        sys->true_state, p, n, 1);
    cps_gaussian_sample(noise, p, sqrt(sys->plant->R_scale));

    if (sys->attack->is_active
        && (sys->attack->target == CPS_TARGET_SENSOR
            || sys->attack->target == CPS_TARGET_NETWORK)) {
        double att = (sys->attack->attack_signal
            && sys->attack->signal_length > 0)
            ? sys->attack->attack_signal[sys->step_count
                % sys->attack->signal_length]
            : sys->attack->magnitude;
        int idx = sys->attack->target_index;
        if (idx >= 0 && idx < p) {
            sys->measurement[idx] += att;
        } else {
            for (int i = 0; i < p; i++)
                sys->measurement[i] += att;
        }
    }

    for (int i = 0; i < p; i++)
        sys->measurement[i] += noise[i];

    free(noise);
}

void cps_step(CPSSecuritySystem* sys, const double* u) {
    if (!sys) return;

    if (u && sys->n_inputs > 0)
        memcpy(sys->control_input, u,
               (size_t)sys->n_inputs * sizeof(double));

    cps_measure(sys);

    double* Ax = vec_alloc(sys->n_states);
    double* Bu = vec_alloc(sys->n_states);
    cps_matrix_multiply(Ax, sys->plant->A.data, sys->state,
                        sys->n_states, sys->n_states, 1);
    cps_matrix_multiply(Bu, sys->plant->B.data, sys->control_input,
                        sys->n_states, sys->n_inputs, 1);
    for (int i = 0; i < sys->n_states; i++)
        sys->state[i] = Ax[i] + Bu[i];
    free(Ax); free(Bu);

    double* y_hat = vec_alloc(sys->n_measurements);
    cps_matrix_multiply(y_hat, sys->plant->C.data, sys->state,
                        sys->n_measurements, sys->n_states, 1);
    double residual_norm = 0.0;
    for (int i = 0; i < sys->n_measurements; i++) {
        double r = sys->measurement[i] - y_hat[i];
        residual_norm += r * r;
    }
    residual_norm = sqrt(residual_norm);

    if (sys->log_length < sys->log_capacity) {
        sys->residual_log[sys->log_length] = residual_norm;
        sys->state_log[sys->log_length * sys->n_states] = sys->state[0];
    }

    cps_evolve_true(sys);

    sys->current_time += sys->time_step;
    sys->step_count++;
    if (sys->log_length < sys->log_capacity) sys->log_length++;

    if (sys->attack->is_active
        && sys->current_time >= sys->attack->end_time)
        sys->attack->is_active = 0;

    free(y_hat);
}

/* ============================================================================
 * Attack Configuration (L1: Definitions, L2: Core Concepts)
 * ============================================================================ */

void cps_attack_configure(CPSSecuritySystem* sys, CPSAttackType type,
                           CPSAttackTarget target, double magnitude,
                           double start_time) {
    if (!sys || !sys->attack) return;
    sys->attack->type = type;
    sys->attack->target = target;
    sys->attack->magnitude = magnitude;
    sys->attack->start_time = start_time;
    sys->attack->end_time = CPS_INF;
    sys->attack->stealthiness = 0.0;
}

void cps_attack_start(CPSSecuritySystem* sys) {
    if (!sys || !sys->attack) return;
    sys->attack->is_active = 1;
    sys->attack->start_time = sys->current_time;
    sys->security_state = CPS_SECURE_ATTACKED;
}

void cps_attack_stop(CPSSecuritySystem* sys) {
    if (!sys || !sys->attack) return;
    sys->attack->is_active = 0;
    sys->attack->end_time = sys->current_time;
    if (sys->security_state == CPS_SECURE_ATTACKED)
        sys->security_state = CPS_SECURE_RECOVERING;
}

void cps_attack_set_target_index(CPSSecuritySystem* sys, int index) {
    if (!sys || !sys->attack) return;
    sys->attack->target_index = index;
}

void cps_attack_generate_signal(CPSSecuritySystem* sys) {
    if (!sys || !sys->attack) return;
    CPSAttackModel* att = sys->attack;
    if (att->signal_capacity <= 0) {
        att->signal_capacity = 256;
        free(att->attack_signal);
        att->attack_signal = vec_alloc(256);
    }
    int len = att->signal_capacity;
    att->signal_length = len;

    for (int i = 0; i < len; i++) {
        double t = (double)i;
        switch (att->type) {
            case CPS_ATTACK_FDI:
                att->attack_signal[i] = att->magnitude
                    * (1.0 + 0.05 * sin(0.05 * t));
                break;
            case CPS_ATTACK_BIAS:
                att->attack_signal[i] = att->magnitude;
                break;
            case CPS_ATTACK_SURGE:
                att->attack_signal[i] = (t < len * 0.2)
                    ? att->magnitude : 0.0;
                break;
            case CPS_ATTACK_REPLAY:
                att->attack_signal[i] = att->magnitude
                    * sin(0.1 * (t - len * 0.25));
                break;
            case CPS_ATTACK_COVERT:
                att->attack_signal[i] = att->magnitude * 0.1
                    * cos(0.3 * t) * exp(-0.005 * t);
                break;
            case CPS_ATTACK_ZERO_DYNAMICS:
                att->attack_signal[i] = att->magnitude
                    * exp(-0.1 * t);
                break;
            case CPS_ATTACK_DOS:
                att->attack_signal[i] = (fmod(t, 20.0) < 10.0)
                    ? att->magnitude : 0.0;
                break;
            default:
                att->attack_signal[i] = 0.0;
                break;
        }
    }
}

bool cps_attack_is_active(CPSSecuritySystem* sys) {
    return sys && sys->attack && sys->attack->is_active;
}

void cps_attack_set_stealthiness(CPSSecuritySystem* sys, double stealth) {
    if (!sys || !sys->attack) return;
    if (stealth < 0.0) stealth = 0.0;
    if (stealth > 1.0) stealth = 1.0;
    sys->attack->stealthiness = stealth;
    sys->attack->magnitude *= (1.0 - 0.9 * stealth);
}

CPSSecurityState cps_security_get_state(CPSSecuritySystem* sys) {
    return sys ? sys->security_state : CPS_SECURE_COMPROMISED;
}

const char* cps_attack_type_name(CPSAttackType type) {
    switch (type) {
        case CPS_ATTACK_NONE:          return "None";
        case CPS_ATTACK_DOS:           return "Denial-of-Service";
        case CPS_ATTACK_FDI:           return "False Data Injection";
        case CPS_ATTACK_REPLAY:        return "Replay";
        case CPS_ATTACK_BIAS:          return "Bias Injection";
        case CPS_ATTACK_COVERT:        return "Covert";
        case CPS_ATTACK_SURGE:         return "Surge";
        case CPS_ATTACK_ZERO_DYNAMICS: return "Zero Dynamics";
        default:                       return "Unknown";
    }
}

const char* cps_security_state_name(CPSSecurityState state) {
    switch (state) {
        case CPS_SECURE_NORMAL:     return "Normal";
        case CPS_SECURE_SUSPICIOUS: return "Suspicious";
        case CPS_SECURE_ATTACKED:   return "Under Attack";
        case CPS_SECURE_DEGRADED:   return "Degraded";
        case CPS_SECURE_RECOVERING: return "Recovering";
        case CPS_SECURE_COMPROMISED:return "Compromised";
        default:                    return "Unknown";
    }
}
