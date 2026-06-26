/* cloud_control_core.c - Cloud Control System Core Implementation */
#include "cloud_control_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* LCG pseudo-random generator */
static unsigned int lcg_state = 1;
static void lcg_seed(unsigned int s) { lcg_state = s; }
static double lcg_rand(void) {
    lcg_state = (1103515245 * lcg_state + 12345) & 0x7fffffff;
    return (double)lcg_state / 0x7fffffff;
}

/* mat_vec_mul: y = A*x where A is rows x cols, row-major */
static void mat_vec_mul(const double *A, const double *x, double *y,
                         int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        y[i] = 0.0;
        for (int j = 0; j < cols; j++) y[i] += A[i*cols + j] * x[j];
    }
}
static double vec_norm(const double *v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i]*v[i];
    return sqrt(s);
}
static double vec_dot(const double *a, const double *b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i]*b[i];
    return s;
}
/* Eigenvalue via power iteration with Rayleigh quotient.
 * Available as alternative to QR for large sparse systems. */
__attribute__((unused))
static double power_iteration(const double *A, int n, double *evec,
                               int max_iter, double tol) {
    double *v = (double*)calloc((size_t)n, sizeof(double));
    double *Av = (double*)calloc((size_t)n, sizeof(double));
    if (!v || !Av) { free(v); free(Av); return 0.0; }
    for (int i = 0; i < n; i++) v[i] = 1.0 / sqrt((double)n);
    double lambda = 0.0, lambda_old = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        mat_vec_mul(A, v, Av, n, n);
        lambda_old = lambda;
        lambda = vec_dot(v, Av, n);
        double norm = vec_norm(Av, n);
        if (norm < 1e-15) break;
        for (int i = 0; i < n; i++) v[i] = Av[i] / norm;
        if (fabs(lambda - lambda_old) < tol) break;
    }
    if (evec) for (int i = 0; i < n; i++) evec[i] = v[i];
    free(v); free(Av);
    return lambda;
}

/* QR algorithm with Wilkinson shift for max real eigenvalue */
static double max_eigenvalue_real(const double *A_in, int n) {
    if (n <= 0) return 0.0;
    if (n == 1) return A_in[0];
    double *A = (double*)calloc((size_t)(n*n), sizeof(double));
    double *Q = (double*)calloc((size_t)(n*n), sizeof(double));
    double *R = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!A || !Q || !R) { free(A); free(Q); free(R); return 0.0; }
    for (int i = 0; i < n*n; i++) A[i] = A_in[i];
    for (int iter = 0; iter < 200; iter++) {
        double shift = A[(n-1)*n + (n-1)];
        if (n >= 2) {
            double a = A[(n-2)*n+(n-2)], b = A[(n-2)*n+(n-1)];
            double c_ = A[(n-1)*n+(n-2)], d = A[(n-1)*n+(n-1)];
            double tr = a + d, det = a*d - b*c_;
            double disc = tr*tr - 4.0*det;
            if (disc >= 0) {
                double r1 = (tr + sqrt(disc))/2.0, r2 = (tr - sqrt(disc))/2.0;
                shift = (fabs(r1-d) < fabs(r2-d)) ? r1 : r2;
            }
        }
        for (int i = 0; i < n; i++) A[i*n+i] -= shift;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) R[j*n+i] = 0.0;
            for (int i = 0; i < n; i++) Q[i*n+j] = A[i*n+j];
            for (int k = 0; k < j; k++) {
                double dot = 0.0;
                for (int i = 0; i < n; i++) dot += Q[i*n+k] * A[i*n+j];
                R[k*n+j] = dot;
                for (int i = 0; i < n; i++) Q[i*n+j] -= dot * Q[i*n+k];
            }
            double norm = 0.0;
            for (int i = 0; i < n; i++) norm += Q[i*n+j]*Q[i*n+j];
            norm = sqrt(norm);
            if (norm > 1e-12) {
                R[j*n+j] = norm;
                for (int i = 0; i < n; i++) Q[i*n+j] /= norm;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) sum += R[i*n+k] * Q[k*n+j];
                A[i*n+j] = sum;
            }
            A[i*n+i] += shift;
        }
        double off_diag = 0.0;
        for (int i = 1; i < n; i++) off_diag += fabs(A[i*n+(i-1)]);
        if (off_diag < 1e-12) break;
    }
    double max_re = -1e300;
    for (int i = 0; i < n; i++)
        if (A[i*n+i] > max_re) max_re = A[i*n+i];
    free(A); free(Q); free(R);
    return max_re;
}
/* -------- Lifecycle API -------- */

CloudControlSystem* ccs_create(const char *id, int n, int m, int p,
                                CloudControlMode mode) {
    if (!id || n <= 0 || n > CCS_MAX_STATES ||
        m <= 0 || m > CCS_MAX_INPUTS || p <= 0 || p > CCS_MAX_OUTPUTS)
        return NULL;
    CloudControlSystem *ccs = (CloudControlSystem*)calloc(1, sizeof(CloudControlSystem));
    if (!ccs) return NULL;
    strncpy(ccs->id, id, sizeof(ccs->id)-1);
    ccs->n = n; ccs->m = m; ccs->p = p;
    ccs->mode = mode;
    ccs->state = CC_STATE_INITIALIZING;
    ccs->plant_state.dim = n;
    ccs->plant_state.input_dim = m;
    ccs->plant_state.output_dim = p;
    for (int i = 0; i < n; i++) ccs->A[i][i] = -1.0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) ccs->B[j][i] = 1.0;
    for (int i = 0; i < p; i++) ccs->C[i][i] = 1.0;
    clock_gettime(CLOCK_REALTIME, &ccs->created_at);
    ccs->last_update = ccs->created_at;
    return ccs;
}

void ccs_free(CloudControlSystem *ccs) { free(ccs); }

int ccs_set_plant_model(CloudControlSystem *ccs,
                         const double *A, const double *B,
                         const double *C, const double *D) {
    if (!ccs || !A || !B || !C) return -1;
    int n = ccs->n, m = ccs->m, p = ccs->p;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ccs->A[i][j] = A[i*n + j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            ccs->B[i][j] = B[i*m + j];
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            ccs->C[i][j] = C[i*n + j];
    if (D)
        for (int i = 0; i < p; i++)
            for (int j = 0; j < m; j++)
                ccs->D[i][j] = D[i*m + j];
    return 0;
}

int ccs_set_controller(CloudControlSystem *ccs,
                        const double *K, const double *L) {
    if (!ccs || !K || !L) return -1;
    int n = ccs->n, m = ccs->m, p = ccs->p;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            ccs->K[i][j] = K[i*n + j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            ccs->L[i][j] = L[i*p + j];
    return 0;
}

int ccs_set_qos(CloudControlSystem *ccs, const QoSProfile *qos) {
    if (!ccs || !qos) return -1;
    ccs->qos = *qos; return 0;
}

int ccs_set_reference(CloudControlSystem *ccs, const double *ref) {
    if (!ccs || !ref) return -1;
    for (int i = 0; i < ccs->n; i++) ccs->reference.x[i] = ref[i];
    return 0;
}
/* -------- Control Loop Execution -------- */

int ccs_compute_control(CloudControlSystem *ccs, double *u_out) {
    if (!ccs || !u_out) return -1;
    int n = ccs->n, m = ccs->m;
    double err[CCS_MAX_STATES];
    for (int i = 0; i < n; i++)
        err[i] = ccs->x_hat[i] - ccs->reference.x[i];
    for (int i = 0; i < m; i++) {
        u_out[i] = 0.0;
        for (int j = 0; j < n; j++)
            u_out[i] -= ccs->K[i][j] * err[j];
    }
    return 0;
}

int ccs_update_observer(CloudControlSystem *ccs,
                         const double *y, double timestamp) {
    if (!ccs || !y) return -1;
    if (timestamp < ccs->plant_state.t - 1.0) return -1;
    int n = ccs->n, p = ccs->p;
    double y_pred[CCS_MAX_OUTPUTS];
    for (int i = 0; i < p; i++) {
        y_pred[i] = 0.0;
        for (int j = 0; j < n; j++)
            y_pred[i] += ccs->C[i][j] * ccs->x_hat[j];
    }
    double innov[CCS_MAX_OUTPUTS];
    for (int i = 0; i < p; i++) innov[i] = y[i] - y_pred[i];
    for (int i = 0; i < n; i++) {
        double corr = 0.0;
        for (int j = 0; j < p; j++) corr += ccs->L[i][j] * innov[j];
        ccs->x_hat[i] += corr;
    }
    return 0;
}

int ccs_apply_control(CloudControlSystem *ccs, const double *u, double dt) {
    if (!ccs || !u || dt <= 0.0) return -1;
    int n = ccs->n, m = ccs->m, p = ccs->p;
    double *x = ccs->plant_state.x;
    double dx[CCS_MAX_STATES] = {0};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) dx[i] += ccs->A[i][j] * x[j];
        for (int j = 0; j < m; j++) dx[i] += ccs->B[i][j] * u[j];
    }
    for (int i = 0; i < n; i++) x[i] += dt * dx[i];
    /* Observer prediction step */
    double obs_pred[CCS_MAX_STATES] = {0};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            obs_pred[i] += ccs->A[i][j] * ccs->x_hat[j];
        for (int j = 0; j < m; j++)
            obs_pred[i] += ccs->B[i][j] * u[j];
    }
    for (int i = 0; i < n; i++) ccs->x_hat[i] += dt * obs_pred[i];
    /* Update output */
    for (int i = 0; i < p; i++) {
        ccs->plant_state.y[i] = 0.0;
        for (int j = 0; j < n; j++)
            ccs->plant_state.y[i] += ccs->C[i][j] * x[j];
    }
    ccs->plant_state.t += dt;
    ccs->plant_state.seq++;
    if (ccs->history_count < CCS_MAX_HISTORY) {
        ccs->history_t[ccs->history_count] = ccs->plant_state.t;
        ccs->history_y[ccs->history_count] = x[0];
        ccs->history_count++;
    }
    return 0;
}

int ccs_step(CloudControlSystem *ccs, const double *measurement,
              double ts, double dt, double delay_us) {
    if (!ccs) return -1;
    (void)ccs->n; (void)ccs->m;
    double u[CCS_MAX_INPUTS] = {0};
    if (measurement) ccs_update_observer(ccs, measurement, ts);
    ccs_compute_control(ccs, u);
    ccs_apply_control(ccs, u, dt);
    ccs->avg_rtt_us = 0.95 * ccs->avg_rtt_us + 0.05 * delay_us;
    if (delay_us > ccs->max_rtt_us) ccs->max_rtt_us = delay_us;
    ccs->packets_sent++;
    ccs->total_cycles++;
    return 0;
}
/* -------- Performance Metrics -------- */

int ccs_compute_performance(CloudControlSystem *ccs) {
    if (!ccs || ccs->history_count < 2) return -1;
    int hc = ccs->history_count;
    double *t = ccs->history_t, *y = ccs->history_y;
    double ref = ccs->reference.x[0];
    double y_final = y[hc-1];
    double band = 0.05 * fabs(ref > 0.01 ? ref : 1.0);
    ccs->settling_time = 0.0;
    for (int i = hc-1; i >= 0; i--) {
        if (fabs(y[i] - y_final) > band) {
            ccs->settling_time = (i+1 < hc) ? t[i+1] : t[hc-1]; break;
        }
    }
    double y_max = y[0], step_mag = fabs(ref) > 0.01 ? fabs(ref) : 1.0;
    for (int i = 1; i < hc; i++) if (y[i] > y_max) y_max = y[i];
    ccs->overshoot_percent = (y_max - y_final) / step_mag * 100.0;
    if (ccs->overshoot_percent < 0) ccs->overshoot_percent = 0;
    ccs->ise = 0.0; ccs->iae = 0.0;
    for (int i = 0; i < hc; i++) {
        double e = y[i] - ref, dt_i = (i > 0) ? (t[i]-t[i-1]) : t[0];
        ccs->ise += e*e*dt_i; ccs->iae += fabs(e)*dt_i;
    }
    ccs->steady_state_error = y_final - ref;
    return 0;
}

int ccs_is_stable(const CloudControlSystem *ccs, double delay_us) {
    if (!ccs) return -1;
    int n = ccs->n, m = ccs->m;
    double *A_cl = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!A_cl) return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            A_cl[i*n+j] = ccs->A[i][j];
            for (int k = 0; k < m; k++)
                A_cl[i*n+j] -= ccs->B[i][k] * ccs->K[k][j];
        }
    double max_re = max_eigenvalue_real(A_cl, n);
    free(A_cl);
    (void)delay_us;
    return (max_re < -1e-9) ? 1 : 0;
}

double ccs_max_allowable_delay(const CloudControlSystem *ccs) {
    if (!ccs || ccs->n <= 0) return -1.0;
    int n = ccs->n, m = ccs->m;
    double *A_cl = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!A_cl) return -1.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            A_cl[i*n+j] = ccs->A[i][j];
            for (int k = 0; k < m; k++)
                A_cl[i*n+j] -= ccs->B[i][k] * ccs->K[k][j];
        }
    double max_re = max_eigenvalue_real(A_cl, n);
    free(A_cl);
    if (max_re >= 0) return 0.0;
    double tau_max_s = -1.0 / max_re;
    return tau_max_s * 1e6 * 0.1;
}

/* -------- Mode & String API -------- */

int ccs_switch_mode(CloudControlSystem *ccs, CloudControlMode new_mode) {
    if (!ccs) return -1;
    if (ccs->mode == new_mode) return 0;
    switch (new_mode) {
    case CC_MODE_EDGE_ONLY:
        ccs->edge_node.is_active = 1; ccs->edge_node.cloud_connected = 0;
        ccs->state = CC_STATE_FALLBACK; break;
    case CC_MODE_CLOUD_ONLY:
        ccs->edge_node.is_active = 0; ccs->edge_node.cloud_connected = 1;
        ccs->state = CC_STATE_RUNNING; break;
    case CC_MODE_COLLABORATIVE:
        ccs->edge_node.is_active = 1; ccs->edge_node.cloud_connected = 1;
        ccs->state = CC_STATE_RUNNING; break;
    default: break;
    }
    ccs->mode = new_mode;
    return 0;
}

const char* ccs_mode_string(CloudControlMode mode) {
    switch (mode) {
    case CC_MODE_EDGE_ONLY: return "EDGE_ONLY";
    case CC_MODE_CLOUD_ONLY: return "CLOUD_ONLY";
    case CC_MODE_COLLABORATIVE: return "COLLABORATIVE";
    case CC_MODE_HIERARCHICAL: return "HIERARCHICAL";
    case CC_MODE_HYBRID: return "HYBRID";
    case CC_MODE_REDUNDANT: return "REDUNDANT";
    default: return "UNKNOWN";
    }
}

const char* ccs_state_string(CloudControlState state) {
    switch (state) {
    case CC_STATE_UNINITIALIZED: return "UNINITIALIZED";
    case CC_STATE_INITIALIZING: return "INITIALIZING";
    case CC_STATE_RUNNING: return "RUNNING";
    case CC_STATE_DEGRADED: return "DEGRADED";
    case CC_STATE_FALLBACK: return "FALLBACK";
    case CC_STATE_RECOVERING: return "RECOVERING";
    case CC_STATE_MAINTENANCE: return "MAINTENANCE";
    case CC_STATE_ERROR: return "ERROR";
    default: return "UNKNOWN";
    }
}
/* -------- Node Management -------- */

int ccs_register_cloud_node(CloudControlSystem *ccs, const CloudNode *node) {
    if (!ccs || !node) return -1;
    ccs->cloud_node = *node; return 0;
}
int ccs_register_edge_node(CloudControlSystem *ccs, const EdgeNode *node) {
    if (!ccs || !node) return -1;
    ccs->edge_node = *node; return 0;
}
int ccs_update_cloud_metrics(CloudControlSystem *ccs,
                              double cpu_util, double mem_util,
                              double resp_time_us, int queue_len) {
    if (!ccs) return -1;
    ccs->cloud_node.cpu_utilization = cpu_util;
    ccs->cloud_node.memory_utilization = mem_util;
    ccs->cloud_node.avg_response_time_us =
        0.7 * ccs->cloud_node.avg_response_time_us + 0.3 * resp_time_us;
    if (resp_time_us > ccs->cloud_node.max_response_time_us)
        ccs->cloud_node.max_response_time_us = resp_time_us;
    ccs->cloud_node.task_queue_length = queue_len;
    return 0;
}
double ccs_get_cloud_load(const CloudControlSystem *ccs) {
    if (!ccs) return 0.0;
    double cpu_l = ccs->cloud_node.cpu_utilization;
    double mem_l = ccs->cloud_node.memory_utilization;
    double q_l = (double)ccs->cloud_node.task_queue_length / CCS_MAX_TASK_QUEUE;
    double max_l = cpu_l;
    if (mem_l > max_l) max_l = mem_l;
    if (q_l > max_l) max_l = q_l;
    return max_l;
}

/* -------- Random System Initialization -------- */

int ccs_random_init(CloudControlSystem *ccs, int n, int m, int p,
                     unsigned int seed) {
    if (!ccs || n <= 0 || m <= 0 || p <= 0) return -1;
    if (n > CCS_MAX_STATES || m > CCS_MAX_INPUTS || p > CCS_MAX_OUTPUTS) return -1;
    ccs->n = n; ccs->m = m; ccs->p = p;
    lcg_seed(seed);
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) { ccs->A[i][j] = (lcg_rand() - 0.5) * 2.0; row_sum += fabs(ccs->A[i][j]); }
        }
        ccs->A[i][i] = -(row_sum + lcg_rand() + 0.1);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) ccs->B[i][j] = lcg_rand();
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++) ccs->C[i][j] = (i==j) ? 1.0 : lcg_rand() * 0.1;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) ccs->K[i][j] = (i == j%m) ? 1.0 : 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++) ccs->L[i][j] = (i==j) ? 1.0 : 0.0;
    ccs->plant_state.dim = n; ccs->plant_state.input_dim = m; ccs->plant_state.output_dim = p;
    return 0;
}

/* -------- Comparison & Reset -------- */

int ccs_compare(const CloudControlSystem *a, const CloudControlSystem *b) {
    if (!a || !b) return 0;
    if (a->n != b->n || a->m != b->m || a->p != b->p) return 0;
    if (a->mode != b->mode) return 0;
    for (int i = 0; i < a->n; i++)
        for (int j = 0; j < a->n; j++)
            if (fabs(a->A[i][j] - b->A[i][j]) > 1e-12) return 0;
    return 1;
}

void ccs_reset_metrics(CloudControlSystem *ccs) {
    if (!ccs) return;
    ccs->settling_time = 0; ccs->overshoot_percent = 0; ccs->steady_state_error = 0;
    ccs->ise = 0; ccs->iae = 0; ccs->control_effort = 0;
    ccs->total_cycles = 0; ccs->avg_rtt_us = 0; ccs->max_rtt_us = 0;
    ccs->avg_jitter_us = 0; ccs->packet_loss_rate = 0;
    ccs->packets_sent = 0; ccs->packets_lost = 0; ccs->history_count = 0;
}

/* -------- Debug Printing -------- */



void ccs_print_state(const CloudControlSystem *ccs){
    if(!ccs){printf("CCS: NULL\n");return;}
    printf("CloudControlSystem[%s] mode=%s state=%s n=%d m=%d p=%d\n",ccs->id,ccs_mode_string(ccs->mode),ccs_state_string(ccs->state),ccs->n,ccs->m,ccs->p);
    printf("  Plant x: [");
    for(int i=0;i<ccs->n&&i<6;i++)printf("%.4f%s",ccs->plant_state.x[i],i<ccs->n-1&&i<5?", ":"");
    if(ccs->n>6){printf("...");} printf("] t=%.4f\n",ccs->plant_state.t);
    printf("  Estimate: [");
    for(int i=0;i<ccs->n&&i<6;i++)printf("%.4f%s",ccs->x_hat[i],i<ccs->n-1&&i<5?", ":"");
    if(ccs->n>6){printf("...");} printf("]\n");
    printf("  RTT avg=%.0fus max=%.0fus cycles=%llu\n",ccs->avg_rtt_us,ccs->max_rtt_us,(unsigned long long)ccs->total_cycles);
}
void ccs_print_metrics(const CloudControlSystem *ccs){
    if(!ccs){printf("No metrics.\n");return;}
    printf("Performance for CCS[%s]:\n",ccs->id);
    printf("  Settling: %.4fs Overshoot: %.2f%% SS-error: %.6f\n",ccs->settling_time,ccs->overshoot_percent,ccs->steady_state_error);
    printf("  ISE: %.6f IAE: %.6f\n",ccs->ise,ccs->iae);
    printf("  Cycles: %llu Avg-RTT: %.0fus Max-RTT: %.0fus\n",(unsigned long long)ccs->total_cycles,ccs->avg_rtt_us,ccs->max_rtt_us);
}
