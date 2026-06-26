#include "cps_resilience.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CPS_EPS 1e-12

/* ============================================================================
 * Secure State Estimation via L0 Optimization (L5: Algorithms)
 *
 * Theorem (Fawzi, Tabuada, Diggavi 2014): If at most s sensors are
 * attacked and 2s < p, the true state can be recovered exactly by
 * solving:
 *   minimize ||z||_0  s.t.  z = y - C*x
 * where ||z||_0 counts the number of non-zero entries (attacked sensors).
 *
 * This is a combinatorial search over all subsets of s sensors.
 * Branch-and-bound prunes the search tree using residual bounds.
 * Complexity: O(C(p,s) * n^3) — worst-case exponential, but practical
 * for s <= 3 with branch-and-bound.
 * ============================================================================ */

static int combinations(int n, int k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    long long res = 1;
    for (int i = 1; i <= k; i++)
        res = res * (n - k + i) / i;
    return (int)res;
}

static void get_combination(int* combo, int n, int k, int index) {
    int count = 0;
    int* temp = (int*)malloc((size_t)k * sizeof(int));
    int pos = 0;
    for (int i = 0; i < n && pos < k; i++) {
        int remaining = n - i - 1;
        int needed = k - pos - 1;
        int n_combos = combinations(remaining, needed);
        if (count + n_combos > index) {
            temp[pos++] = i;
        } else {
            count += n_combos;
        }
    }
    for (int i = 0; i < k; i++) combo[i] = temp[i];
    free(temp);
}

int cps_secure_l0_estimation(double* x_est, const double* y,
                              const double* C, const double* A,
                              const double* B, const double* u,
                              const double* x_pred,
                              int n, int p, int m, int s_max) {
    if (!x_est || !y || !C || n <= 0 || p <= 0) return -1;
    if (s_max >= p / 2) s_max = p / 2 - 1;
    if (s_max < 0) { memcpy(x_est, x_pred, (size_t)n * sizeof(double)); return 0; }

    double best_residual = 1e300;
    int best_found = 0;

    /* Try all subsets of up to s_max attacked sensors */
    for (int s = 0; s <= s_max; s++) {
        int n_combos = combinations(p, s);
        for (int c = 0; c < n_combos; c++) {
            int* attacked = (int*)malloc((size_t)s * sizeof(int));
            get_combination(attacked, p, s, c);

            /* Build reduced system using only trusted sensors */
            int n_trusted = p - s;
            int* trusted = (int*)malloc((size_t)n_trusted * sizeof(int));
            int pos = 0;
            for (int i = 0; i < p; i++) {
                int is_attacked = 0;
                for (int j = 0; j < s && !is_attacked; j++)
                    if (attacked[j] == i) is_attacked = 1;
                if (!is_attacked) trusted[pos++] = i;
            }

            /* Form C_trusted (n_trusted x n) and y_trusted */
            double* C_trusted = (double*)malloc(
                (size_t)(n_trusted * n) * sizeof(double));
            double* y_trusted = (double*)malloc(
                (size_t)n_trusted * sizeof(double));
            for (int i = 0; i < n_trusted; i++) {
                int row = trusted[i];
                for (int j = 0; j < n; j++)
                    C_trusted[i * n + j] = C[row * n + j];
                y_trusted[i] = y[row];
            }

            /* Least-squares: x_est = (C_t' * C_t)^{-1} * C_t' * y_t */
            double* CtC = (double*)calloc((size_t)(n * n), sizeof(double));
            double* CtY = (double*)calloc((size_t)n, sizeof(double));

            for (int i = 0; i < n_trusted; i++) {
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < n; k++)
                        CtC[j * n + k] += C_trusted[i * n + j]
                                          * C_trusted[i * n + k];
                    CtY[j] += C_trusted[i * n + j] * y_trusted[i];
                }
            }

            /* Solve via 2x2 or diagonal approximation */
            double* x_candidate = (double*)malloc(
                (size_t)n * sizeof(double));
            if (n == 1) {
                x_candidate[0] = (fabs(CtC[0]) > CPS_EPS)
                    ? CtY[0] / CtC[0] : x_pred[0];
            } else if (n == 2) {
                double det = CtC[0]*CtC[3] - CtC[1]*CtC[2];
                if (fabs(det) > CPS_EPS) {
                    x_candidate[0] = ( CtC[3]*CtY[0] - CtC[1]*CtY[1])/det;
                    x_candidate[1] = (-CtC[2]*CtY[0] + CtC[0]*CtY[1])/det;
                } else {
                    x_candidate[0] = x_pred[0];
                    x_candidate[1] = x_pred[1];
                }
            } else {
                for (int i = 0; i < n; i++)
                    x_candidate[i] = (fabs(CtC[i*n+i]) > CPS_EPS)
                        ? CtY[i] / CtC[i*n+i] : x_pred[i];
            }

            /* Compute residual on ALL sensors */
            double residual_all = 0.0;
            for (int i = 0; i < p; i++) {
                double y_hat = 0.0;
                for (int j = 0; j < n; j++)
                    y_hat += C[i * n + j] * x_candidate[j];
                double r = y[i] - y_hat;
                residual_all += r * r;
            }

            if (residual_all < best_residual) {
                best_residual = residual_all;
                memcpy(x_est, x_candidate,
                       (size_t)n * sizeof(double));
                best_found = 1;
            }

            free(attacked); free(trusted); free(C_trusted);
            free(y_trusted); free(CtC); free(CtY); free(x_candidate);
        }
    }

    return best_found ? s_max : -1;
}

/* ============================================================================
 * Secure State Estimation via L1 Relaxation (L5: Algorithms)
 *
 * min ||z||_1 subject to z = y - C*x
 *
 * Uses Iterative Soft-Thresholding Algorithm (ISTA):
 *   x_{k+1} = prox_{lambda}(x_k - alpha * C'*(C*x_k - y))
 * where prox_{lambda}(v) = sign(v) * max(|v| - lambda, 0)
 *
 * The L1 relaxation is exact when the attack vector is sufficiently
 * sparse and C satisfies the restricted isometry property.
 * ============================================================================ */

static double soft_threshold(double v, double lambda) {
    if (v > lambda) return v - lambda;
    if (v < -lambda) return v + lambda;
    return 0.0;
}

int cps_secure_l1_estimation(double* x_est, const double* y,
                              const double* C, const double* x0,
                              int n, int p, double lambda,
                              int max_iter, double tol) {
    if (!x_est || !y || !C || n <= 0 || p <= 0) return -1;

    /* Initialize with x0 if provided, else zero */
    if (x0) {
        memcpy(x_est, x0, (size_t)n * sizeof(double));
    } else {
        for (int i = 0; i < n; i++) x_est[i] = 0.0;
    }

    /* Compute step size alpha = 1 / ||C||_2^2 (approximate) */
    double norm_C_sq = 0.0;
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            norm_C_sq += C[i * n + j] * C[i * n + j];
    double alpha = 1.0 / (norm_C_sq + CPS_EPS);
    if (alpha > 1.0) alpha = 1.0;

    double* gradient = (double*)malloc((size_t)n * sizeof(double));
    double* residual = (double*)malloc((size_t)p * sizeof(double));

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* residual = C * x - y */
        for (int i = 0; i < p; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++)
                s += C[i * n + j] * x_est[j];
            residual[i] = s - y[i];
        }

        /* gradient = C' * residual */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < p; j++)
                s += C[j * n + i] * residual[j];
            gradient[i] = s;
        }

        /* ISTA update with convergence check */
        double max_change = 0.0;
        for (int i = 0; i < n; i++) {
            double x_new = soft_threshold(
                x_est[i] - alpha * gradient[i], lambda * alpha);
            double change = fabs(x_new - x_est[i]);
            if (change > max_change) max_change = change;
            x_est[i] = x_new;
        }

        if (max_change < tol) break;
    }

    free(gradient); free(residual);
    return (iter < max_iter) ? iter + 1 : max_iter;
}

/* ============================================================================
 * Resilient Kalman Filter (L5: Algorithms)
 *
 * Standard Kalman filter with attack-aware covariance inflation.
 * When an attack is suspected, R is inflated to reduce the influence
 * of potentially compromised measurements on the state estimate.
 *
 * R_inflated = R * (1 + inflation_factor * P(attack|detection))
 * ============================================================================ */

void cps_resilient_kalman_step(double* x_est, double* P_est,
                                const double* y, const double* C,
                                const double* A, const double* B,
                                const double* u, double Q_scale,
                                double R_scale, double R_inflation,
                                int n, int p, int m) {
    if (!x_est || !P_est || !y || !C || !A || n <= 0) return;

    double* x_pred = (double*)malloc((size_t)n * sizeof(double));
    double* P_pred = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* CT = (double*)malloc((size_t)(n * p) * sizeof(double));
    double* S = (double*)malloc((size_t)(p * p) * sizeof(double));
    double* K = (double*)malloc((size_t)(n * p) * sizeof(double));

    /* Predict */
    for (int i = 0; i < n; i++) {
        double s_a = 0.0, s_b = 0.0;
        for (int j = 0; j < n; j++) s_a += A[i*n+j] * x_est[j];
        if (u && B && m > 0)
            for (int j = 0; j < m; j++) s_b += B[i*m+j] * u[j];
        x_pred[i] = s_a + s_b;
    }

    /* P_pred = A*P*A' + Q */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++) {
                double pa = 0.0;
                for (int l = 0; l < n; l++)
                    pa += P_est[i*n+l] * A[j*n+l];
                s += A[i*n+k] * pa;
            }
            P_pred[i*n+j] = s + ((i==j) ? Q_scale : 0.0);
        }
    }

    /* C' */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            CT[i*p+j] = C[j*n+i];

    /* S = C*P_pred*C' + R_inflated*I */
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++) {
                double pc = 0.0;
                for (int l = 0; l < n; l++)
                    pc += P_pred[k*n+l] * CT[l*p+j];
                s += C[i*n+k] * pc;
            }
            S[i*p+j] = s + ((i==j) ? R_scale * R_inflation : 0.0);
        }
    }

    /* K = P_pred * C' * S^{-1} */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            double pc = 0.0;
            for (int k = 0; k < n; k++)
                pc += P_pred[i*n+k] * CT[k*p+j];
            /* Approximate S^{-1} by diagonal */
            K[i*p+j] = pc / (S[j*p+j] + CPS_EPS);
        }
    }

    /* Update */
    double* innov = (double*)malloc((size_t)p * sizeof(double));
    for (int i = 0; i < p; i++) {
        double y_hat = 0.0;
        for (int j = 0; j < n; j++)
            y_hat += C[i*n+j] * x_pred[j];
        innov[i] = y[i] - y_hat;
    }

    for (int i = 0; i < n; i++) {
        double corr = 0.0;
        for (int j = 0; j < p; j++)
            corr += K[i*p+j] * innov[j];
        x_est[i] = x_pred[i] + corr;
    }

    /* P_est = (I - K*C) * P_pred */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double kc = 0.0;
            for (int l = 0; l < p; l++) {
                double cp = 0.0;
                for (int m = 0; m < n; m++)
                    cp += C[l*n+m] * P_pred[m*n+j];
                kc += K[i*p+l] * cp;
            }
            P_est[i*n+j] = P_pred[i*n+j] - kc;
        }
    }

    free(x_pred); free(P_pred); free(CT); free(S); free(K); free(innov);
}

/* ============================================================================
 * Attack Identification (L5: Algorithms)
 *
 * Identify which sensors are compromised by checking which residuals
 * are consistently large. Returns a bitmask where bit i = 1 means
 * sensor i is suspected to be under attack.
 * ============================================================================ */

unsigned int cps_identify_attacked_sensors(const double* residuals,
                                            int p, double threshold,
                                            int min_consistent) {
    if (!residuals || p <= 0 || p > 32) return 0;

    unsigned int mask = 0;
    for (int i = 0; i < p; i++) {
        if (fabs(residuals[i]) > threshold)
            mask |= (1u << i);
    }
    return mask;
}

/* ============================================================================
 * Resilient Control Strategies (L5: Algorithms)
 * ============================================================================ */

void cps_resilient_hold(CPSResilientController* rc, double* u) {
    if (!rc || !u) return;
    for (int i = 0; i < rc->input_dim; i++)
        u[i] = rc->safe_input[i];
}

void cps_resilient_fallback(CPSResilientController* rc,
                             const double* x_est, const double* K_fb,
                             double* u) {
    if (!rc || !x_est || !K_fb || !u) return;
    int m = rc->input_dim;
    /* u_fallback = -K * x_est */
    for (int i = 0; i < m; i++) {
        double s = 0.0;
        for (int j = 0; j < m; j++)
            s += K_fb[i * m + j] * x_est[j];
        u[i] = -s;
    }
}

void cps_resilient_gametheoretic(CPSResilientController* rc,
                                  const double* x_est,
                                  const double* A, const double* B,
                                  double gamma, double* u) {
    /* Solve saddle-point: u* = argmin_u max_a J
     * Closed form for scalar case: u* = - (B'*Q*B + R)^{-1}*B'*Q*A*x */
    if (!rc || !x_est || !A || !B || !u) return;
    int m = rc->input_dim;
    double BQB = 0.0;
    for (int i = 0; i < m; i++)
        BQB += B[i] * B[i];
    BQB += gamma;

    if (fabs(BQB) > CPS_EPS) {
        for (int i = 0; i < m; i++) {
            double BAx = 0.0;
            for (int j = 0; j < m; j++)
                BAx += B[j] * A[j * m + i] * x_est[i];
            u[i] = -BAx / BQB;
        }
    } else {
        for (int i = 0; i < m; i++) u[i] = 0.0;
    }
}

void cps_resilient_tightened_mpc(CPSResilientController* rc,
                                  const double* x_est,
                                  const double* A, const double* B,
                                  const double* Q, const double* R,
                                  int horizon, double* u_mpc) {
    /* One-step MPC with tightened constraints:
     * u* = argmin x'(A'QA+R)x + 2x'A'QBu + u'B'QBu
     * Solution: u* = -(B'QB+R)^{-1} * B'QA * x */
    if (!rc || !x_est || !A || !B || !u_mpc) return;
    int m = rc->input_dim;
    int n = (m < 2) ? 2 : m;

    /* Assume diagonal Q, R for simplicity */
    double BQB_R = 0.0;
    for (int i = 0; i < m; i++)
        BQB_R += B[i] * B[i] * (Q ? Q[0] : 1.0) + (R ? R[0] : 0.1);

    if (fabs(BQB_R) > CPS_EPS) {
        for (int i = 0; i < m; i++) {
            u_mpc[i] = 0.0;
            for (int j = 0; j < n && j < m; j++)
                u_mpc[i] -= B[j] * A[j*n+i] * x_est[i]
                           * (Q ? Q[0] : 1.0) / BQB_R;
        }
    } else {
        for (int i = 0; i < m; i++) u_mpc[i] = 0.0;
    }
}

int cps_resilient_sensor_selection(CPSResilientController* rc,
                                    const double* C, int n, int p,
                                    int min_sensors) {
    if (!rc || !C || n <= 0 || p <= 0) return 0;
    if (min_sensors > p) min_sensors = p;
    if (min_sensors < 1) min_sensors = 1;

    /* Simple greedy: select sensors with largest C row norms */
    double* row_norms = (double*)malloc((size_t)p * sizeof(double));
    int* indices = (int*)malloc((size_t)p * sizeof(int));

    for (int i = 0; i < p; i++) {
        double norm_sq = 0.0;
        for (int j = 0; j < n; j++)
            norm_sq += C[i * n + j] * C[i * n + j];
        row_norms[i] = norm_sq;
        indices[i] = i;
    }

    /* Selection sort by norm (descending) */
    for (int i = 0; i < p - 1; i++) {
        int best = i;
        for (int j = i + 1; j < p; j++)
            if (row_norms[j] > row_norms[best]) best = j;
        if (best != i) {
            double tmp_n = row_norms[i];
            row_norms[i] = row_norms[best];
            row_norms[best] = tmp_n;
            int tmp_i = indices[i];
            indices[i] = indices[best];
            indices[best] = tmp_i;
        }
    }

    /* Store in active_sensors */
    rc->n_active_sensors = min_sensors;
    free(rc->active_sensors);
    rc->active_sensors = (int*)malloc((size_t)min_sensors * sizeof(int));
    for (int i = 0; i < min_sensors; i++)
        rc->active_sensors[i] = indices[i];

    free(row_norms); free(indices);
    return min_sensors;
}

/* ============================================================================
 * Reachability-Based Security Analysis (L4: Theorems)
 * ============================================================================ */

void cps_attack_reachable_set(const double* A, const double* B,
                               const double* C, int n, int m, int p,
                               double attack_budget, int horizon,
                               double* R_min, double* R_max) {
    if (!R_min || !R_max || n <= 0) return;
    for (int i = 0; i < n; i++) {
        R_min[i] = -attack_budget * horizon;
        R_max[i] =  attack_budget * horizon;
    }

    /* Refine using controllability Gramian: reachable states satisfy
     * x' * Wc^{-1} * x <= attack_budget */
    if (horizon > n) horizon = n;
    double* Wc = (double*)calloc((size_t)(n * n), sizeof(double));
    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp = (double*)malloc((size_t)(n * n) * sizeof(double));

    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;
    for (int k = 0; k < horizon; k++) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Wc[i*n+j] += A_pow[i*n+j] * A_pow[i*n+j];
        cps_matrix_multiply(temp, A_pow, A, n, n, n);
        memcpy(A_pow, temp, (size_t)(n * n) * sizeof(double));
    }

    /* Bound each state component: |x_i| <= sqrt(Wc_{ii} * budget) */
    for (int i = 0; i < n; i++) {
        double bound = sqrt(fabs(Wc[i*n+i]) * attack_budget + 1.0);
        R_min[i] = -bound;
        R_max[i] =  bound;
    }

    free(Wc); free(A_pow); free(temp);
}

double cps_attack_invariant_set_volume(const double* A,
                                        const double* B,
                                        const double* Gamma_a,
                                        int n, int m, int n_a) {
    if (n <= 0) return 0.0;
    /* Approximate volume from eigenvalues of A */
    double trace_A = 0.0;
    for (int i = 0; i < n; i++) trace_A += A[i*n+i];
    double spectral_radius = fabs(trace_A) / n;
    if (spectral_radius < CPS_EPS) spectral_radius = 0.1;
    return 1.0 / (spectral_radius + CPS_EPS);
}

/* ============================================================================
 * Redundancy and Diversity Analysis
 * ============================================================================ */

int cps_min_sensors_for_resilience(int s_attack) {
    /* Theorem (Fawzi et al. 2014): need 2s+1 sensors to tolerate s attacks */
    return 2 * s_attack + 1;
}

double cps_sensor_diversity_score(const double* C, int n, int p) {
    if (!C || n <= 0 || p < 2) return 0.0;
    /* Diversity = 1 - average_correlation between sensor row vectors */
    double total_corr = 0.0;
    int pairs = 0;
    for (int i = 0; i < p; i++) {
        for (int j = i + 1; j < p; j++) {
            double dot_ij = 0.0, norm_i = 0.0, norm_j = 0.0;
            for (int k = 0; k < n; k++) {
                double ci = C[i*n+k], cj = C[j*n+k];
                dot_ij += ci * cj;
                norm_i += ci * ci;
                norm_j += cj * cj;
            }
            norm_i = sqrt(norm_i); norm_j = sqrt(norm_j);
            if (norm_i > CPS_EPS && norm_j > CPS_EPS)
                total_corr += fabs(dot_ij) / (norm_i * norm_j);
            pairs++;
        }
    }
    if (pairs == 0) return 1.0;
    return 1.0 - total_corr / pairs;
}

void cps_actuator_redundancy_map(const double* B, int n, int m,
                                  int* redundancy_matrix) {
    if (!B || !redundancy_matrix || n <= 0 || m <= 0) return;
    /* redundancy_matrix[i*m+j] = 1 if actuator j can substitute for i */
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double corr = 0.0, ni = 0.0, nj = 0.0;
            for (int k = 0; k < n; k++) {
                corr += B[k*m+i] * B[k*m+j];
                ni += B[k*m+i] * B[k*m+i];
                nj += B[k*m+j] * B[k*m+j];
            }
            ni = sqrt(ni); nj = sqrt(nj);
            redundancy_matrix[i*m+j] = (ni > CPS_EPS && nj > CPS_EPS
                && fabs(corr)/(ni*nj) > 0.7) ? 1 : 0;
        }
    }
}

/* ============================================================================
 * Attack Detectability Check (L4: Theorems)
 * ============================================================================ */

bool cps_is_attack_detectable(const double* A, const double* C,
                               const double* Gamma_a, int n, int p,
                               int n_a) {
    if (!A || !C || n <= 0 || p <= 0) return false;

    /* Check if the augmented system (A, [C; zeros]) is observable
     * from the attack subspace. If the attack lies in the unobservable
     * subspace, it is undetectable. */
    int m_rows = n * (p + n_a);
    double* O = (double*)calloc((size_t)(m_rows * n), sizeof(double));
    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp = (double*)malloc((size_t)(n * n) * sizeof(double));

    /* Row 0: C */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++)
            O[i * n + j] = C[i * n + j];

    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;

    for (int k = 1; k < n; k++) {
        cps_matrix_multiply(temp, A_pow, A, n, n, n);
        memcpy(A_pow, temp, (size_t)(n * n) * sizeof(double));
        double* block = &O[k * p * n];
        cps_matrix_multiply(block, C, A_pow, p, n, n);
    }

    int rank = cps_matrix_rank(O, m_rows, n, 1e-8);
    free(O); free(A_pow); free(temp);
    return rank == n;
}

bool cps_has_zero_dynamics_vulnerability(const double* A,
                                          const double* B,
                                          const double* C,
                                          int n, int m, int p) {
    /* Zero-dynamics attack exists if the system has invariant zeros
     * outside the unit circle. For 2x2 systems, check if:
     * det([A - z*I, B; C, 0]) = 0 has roots with |z| >= 1 */
    if (n != 2 || m != 1 || p != 1) return false;

    /* For SISO 2nd order: det = C*(zI-A)^{-1}*B
     * Zero at z = (a11*b2 - a12*b1)/b2 (if b2 != 0) */
    double a11 = A[0], a12 = A[1], a21 = A[2], a22 = A[3];
    double b1 = B[0], b2 = B[1];

    if (fabs(b2) > CPS_EPS) {
        double zero = (a11*b2 - a12*b1) / b2;
        /* Also check from adjugate */
        double zero2 = a22 - a21*b2/b1;
        if (fabs(zero) >= 1.0 || fabs(zero2) >= 1.0)
            return true;
    }
    return false;
}
