#include "packet_loss_controller.h"
#include "packet_loss_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Internal Matrix Helpers
 * ============================================================================ */

static double* _mat_alloc(int rows, int cols) {
    return (double*)calloc(rows * cols, sizeof(double));
}

/* ============================================================================
 * Matrix Multiplication: C = A (m×k) · B (k×n)
 * Row-major storage. Complexity: O(m·k·n).
 * ============================================================================ */

void pl_mat_mul(const double* A, const double* B, double* C,
                int m, int k, int n) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < k; l++)
                sum += A[i * k + l] * B[l * n + j];
            C[i * n + j] = sum;
        }
}

/* ============================================================================
 * Matrix Transpose: AT[j·rows+i] = A[i·cols+j]
 * Complexity: O(rows·cols).
 * ============================================================================ */

void pl_mat_transpose(const double* A, double* AT, int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            AT[j * rows + i] = A[i * cols + j];
}

/* ============================================================================
 * Gaussian Elimination with Partial Pivoting for Matrix Inversion
 *
 * Forms augmented matrix [A | I], reduces to [I | A^{-1}].
 * Partial pivoting ensures numerical stability.
 *
 * Complexity: O(n³). Returns true on success, false if singular.
 * Reference: Golub & Van Loan, "Matrix Computations" (2013), §3.2.
 * ============================================================================ */

bool pl_mat_invert(double* A, int n) {
    if (n <= 0 || !A) return false;
    int tc = 2 * n;
    double* aug = (double*)malloc(n * tc * sizeof(double));
    if (!aug) return false;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i * tc + j] = A[i * n + j];
        for (int j = n; j < tc; j++)
            aug[i * tc + j] = (i == (j - n)) ? 1.0 : 0.0;
    }

    for (int col = 0; col < n; col++) {
        int prow = col;
        double mv = fabs(aug[col * tc + col]);
        for (int r = col + 1; r < n; r++) {
            double v = fabs(aug[r * tc + col]);
            if (v > mv) { mv = v; prow = r; }
        }
        if (mv < 1e-14) { free(aug); return false; }

        if (prow != col)
            for (int j = 0; j < tc; j++) {
                double t = aug[col * tc + j];
                aug[col * tc + j] = aug[prow * tc + j];
                aug[prow * tc + j] = t;
            }

        double piv = aug[col * tc + col];
        for (int j = 0; j < tc; j++) aug[col * tc + j] /= piv;

        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double fac = aug[r * tc + col];
            for (int j = 0; j < tc; j++)
                aug[r * tc + j] -= fac * aug[col * tc + j];
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i * n + j] = aug[i * tc + n + j];
    free(aug);
    return true;
}

/* ============================================================================
 * Cholesky Decomposition: A = L·L' for Symmetric Positive Definite A
 *
 * Algorithm: for j=0..n-1, L_jj = sqrt(A_jj - Σ_{k<j} L_jk²)
 *                       L_ij = (A_ij - Σ_{k<j} L_ik·L_jk) / L_jj
 *
 * Lower triangle of A is overwritten with L.
 * Complexity: O(n³/3). No pivoting needed for SPD matrices.
 *
 * Reference: Golub & Van Loan (2013), §4.2.
 * ============================================================================ */

bool pl_mat_cholesky(double* A, int n) {
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < j; k++)
            sum += A[j * n + k] * A[j * n + k];
        double diag = A[j * n + j] - sum;
        if (diag <= 0.0) return false;
        A[j * n + j] = sqrt(diag);
        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++)
                sum += A[i * n + k] * A[j * n + k];
            A[i * n + j] = (A[i * n + j] - sum) / A[j * n + j];
        }
    }
    return true;
}

/* ============================================================================
 * Spectral Radius via Power Iteration (von Mises)
 *
 * ρ(A) = max_i |λ_i(A)| — the magnitude of the dominant eigenvalue.
 * Rayleigh quotient: λ = (v'·A·v) / (v'·v).
 *
 * Convergence rate: |λ_2/λ_1|. Fails for complex dominant eigenvalues
 * of equal magnitude; in that case returns an approximation.
 *
 * Used as the key parameter in critical loss probability formulas.
 * ============================================================================ */

double pl_spectral_radius_power(const double* A, int n,
                                 int max_iter, double tol) {
    if (n <= 0 || !A) return 0.0;
    double* v  = (double*)malloc(n * sizeof(double));
    double* Av = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) v[i] = 1.0;

    double lambda = 0.0, lambda_old = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++) Av[i] += A[i * n + j] * v[j];
        }
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; i++) { num += v[i] * Av[i]; den += v[i] * v[i]; }
        lambda = (den > 1e-15) ? num / den : 0.0;
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += Av[i] * Av[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (int i = 0; i < n; i++) v[i] = Av[i] / norm;
        if (fabs(lambda - lambda_old) < tol) break;
        lambda_old = lambda;
    }
    free(v); free(Av);
    return fabs(lambda);
}

/* ============================================================================
 * LTI System: Lifecycle, Configuration, Simulation
 * ============================================================================ */

LTISystem* pl_lti_create(int n, int m, int p, int q) {
    LTISystem* s = (LTISystem*)calloc(1, sizeof(LTISystem));
    if (!s) return NULL;
    s->n = n; s->m = m; s->p = p;
    s->q = (q > 0) ? q : n;
    s->A  = _mat_alloc(n, n);
    s->B  = _mat_alloc(n, m);
    s->C  = _mat_alloc(p, n);
    s->Cz = _mat_alloc(s->q, n);
    s->Q  = _mat_alloc(n, n);
    s->R  = _mat_alloc(p, p);
    for (int i = 0; i < n; i++) s->A[i * n + i] = 1.0;
    for (int i = 0; i < n; i++) s->Q[i * n + i] = 0.01;
    for (int i = 0; i < p; i++) s->R[i * p + i] = 0.01;
    return s;
}

void pl_lti_free(LTISystem* s) {
    if (!s) return;
    free(s->A); free(s->B); free(s->C); free(s->Cz);
    free(s->Q); free(s->R); free(s);
}

double pl_lti_spectral_radius(const LTISystem* s) {
    return pl_spectral_radius_power(s->A, s->n, 1000, 1e-10);
}

void pl_lti_set_A(LTISystem* s, const double* d, int n) {
    if (s && d && n == s->n) memcpy(s->A, d, n * n * sizeof(double));
}

void pl_lti_set_B(LTISystem* s, const double* d) {
    if (s && d) memcpy(s->B, d, s->n * s->m * sizeof(double));
}

void pl_lti_set_C(LTISystem* s, const double* d) {
    if (s && d) memcpy(s->C, d, s->p * s->n * sizeof(double));
}

void pl_lti_set_Q(LTISystem* s, const double* d) {
    if (s && d) memcpy(s->Q, d, s->n * s->n * sizeof(double));
}

void pl_lti_set_R(LTISystem* s, const double* d) {
    if (s && d) memcpy(s->R, d, s->p * s->p * sizeof(double));
}

/* ============================================================================
 * LTI System Step and Measurement Simulation
 * ============================================================================ */

void pl_lti_step(const LTISystem* s, const double* x, const double* u,
                 bool add_noise, unsigned long* rng, double* xn) {
    int n = s->n, m = s->m;
    for (int i = 0; i < n; i++) {
        xn[i] = 0.0;
        for (int j = 0; j < n; j++) xn[i] += s->A[i * n + j] * x[j];
        for (int j = 0; j < m; j++) xn[i] += s->B[i * m + j] * u[j];
    }
    if (add_noise && rng) {
        for (int i = 0; i < n; i++) {
            double sd = sqrt(s->Q[i * n + i]);
            double u1 = pl_uniform(rng), u2 = pl_uniform(rng);
            if (u1 < 1e-15) u1 = 1e-15;
            xn[i] += sd * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        }
    }
}

void pl_lti_measure(const LTISystem* s, const double* x, bool add_noise,
                    unsigned long* rng, double* y) {
    int n = s->n, p = s->p;
    for (int i = 0; i < p; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) y[i] += s->C[i * n + j] * x[j];
    }
    if (add_noise && rng) {
        for (int i = 0; i < p; i++) {
            double sd = sqrt(s->R[i * p + i]);
            double u1 = pl_uniform(rng), u2 = pl_uniform(rng);
            if (u1 < 1e-15) u1 = 1e-15;
            y[i] += sd * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        }
    }
}

/* ============================================================================
 * Controllability & Observability (Kalman Rank Conditions)
 *
 * Controllability: C = [B, AB, A²B, ..., A^{n-1}B], rank must be n.
 * Observability:   O = [C; CA; ...; CA^{n-1}], rank must be n.
 *
 * Rank computed via Gaussian elimination (row reduction with pivoting).
 * ============================================================================ */

double* pl_lti_controllability_matrix(const LTISystem* s, int* rank_out) {
    int n = s->n, m = s->m, nm = n * m;
    double* Cc = _mat_alloc(n, nm);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Cc[i * nm + j] = s->B[i * m + j];

    double* Ap = _mat_alloc(n, n), *tmp = _mat_alloc(n, n);
    memcpy(Ap, s->A, n * n * sizeof(double));
    for (int k = 1; k < n; k++) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                double v = 0.0;
                for (int l = 0; l < n; l++) v += Ap[i * n + l] * s->B[l * m + j];
                Cc[i * nm + k * m + j] = v;
            }
        pl_mat_mul(Ap, s->A, tmp, n, n, n);
        memcpy(Ap, tmp, n * n * sizeof(double));
    }

    double* cp = (double*)malloc(n * nm * sizeof(double));
    memcpy(cp, Cc, n * nm * sizeof(double));
    int rk = 0;
    for (int col = 0; col < nm && rk < n; col++)
        for (int row = rk; row < n; row++)
            if (fabs(cp[row * nm + col]) > 1e-10) {
                for (int j = 0; j < nm; j++) {
                    double t = cp[rk * nm + j];
                    cp[rk * nm + j] = cp[row * nm + j]; cp[row * nm + j] = t;
                }
                double piv = cp[rk * nm + col];
                for (int r = rk + 1; r < n; r++) {
                    double fac = cp[r * nm + col] / piv;
                    if (fabs(fac) > 1e-15)
                        for (int j = 0; j < nm; j++)
                            cp[r * nm + j] -= fac * cp[rk * nm + j];
                }
                rk++; break;
            }
    *rank_out = rk;
    free(cp); free(Ap); free(tmp);
    return Cc;
}

double* pl_lti_observability_matrix(const LTISystem* s, int* rank_out) {
    int n = s->n, p = s->p, np = n * p;
    double* Ob = _mat_alloc(np, n);
    memcpy(Ob, s->C, p * n * sizeof(double));
    double* Ap = _mat_alloc(n, n), *tmp = _mat_alloc(n, n);
    memcpy(Ap, s->A, n * n * sizeof(double));
    for (int k = 1; k < n; k++) {
        for (int i = 0; i < p; i++)
            for (int j = 0; j < n; j++) {
                double v = 0.0;
                for (int l = 0; l < n; l++) v += s->C[i * n + l] * Ap[l * n + j];
                Ob[(k * p + i) * n + j] = v;
            }
        pl_mat_mul(Ap, s->A, tmp, n, n, n);
        memcpy(Ap, tmp, n * n * sizeof(double));
    }
    double* cp = (double*)malloc(np * n * sizeof(double));
    memcpy(cp, Ob, np * n * sizeof(double));
    int rk = 0;
    for (int col = 0; col < n && rk < np; col++)
        for (int row = rk; row < np; row++)
            if (fabs(cp[row * n + col]) > 1e-10) {
                for (int j = 0; j < n; j++) {
                    double t = cp[rk * n + j];
                    cp[rk * n + j] = cp[row * n + j]; cp[row * n + j] = t;
                }
                double piv = cp[rk * n + col];
                for (int r = rk + 1; r < np; r++) {
                    double fac = cp[r * n + col] / piv;
                    if (fabs(fac) > 1e-15)
                        for (int j = 0; j < n; j++)
                            cp[r * n + j] -= fac * cp[rk * n + j];
                }
                rk++; break;
            }
    *rank_out = rk;
    free(cp); free(Ap); free(tmp);
    return Ob;
}

bool pl_lti_is_stabilizable(const LTISystem* s) {
    int rk; double* Cc = pl_lti_controllability_matrix(s, &rk); free(Cc);
    if (rk == s->n) return true;
    return pl_lti_spectral_radius(s) < 1.0;
}

bool pl_lti_is_detectable(const LTISystem* s) {
    int rk; double* Ob = pl_lti_observability_matrix(s, &rk); free(Ob);
    if (rk == s->n) return true;
    return pl_lti_spectral_radius(s) < 1.0;
}

/* ============================================================================
 * DARE Solver — Discrete Algebraic Riccati Equation
 *
 * P = A'PA - A'PB(R + B'PB)^{-1}B'PA + Q
 *
 * Solved by value iteration (Bertsekas 2012, §4.1):
 *   P_{k+1} = A'P_k A - A'P_k B(R+B'P_k B)^{-1} B'P_k A + Q,  P_0 = I.
 *
 * Optimal LQR gain: L = (R + B'PB)^{-1} B'PA.
 * Convergence guaranteed if (A,B) stabilizable, (A,sqrt(Q)) detectable.
 * Rate: linear with ρ(A-BL)².
 * ============================================================================ */

LQRSolution* pl_lqr_solve(const LTISystem* s, double* Rc,
                           int max_iter, double tol) {
    int n = s->n, m = s->m;
    LQRSolution* sol = (LQRSolution*)calloc(1, sizeof(LQRSolution));
    sol->n = n; sol->m = m; sol->tolerance = tol;
    sol->P = _mat_alloc(n, n);
    sol->L = _mat_alloc(m, n);
    sol->K_cl = _mat_alloc(n, n);
    for (int i = 0; i < n; i++) sol->P[i * n + i] = 1.0;

    double *Pn = _mat_alloc(n, n), *AT = _mat_alloc(n, n);
    double *T1 = _mat_alloc(n, n), *T2 = _mat_alloc(n, m), *T3 = _mat_alloc(m, m);
    double *Tn = _mat_alloc(n, m), *BT = _mat_alloc(m, n);
    pl_mat_transpose(s->A, AT, n, n);
    pl_mat_transpose(s->B, BT, n, m);

    for (int iter = 0; iter < max_iter; iter++) {
        double *PA = _mat_alloc(n, n);
        pl_mat_mul(s->A, sol->P, PA, n, n, n);
        pl_mat_mul(AT, PA, T1, n, n, n);
        free(PA);

        double *PB = _mat_alloc(n, m);
        pl_mat_mul(sol->P, s->B, PB, n, n, m);
        pl_mat_mul(AT, PB, T2, n, n, m);
        free(PB);

        pl_mat_mul(sol->P, s->B, Tn, n, n, m);
        pl_mat_mul(BT, Tn, T3, m, n, m);
        for (int i = 0; i < m * m; i++) T3[i] += Rc[i];

        double *T3i = (double*)malloc(m * m * sizeof(double));
        memcpy(T3i, T3, m * m * sizeof(double));
        if (!pl_mat_invert(T3i, m)) { free(T3i); sol->iterations = iter; break; }

        double *T2T = _mat_alloc(m, n);
        pl_mat_transpose(T2, T2T, n, m);
        double *Tm = _mat_alloc(n, m);
        pl_mat_mul(T2, T3i, Tm, n, m, m);
        double *Tf = _mat_alloc(n, n);
        pl_mat_mul(Tm, T2T, Tf, n, m, n);
        for (int i = 0; i < n * n; i++) Pn[i] = T1[i] - Tf[i] + s->Q[i];

        double diff = 0.0;
        for (int i = 0; i < n * n; i++) {
            double d = Pn[i] - sol->P[i]; diff += d * d;
        }
        memcpy(sol->P, Pn, n * n * sizeof(double));
        sol->iterations = iter + 1;
        free(T3i); free(T2T); free(Tm); free(Tf);
        if (sqrt(diff) < tol) break;
    }

    double *PA2 = _mat_alloc(n, n), *BTPA = _mat_alloc(m, n);
    pl_mat_mul(sol->P, s->A, PA2, n, n, n);
    pl_mat_mul(BT, PA2, BTPA, m, n, n);
    double *Rp = _mat_alloc(m, m);
    pl_mat_mul(sol->P, s->B, Tn, n, n, m);
    pl_mat_mul(BT, Tn, Rp, m, n, m);
    for (int i = 0; i < m * m; i++) Rp[i] += Rc[i];
    if (pl_mat_invert(Rp, m)) pl_mat_mul(Rp, BTPA, sol->L, m, m, n);

    double *BL = _mat_alloc(n, n);
    pl_mat_mul(s->B, sol->L, BL, n, m, n);
    for (int i = 0; i < n * n; i++) sol->K_cl[i] = s->A[i] - BL[i];

    free(Pn); free(AT); free(T1); free(T2); free(T3);
    free(Tn); free(BT); free(PA2); free(BTPA); free(Rp); free(BL);
    return sol;
}

void pl_lqr_free(LQRSolution* sol) {
    if (sol) { free(sol->P); free(sol->L); free(sol->K_cl); free(sol); }
}

double pl_lqr_spectral_radius(const LQRSolution* sol, const LTISystem* s) {
    (void)s;
    return pl_spectral_radius_power(sol->K_cl, sol->n, 1000, 1e-10);
}

double pl_lqr_cost(const LTISystem* s, const LQRSolution* sol,
                   const double* xt, const double* ut, int T) {
    double c = 0.0;
    for (int t = 0; t < T; t++) {
        for (int i = 0; i < s->n; i++)
            for (int j = 0; j < s->n; j++)
                c += xt[t * s->n + i] * s->Q[i * s->n + j] * xt[t * s->n + j];
        for (int i = 0; i < s->m; i++)
            for (int j = 0; j < s->m; j++)
                c += ut[t * s->m + i] * 0.1 * ut[t * s->m + j];
    }
    for (int i = 0; i < s->n; i++)
        for (int j = 0; j < s->n; j++)
            c += xt[(T-1)*s->n + i] * sol->P[i * s->n + j] * xt[(T-1)*s->n + j];
    return c;
}

/* ============================================================================
 * TCP-like Controller
 *
 * Under TCP-like protocols with ACKs, the controller knows which control
 * signals were applied. Separation principle holds (Schenato et al., 2007):
 *
 *   u_k = -L · x̂_{k|k}
 *
 * If actuator ACK = loss, apply hold strategy:
 *   ZERO_INPUT → u = 0            (safety)
 *   ZERO_ORDER → u = u_{last}     (simple)
 *   PREDICTIVE → u = -L·(A-BL)·x̂ (one-step prediction)
 * ============================================================================ */

TCPController* pl_tcp_controller_create(const LTISystem* s, const LQRSolution* lqr) {
    int n = s->n, m = s->m;
    TCPController* c = (TCPController*)calloc(1, sizeof(TCPController));
    if (!c) return NULL;
    c->sys = s;
    c->L = _mat_alloc(m, n); memcpy(c->L, lqr->L, m * n * sizeof(double));
    c->x_hat = (double*)calloc(n, sizeof(double));
    c->P = _mat_alloc(n, n); for (int i = 0; i < n; i++) c->P[i * n + i] = 1.0;
    c->x_pred = (double*)calloc(n, sizeof(double));
    c->P_pred = _mat_alloc(n, n); for (int i = 0; i < n; i++) c->P_pred[i * n + i] = 1.0;
    c->u_last = (double*)calloc(m, sizeof(double));
    c->u_sent = (double*)calloc(m, sizeof(double));
    c->protocol = PROTO_TCP_LIKE;
    c->sensor_hold = HOLD_ZERO_ORDER;
    c->actuator_hold = HOLD_ZERO_ORDER;
    c->cost_history_cap = 1000;
    c->cost_history = (double*)malloc(1000 * sizeof(double));
    return c;
}

void pl_tcp_controller_free(TCPController* c) {
    if (!c) return;
    free(c->L); free(c->x_hat); free(c->P); free(c->x_pred);
    free(c->P_pred); free(c->u_last); free(c->u_sent);
    free(c->cost_history); free(c);
}

void pl_tcp_controller_set_hold(TCPController* c, HoldStrategy sh, HoldStrategy ah) {
    if (c) { c->sensor_hold = sh; c->actuator_hold = ah; }
}

const double* pl_tcp_controller_update_estimate(TCPController* c, const double* y,
                                                  PacketStatus ss) {
    if (!c) return NULL;
    int n = c->sys->n, p = c->sys->p;
    pl_mat_mul(c->sys->A, c->x_hat, c->x_pred, n, n, 1);

    if (ss == PACKET_RECEIVED && y) {
        double *AP = _mat_alloc(n, n), *AT = _mat_alloc(n, n);
        pl_mat_mul(c->sys->A, c->P, AP, n, n, n);
        pl_mat_transpose(c->sys->A, AT, n, n);
        pl_mat_mul(AP, AT, c->P_pred, n, n, n);
        for (int i = 0; i < n * n; i++) c->P_pred[i] += c->sys->Q[i];
        free(AP); free(AT);

        double *CT = _mat_alloc(n, p), *PCT = _mat_alloc(n, p);
        pl_mat_transpose(c->sys->C, CT, p, n);
        pl_mat_mul(c->P_pred, CT, PCT, n, n, p);
        double *CPC = _mat_alloc(p, p), *CP = _mat_alloc(p, n);
        pl_mat_mul(c->sys->C, c->P_pred, CP, p, n, n);
        pl_mat_mul(CP, CT, CPC, p, n, p);
        for (int i = 0; i < p * p; i++) CPC[i] += c->sys->R[i];
        free(CP); free(CT);

        double *K = _mat_alloc(n, p);
        if (pl_mat_invert(CPC, p)) pl_mat_mul(PCT, CPC, K, n, p, p);

        double *inv = (double*)malloc(p * sizeof(double));
        for (int i = 0; i < p; i++) {
            double cp = 0.0;
            for (int j = 0; j < n; j++) cp += c->sys->C[i * n + j] * c->x_pred[j];
            inv[i] = y[i] - cp;
        }
        for (int i = 0; i < n; i++) {
            c->x_hat[i] = c->x_pred[i];
            for (int j = 0; j < p; j++) c->x_hat[i] += K[i * p + j] * inv[j];
        }
        double *KC = _mat_alloc(n, n), *ImKC = _mat_alloc(n, n);
        pl_mat_mul(K, c->sys->C, KC, n, p, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                ImKC[i * n + j] = (i == j) ? 1.0 - KC[i * n + j] : -KC[i * n + j];
        pl_mat_mul(ImKC, c->P_pred, c->P, n, n, n);
        free(K); free(inv); free(PCT); free(CPC); free(KC); free(ImKC);
    } else {
        switch (c->sensor_hold) {
            case HOLD_ZERO_INPUT: memset(c->x_hat, 0, n * sizeof(double)); break;
            default: memcpy(c->x_hat, c->x_pred, n * sizeof(double)); break;
        }
        c->sensor_losses++;
    }
    c->total_sensor_attempts++;
    return c->x_hat;
}

const double* pl_tcp_controller_compute_control(TCPController* c,
                                                  PacketStatus as) {
    if (!c) return NULL;
    int n = c->sys->n, m = c->sys->m;
    for (int i = 0; i < m; i++) {
        c->u_sent[i] = 0.0;
        for (int j = 0; j < n; j++) c->u_sent[i] -= c->L[i * n + j] * c->x_hat[j];
    }
    if (as == PACKET_RECEIVED) {
        memcpy(c->u_last, c->u_sent, m * sizeof(double));
        return c->u_last;
    }
    c->actuator_losses++;
    switch (c->actuator_hold) {
        case HOLD_ZERO_INPUT:
            memset(c->u_last, 0, m * sizeof(double)); break;
        case HOLD_ZERO_ORDER: break;
        case HOLD_PREDICTIVE: {
            double *xn = (double*)malloc(n * sizeof(double));
            double *Acl = _mat_alloc(n, n), *BL = _mat_alloc(n, n);
            pl_mat_mul(c->sys->B, c->L, BL, n, m, n);
            for (int i = 0; i < n * n; i++) Acl[i] = c->sys->A[i] - BL[i];
            pl_mat_mul(Acl, c->x_hat, xn, n, n, 1);
            for (int i = 0; i < m; i++) {
                c->u_last[i] = 0.0;
                for (int j = 0; j < n; j++) c->u_last[i] -= c->L[i*n+j] * xn[j];
            }
            free(xn); free(Acl); free(BL);
        } break;
        default: break;
    }
    c->total_actuator_attempts++;
    return c->u_last;
}

void pl_tcp_controller_step(TCPController* c, const double* y,
                             PacketStatus ss, PacketStatus as, double* uo) {
    pl_tcp_controller_update_estimate(c, y, ss);
    const double* u = pl_tcp_controller_compute_control(c, as);
    if (uo && u) memcpy(uo, u, c->sys->m * sizeof(double));
}

void pl_tcp_controller_reset(TCPController* c) {
    if (!c) return;
    memset(c->x_hat, 0, c->sys->n * sizeof(double));
    memset(c->u_last, 0, c->sys->m * sizeof(double));
    c->sensor_losses = 0; c->actuator_losses = 0;
    c->total_sensor_attempts = 0; c->total_actuator_attempts = 0;
}

void pl_tcp_controller_print(const TCPController* c) {
    if (!c) { printf("TCPController: NULL\n"); return; }
    printf("=== TCP Controller ===\n");
    printf("Sensor loss: %lu/%lu (%.1f%%)  Actuator loss: %lu/%lu (%.1f%%)\n",
           c->sensor_losses, c->total_sensor_attempts,
           100.0*c->sensor_losses/(c->total_sensor_attempts+1),
           c->actuator_losses, c->total_actuator_attempts,
           100.0*c->actuator_losses/(c->total_actuator_attempts+1));
}

/* ============================================================================
 * UDP-like Controller
 *
 * Without ACKs, separation fails (Gupta et al., 2007).
 * Practical approach: pessimistic LQR — scale gains by γ = 1-p.
 * ============================================================================ */

UDPController* pl_udp_controller_create(const LTISystem* s, const LQRSolution* lqr) {
    int n = s->n, m = s->m;
    UDPController* c = (UDPController*)calloc(1, sizeof(UDPController));
    if (!c) return NULL;
    c->sys = s;
    c->L = _mat_alloc(m, n); memcpy(c->L, lqr->L, m * n * sizeof(double));
    c->x_hat = (double*)calloc(n, sizeof(double));
    c->P = _mat_alloc(n, n); for (int i = 0; i < n; i++) c->P[i*n+i] = 1.0;
    c->u_last_sent = (double*)calloc(m, sizeof(double));
    c->prob_received = 1.0;
    c->cost_history_cap = 1000;
    c->cost_history = (double*)malloc(1000 * sizeof(double));
    return c;
}

void pl_udp_controller_free(UDPController* c) {
    if (!c) return;
    free(c->L); free(c->x_hat); free(c->P); free(c->u_last_sent);
    free(c->cost_history); free(c);
}

const double* pl_udp_controller_update_estimate(UDPController* c,
                                                  const double* y, PacketStatus ss) {
    if (!c) return NULL;
    int n = c->sys->n, p = c->sys->p;
    double *xp = (double*)malloc(n * sizeof(double));
    pl_mat_mul(c->sys->A, c->x_hat, xp, n, n, 1);
    if (ss == PACKET_RECEIVED && y) {
        double *inv = (double*)malloc(p * sizeof(double));
        for (int i = 0; i < p; i++) {
            double cp = 0.0;
            for (int j = 0; j < n; j++) cp += c->sys->C[i*n+j] * xp[j];
            inv[i] = y[i] - cp;
        }
        for (int i = 0; i < n; i++) c->x_hat[i] = xp[i] + 0.3 * inv[0];
        free(inv);
    } else {
        memcpy(c->x_hat, xp, n * sizeof(double));
    }
    free(xp);
    return c->x_hat;
}

const double* pl_udp_controller_compute_control(UDPController* c, double lp) {
    if (!c) return NULL;
    double g = 1.0 - lp; if (g < 0.0) g = 0.0; if (g > 1.0) g = 1.0;
    for (int i = 0; i < c->sys->m; i++) {
        c->u_last_sent[i] = 0.0;
        for (int j = 0; j < c->sys->n; j++)
            c->u_last_sent[i] -= g * c->L[i * c->sys->n + j] * c->x_hat[j];
    }
    return c->u_last_sent;
}

void pl_udp_controller_step(UDPController* c, const double* y,
                             PacketStatus ss, double alp, double* uo) {
    pl_udp_controller_update_estimate(c, y, ss);
    const double* u = pl_udp_controller_compute_control(c, alp);
    if (uo && u) memcpy(uo, u, c->sys->m * sizeof(double));
}

void pl_udp_controller_reset(UDPController* c) {
    if (!c) return;
    memset(c->x_hat, 0, c->sys->n * sizeof(double));
    memset(c->u_last_sent, 0, c->sys->m * sizeof(double));
}

void pl_udp_controller_print(const UDPController* c) {
    if (!c) return;
    printf("=== UDP Controller ===  Prob(received)=%.4f\n", c->prob_received);
}

/* ============================================================================
 * Lyapunov Equation Solver: A·X + X·A' + Q = 0
 *
 * Solved by Kronecker linearization:
 *   (I ⊗ A + A ⊗ I) · vec(X) = -vec(Q)
 *
 * Forms an n²×n² system. Tractable for n ≤ 16.
 * Uses Gaussian elimination with partial pivoting.
 *
 * Reference: Bartels & Stewart (1972), "Solution of the Matrix
 * Equation AX + XB = C", Comm. ACM 15(9):820-826.
 * ============================================================================ */

bool pl_solve_lyapunov(const double* A, const double* Q, double* X, int n) {
    if (!A || !Q || !X || n <= 0 || n > 16) return false;
    int n2 = n * n;
    double* M = _mat_alloc(n2, n2);
    double* b = (double*)malloc(n2 * sizeof(double));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                for (int l = 0; l < n; l++) {
                    int row = i * n + k, col = j * n + l;
                    double val = 0.0;
                    if (i == j) val += A[k * n + l];
                    if (k == l) val += A[i * n + j];
                    M[row * n2 + col] = val;
                }

    for (int i = 0; i < n2; i++) b[i] = -Q[i];

    double* aug = (double*)malloc(n2 * (n2 + 1) * sizeof(double));
    for (int i = 0; i < n2; i++) {
        memcpy(aug + i * (n2 + 1), M + i * n2, n2 * sizeof(double));
        aug[i * (n2 + 1) + n2] = b[i];
    }

    bool ok = true;
    for (int col = 0; col < n2 && ok; col++) {
        int prow = col;
        double mv = fabs(aug[col * (n2 + 1) + col]);
        for (int r = col + 1; r < n2; r++) {
            double v = fabs(aug[r * (n2 + 1) + col]);
            if (v > mv) { mv = v; prow = r; }
        }
        if (mv < 1e-14) { ok = false; break; }
        if (prow != col)
            for (int j = 0; j <= n2; j++) {
                double t = aug[col * (n2 + 1) + j];
                aug[col * (n2 + 1) + j] = aug[prow * (n2 + 1) + j];
                aug[prow * (n2 + 1) + j] = t;
            }
        double piv = aug[col * (n2 + 1) + col];
        for (int j = col; j <= n2; j++) aug[col * (n2 + 1) + j] /= piv;
        for (int r = 0; r < n2; r++) {
            if (r == col) continue;
            double fac = aug[r * (n2 + 1) + col];
            for (int j = col; j <= n2; j++)
                aug[r * (n2 + 1) + j] -= fac * aug[col * (n2 + 1) + j];
        }
    }

    if (ok)
        for (int i = 0; i < n2; i++) X[i] = aug[i * (n2 + 1) + n2];

    free(M); free(b); free(aug);
    return ok;
}
