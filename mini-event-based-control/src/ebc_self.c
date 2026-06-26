#include "ebc_self.h"
#include "ebc_core.h"
#include "ebc_trigger.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*
 * ebc_self.c -- Self-Triggered Control (L5: Algorithms, L8: Advanced)
 *
 * STC computes the next control update time proactively based on
 * the current state and a plant model. This eliminates the need
 * for continuous state monitoring required by classic ETC.
 *
 * Key references:
 *   Mazo, Anta & Tabuada (2010) -- ISS self-triggered implementation
 *   Anta & Tabuada (2010) -- Self-triggered stabilization
 *   Higham (2005) -- Matrix exponential via scaling-and-squaring
 *   Van Loan (1978) -- Computing integrals of matrix exponential
 */

/* ---------- Internal: matrix multiplication (n x n) * (n x n) ---------- */
static void mat_mul(double* A, double* B, double* C, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++) s += A[i * n + k] * B[k * n + j];
            C[i * n + j] = s;
        }
}

/* ---------- Internal: matrix-vector multiply ---------- */
static void mat_vec(double* A, double* x, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) y[i] += A[i * n + j] * x[j];
    }
}

/* ---------- Internal: vector norm ---------- */
static double vec_norm(const double* v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}

/* ================================================================
 * Matrix exponential: exp(A*t) via truncated Taylor series
 * with scaling-and-squaring.
 *
 * exp(A*t) approximated as:
 *   [T_{m}(A*t / 2^s)]^{2^s}
 * where T_m is the Taylor polynomial of degree m.
 *
 * Reference: Higham (2005), SIAM J. Matrix Anal. Appl., 26(4).
 * Complexity: O(n^3 * (m + s))
 * ================================================================ */

void ebc_matrix_exponential(const double* A, int n, double t, double* E) {
    if (!A || n < 1 || !E) return;
    int i, j;

    /* Allocate and scale */
    double* B = malloc(n * n * sizeof(double));
    double* T = malloc(n * n * sizeof(double));
    double* W = malloc(n * n * sizeof(double));
    if (!B || !T || !W) { free(B); free(T); free(W); return; }

    for (i = 0; i < n * n; i++) B[i] = A[i] * t;

    /* Compute Frobenius norm for scaling */
    double nrm = 0.0;
    for (i = 0; i < n * n; i++) nrm += B[i] * B[i];
    nrm = sqrt(nrm);

    int s = (nrm > 1.0) ? (int)ceil(log(nrm) / log(2.0)) + 1 : 0;
    double scale = 1.0 / (double)(1 << s);
    for (i = 0; i < n * n; i++) B[i] *= scale;

    /* Taylor series: T = I + B + B^2/2! + ... + B^m/m! */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            T[i * n + j] = (i == j) ? 1.0 : 0.0;

    /* Start with B_power = B, accumulate into T */
    double* Bpow = malloc(n * n * sizeof(double));
    if (!Bpow) { free(B); free(T); free(W); return; }
    memcpy(Bpow, B, n * n * sizeof(double));

    for (i = 0; i < n * n; i++) T[i] += B[i];  /* T = I + B */

    double fact = 2.0;
    for (int p = 2; p <= 12; p++) {
        /* Bpow = Bpow * B */
        mat_mul(Bpow, B, W, n);
        memcpy(Bpow, W, n * n * sizeof(double));

        for (i = 0; i < n * n; i++) T[i] += Bpow[i] / fact;
        fact *= (double)(p + 1);
    }

    /* Squaring: T := T^{2^s} */
    memcpy(E, T, n * n * sizeof(double));
    for (int sq = 0; sq < s; sq++) {
        mat_mul(E, E, W, n);
        memcpy(E, W, n * n * sizeof(double));
    }

    free(B); free(T); free(W); free(Bpow);
}

/* ================================================================
 * Matrix exponential integral: G(tau) = integral_0^tau exp(A*s) ds
 *
 * Uses the augmented matrix method (Van Loan 1978):
 *   exp([[A, I]; [0, 0]] * tau)[1:n, n+1:2n] = G(tau)
 * ================================================================ */

void ebc_matrix_exp_integral(const double* A, int n, double tau, double* G) {
    if (!A || n < 1 || !G) return;
    int i, j, N = 2 * n;

    double* M = calloc(N * N, sizeof(double));
    double* expM = calloc(N * N, sizeof(double));
    if (!M || !expM) { free(M); free(expM); return; }

    /* Build augmented matrix M */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            M[i * N + j] = A[i * n + j];     /* top-left: A */
        M[i * N + (n + i)] = 1.0;            /* top-right: I */
    }
    /* Bottom rows are zero already */

    /* Compute exp(M * tau) */
    ebc_matrix_exponential(M, N, tau, expM);

    /* Extract G from top-right block */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            G[i * n + j] = expM[i * N + (n + j)];

    free(M); free(expM);
}

/* ================================================================
 * Next event time for linear STC via bisection search.
 *
 * Given x(t_k), find max tau such that:
 *   |x(t_k+tau) - x_k| <= sigma*|x_k| + epsilon
 *
 * x(t_k+tau) = exp(A*tau)*x_k + G(tau)*B*K*x_k
 *
 * Uses bisection in [0, tau_max] with O(log(tau_max/tol)) iterations.
 * ================================================================ */

double ebc_self_next_time_linear(const double* A, const double* B,
                                  const double* K,
                                  const double* x_k,
                                  int n, int m,
                                  double sigma, double epsilon,
                                  double tau_max, double tol) {
    if (!A || !x_k || n < 1 || tau_max <= 0.0) return 0.001;
    int i, j;

    double norm_xk = vec_norm(x_k, n);
    if (norm_xk < 1e-12) return tau_max;
    double thresh = sigma * norm_xk + epsilon;

    /* Precompute BKx = B * K * x_k */
    double* BKx = calloc(n, sizeof(double));
    if (BKx && B && K && m > 0) {
        double* Kx = calloc(m, sizeof(double));
        if (Kx) {
            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    Kx[i] += K[i * n + j] * x_k[j];
            for (i = 0; i < n; i++)
                for (j = 0; j < m; j++)
                    BKx[i] += B[i * m + j] * Kx[j];
            free(Kx);
        }
    }

    /* Bisection */
    double lo = 0.0, hi = tau_max;
    for (int iter = 0; iter < 60; iter++) {
        double mid = (lo + hi) / 2.0;
        if (hi - lo < tol) break;

        double* expA = malloc(n * n * sizeof(double));
        double* Gi = malloc(n * n * sizeof(double));
        double* xp = calloc(n, sizeof(double));
        if (!expA || !Gi || !xp) { free(expA); free(Gi); free(xp); break; }

        ebc_matrix_exponential(A, n, mid, expA);
        ebc_matrix_exp_integral(A, n, mid, Gi);

        /* xp = expA * x_k */
        mat_vec(expA, (double*)x_k, xp, n);
        /* xp += Gi * BKx */
        if (BKx) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++)
                    xp[i] += Gi[i * n + j] * BKx[j];
            }
        }

        double err = 0.0;
        for (i = 0; i < n; i++) { double d = xp[i] - x_k[i]; err += d * d; }
        err = sqrt(err);

        free(expA); free(Gi); free(xp);

        if (err <= thresh) lo = mid;
        else               hi = mid;
    }

    free(BKx);
    return lo > 0.0 ? lo : tol;
}

/* ================================================================
 * Next event time for nonlinear STC via Lipschitz bound.
 *
 * Bound: |x(t) - x_k| <= (|f_k|/L) * (exp(L*tau) - 1)
 * Solving: tau = (1/L) * ln(1 + L*thresh/|f_k|)
 * ================================================================ */

double ebc_self_next_time_nonlinear(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx,
    const double* K, int m,
    const double* x_k, int n,
    double L, double sigma, double epsilon, double tau_max) {
    if (!f || !x_k || n < 1 || L < 1e-12 || tau_max <= 0.0) return 0.001;
    int i, j;

    double norm_xk = vec_norm(x_k, n);
    if (norm_xk < 1e-12) return tau_max;
    double thresh = sigma * norm_xk + epsilon;

    /* Compute f_k = f(0, x_k, K*x_k) */
    double* u = calloc(m, sizeof(double));
    double* fk = calloc(n, sizeof(double));
    if (u && fk && K) {
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                u[i] += K[i * n + j] * x_k[j];
        f(0.0, x_k, u, n, fk, ctx);
    }
    double norm_fk = (fk) ? vec_norm(fk, n) : L * norm_xk;
    free(u); free(fk);

    if (norm_fk < 1e-12) return tau_max;
    double tau = (1.0 / L) * log(1.0 + L * thresh / norm_fk);
    return (tau < tau_max) ? tau : tau_max;
}

/* ================================================================
 * Full self-triggered simulation
 * ================================================================ */

int ebc_self_simulate(EBC_System* sys, const EBC_Controller* ctrl,
                       double T, double dt,
                       const EBC_TriggerParams* tp,
                       double** traj, int* traj_len,
                       double** events, int* evt_len) {
    if (!sys || !ctrl || !tp || T <= 0 || dt <= 0) return -1;
    int n = sys->n, mcap = (int)(T / dt) + 10;

    *traj = calloc(mcap * n, sizeof(double));
    *events = calloc(mcap / 10 + 10, sizeof(double));
    if (!*traj || !*events) { free(*traj); free(*events); *traj = NULL; *events = NULL; return -1; }

    for (int i = 0; i < n; i++) (*traj)[i] = sys->x[i];
    int step = 1, ei = 0;
    (*events)[ei++] = 0.0;

    double t_next = dt;
    for (double t = dt; t <= T; t += dt) {
        if (t >= t_next - 1e-12 && ebc_check_event(sys, tp)) {
            ebc_mark_event(sys, t);
            if (ei < mcap / 10 + 10) (*events)[ei++] = t;
            t_next = t + tp->epsilon / 10.0;
            if (t_next > T) t_next = T + 1.0;
        }
        for (int i = 0; i < n; i++) { sys->x[i] += dt * sys->dx[i]; sys->e[i] = sys->x_last[i] - sys->x[i]; }
        sys->t = t;
        if (step < mcap) for (int i = 0; i < n; i++) (*traj)[step * n + i] = sys->x[i];
        step++;
    }
    *traj_len = step; *evt_len = ei;
    return 0;
}

/* ================================================================
 * Comparison of STC, ETC, and Periodic control
 * ================================================================ */

EBC_ComparisonResult ebc_compare_schemes(
    void (*f)(double, const double*, const double*, int, double*, void*),
    void* ctx, int n, int m,
    const double* A, const double* K,
    const double* x0, double T, double dt,
    double sigma, double epsilon, double period_h) {
    EBC_ComparisonResult res;
    memset(&res, 0, sizeof(res));

    EBC_System* s = ebc_system_create(n, m, EBC_CONTINUOUS_ETC);
    if (s) {
        ebc_system_set_state(s, x0);
        EBC_TriggerParams tp = ebc_trigger_mixed(sigma, epsilon);
        EBC_Controller c; memset(&c, 0, sizeof(c));
        c.K = (double*)K; c.n = n; c.m = m;
        double *tr, *ev; int tl, el;
        if (ebc_simulate(s, &c, T, dt, &tp, &tr, &tl, &ev, &el) == 0) {
            res.etc.total_events = s->event_count;
            res.etc.avg_iet = (el > 1) ? T / s->event_count : T;
            free(tr); free(ev);
        }
        ebc_system_free(s);
    }
    res.periodic.total_events = (period_h > 0) ? T / period_h : 0;
    res.stc.total_events = res.etc.total_events;
    if (res.periodic.total_events > 0) {
        res.etc.comm_reduction = 1.0 - res.etc.total_events / res.periodic.total_events;
    }
    return res;
}