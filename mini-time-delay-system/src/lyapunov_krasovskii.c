#include "lyapunov_krasovskii.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

#define LK_EPS 1e-12

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Symmetric matrix product: C = A * B * Aᵀ (n×n) */
__attribute__((unused))
static void sym_triple_product(const double* A, const double* B,
                                int n, double* C) {
    /* temp = B * Aᵀ */
    double* temp = (double*)calloc((size_t)(n * n), sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                temp[i * n + j] += B[i * n + k] * A[j * n + k];
    /* result = A * temp */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            C[i * n + j] = 0.0;
            for (int k = 0; k < n; k++)
                C[i * n + j] += A[i * n + k] * temp[k * n + j];
        }
    free(temp);
}

/* Matrix addition: C = A + B */
__attribute__((unused))
static void mat_add(const double* A, const double* B, int n, double* C) {
    int n2 = n * n;
    for (int i = 0; i < n2; i++) C[i] = A[i] + B[i];
}

/* Matrix subtract: C = A - B */
__attribute__((unused))
static void mat_sub(const double* A, const double* B, int n, double* C) {
    int n2 = n * n;
    for (int i = 0; i < n2; i++) C[i] = A[i] - B[i];
}

/* Eigenvalues of symmetric 2×2 matrix using closed form.
 * For n>2, use Jacobi iteration. */
static int sym_eigenvalues(const double* A, int n, double* eigvals) {
    if (n == 1) { eigvals[0] = A[0]; return 1; }
    if (n == 2) {
        double a11 = A[0], a12 = A[1], a22 = A[3];
        double trace = a11 + a22;
        double det = a11 * a22 - a12 * a12;
        double disc = trace * trace - 4.0 * det;
        if (disc < 0) disc = 0;
        double sqrt_disc = sqrt(disc);
        eigvals[0] = 0.5 * (trace + sqrt_disc);
        eigvals[1] = 0.5 * (trace - sqrt_disc);
        return 2;
    }

    /* Jacobi eigenvalue algorithm for n>2 */
    int max_iter = 100;
    double* V = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* work = (double*)malloc((size_t)(n * n) * sizeof(double));
    memcpy(work, A, (size_t)(n * n) * sizeof(double));

    /* Initialize V to identity */
    for (int i = 0; i < n * n; i++) V[i] = 0.0;
    for (int i = 0; i < n; i++) V[i * n + i] = 1.0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Find largest off-diagonal element */
        int p = 0, q = 1;
        double max_off = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                double v = fabs(work[i * n + j]);
                if (v > max_off) { max_off = v; p = i; q = j; }
            }
        if (max_off < LK_EPS) break;

        /* Compute rotation */
        double theta;
        double app = work[p * n + p];
        double aqq = work[q * n + q];
        double apq = work[p * n + q];
        if (fabs(app - aqq) < LK_EPS) {
            theta = M_PI / 4.0;
        } else {
            theta = 0.5 * atan2(2.0 * apq, app - aqq);
        }

        double c = cos(theta), s = sin(theta);

        /* Apply rotation to work matrix */
        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                double w_ip = work[i * n + p];
                double w_iq = work[i * n + q];
                work[i * n + p] = c * w_ip - s * w_iq;
                work[p * n + i] = work[i * n + p];
                work[i * n + q] = s * w_ip + c * w_iq;
                work[q * n + i] = work[i * n + q];
            }
        }
        work[p * n + p] = c * c * app + s * s * aqq - 2.0 * c * s * apq;
        work[q * n + q] = s * s * app + c * c * aqq + 2.0 * c * s * apq;
        work[p * n + q] = 0.0;
        work[q * n + p] = 0.0;
    }

    /* Extract eigenvalues */
    for (int i = 0; i < n; i++) eigvals[i] = work[i * n + i];

    /* Simple bubble sort (descending) */
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if (eigvals[i] < eigvals[j]) {
                double tmp = eigvals[i];
                eigvals[i] = eigvals[j];
                eigvals[j] = tmp;
            }

    free(V); free(work);
    return n;
}

/* ============================================================================
 * LK Functional — Construction
 * ============================================================================ */

LKFunctional* lkf_create(int n) {
    LKFunctional* lkf = (LKFunctional*)calloc(1, sizeof(LKFunctional));
    if (!lkf) return NULL;
    lkf->n = n;
    lkf->P = (double*)calloc((size_t)(n * n), sizeof(double));
    lkf->Q = (double*)calloc((size_t)(n * n), sizeof(double));
    lkf->R = (double*)calloc((size_t)(n * n), sizeof(double));
    lkf->S = (double*)calloc((size_t)(n * n), sizeof(double));
    return lkf;
}

void lkf_free(LKFunctional* lkf) {
    if (!lkf) return;
    free(lkf->P); free(lkf->Q); free(lkf->R); free(lkf->S);
    free(lkf);
}

void lkf_set_P(LKFunctional* lkf, const double* P_data) {
    if (!lkf || !P_data) return;
    memcpy(lkf->P, P_data, (size_t)(lkf->n * lkf->n) * sizeof(double));
}

void lkf_set_Q(LKFunctional* lkf, const double* Q_data) {
    if (!lkf || !Q_data) return;
    memcpy(lkf->Q, Q_data, (size_t)(lkf->n * lkf->n) * sizeof(double));
}

void lkf_set_R(LKFunctional* lkf, const double* R_data) {
    if (!lkf || !R_data) return;
    memcpy(lkf->R, R_data, (size_t)(lkf->n * lkf->n) * sizeof(double));
}

void lkf_set_identity(LKFunctional* lkf) {
    if (!lkf) return;
    for (int i = 0; i < lkf->n; i++) {
        lkf->P[i * lkf->n + i] = 1.0;
        lkf->Q[i * lkf->n + i] = 1.0;
        lkf->R[i * lkf->n + i] = 1.0;
    }
}

/* ============================================================================
 * LK Functional — Evaluation
 *
 * V(t, x_t) = x(t)ᵀ P x(t)
 *           + ∫_{t-τ}^{t} x(s)ᵀ Q x(s) ds
 *           + ∫_{-τ}^{0} ∫_{t+θ}^{t} ẋ(s)ᵀ R ẋ(s) ds dθ
 *
 * Discrete approximation:
 *   ∫_{t-τ}^{t} x(s)ᵀ Q x(s) ds ≈ Σᵢ x(t_i)ᵀ Q x(t_i) · Δt
 *   ∫∫ ẋ(s)ᵀ R ẋ(s) ds dθ ≈ Σᵢ Σⱼ ẋ(t_j)ᵀ R ẋ(t_j) · Δθ · Δt
 * ============================================================================ */

double lkf_evaluate(const LKFunctional* lkf,
                    const double* x_current,
                    const double* x_history,
                    const double* xdot_history,
                    int n_hist, double tau, double dt) {
    if (!lkf || !x_current) return 0.0;
    int n = lkf->n;

    /* Term 1: xᵀ P x */
    double V = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V += x_current[i] * lkf->P[i * n + j] * x_current[j];

    /* Term 2: ∫_{t-τ}^{t} xᵀ Q x ds
     * ≈ Σ_{k=0}^{n_hist-1} x_kᵀ Q x_k * dt */
    if (x_history && n_hist > 0) {
        for (int k = 0; k < n_hist; k++) {
            const double* xk = x_history + (size_t)k * n;
            double quad = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    quad += xk[i] * lkf->Q[i * n + j] * xk[j];
            V += quad * dt;
        }
    }

    /* Term 3: ∫_{-τ}^{0} ∫_{t+θ}^{t} ẋᵀ R ẋ ds dθ
     * ≈ Σ_{m=0}^{n_hist-1} Σ_{k=m}^{n_hist-1} ẋ_kᵀ R ẋ_k * dt * dt
     * ≈ (1/2) Σ_{k=0}^{n_hist-1} (k+1)(n_hist-k) ẋ_kᵀ R ẋ_k * dt²
     * Simplified: τ * Σ ẋ_kᵀ R ẋ_k * dt */
    if (xdot_history && n_hist > 0) {
        double sum_quad = 0.0;
        for (int k = 0; k < n_hist; k++) {
            const double* xdk = xdot_history + (size_t)k * n;
            double quad = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    quad += xdk[i] * lkf->R[i * n + j] * xdk[j];
            /* Trapezoidal weight approximation */
            double weight = (k == 0 || k == n_hist - 1) ? 0.5 * dt * dt * tau
                                                        : dt * dt * tau;
            sum_quad += quad * weight;
        }
        V += sum_quad;
    }

    return V;
}

/* ============================================================================
 * LK Functional — Derivative
 *
 * For ẋ = A x + A_d x(t-τ):
 * dV/dt = 2xᵀP(Ax+A_d x_d)
 *       + xᵀQ x - x_dᵀQ x_d
 *       + τ ẋᵀR ẋ - ∫_{t-τ}^{t} ẋ(s)ᵀR ẋ(s) ds
 *
 * Using Jensen's inequality: -∫ ẋᵀR ẋ ≤ -(1/τ)(x-x_d)ᵀR(x-x_d)
 * ============================================================================ */

double lkf_derivative(const LKFunctional* lkf,
                      const TimeDelaySystem* sys,
                      const double* x_current,
                      const double* x_delayed,
                      const double* xdot_current,
                      const double* xdot_history,
                      int n_hist, double tau, double dt) {
    if (!lkf || !sys || !x_current || !x_delayed) return 0.0;
    int n = lkf->n;

    /* Compute xdot = A x + A_d x_d (if not provided) */
    double* xdot = (double*)malloc((size_t)n * sizeof(double));
    if (xdot_current) {
        memcpy(xdot, xdot_current, (size_t)n * sizeof(double));
    } else {
        for (int i = 0; i < n; i++) {
            xdot[i] = 0.0;
            for (int j = 0; j < n; j++) {
                xdot[i] += sys->A[i * n + j] * x_current[j]
                         + sys->A_delayed[i * n + j] * x_delayed[j];
            }
        }
    }

    /* Term 1: 2 xᵀ P (A x + A_d x_d) */
    double dV = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            dV += 2.0 * x_current[i] * lkf->P[i * n + j] * xdot[j];

    /* Term 2: xᵀ Q x - x_dᵀ Q x_d */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            dV += x_current[i] * lkf->Q[i * n + j] * x_current[j];
            dV -= x_delayed[i] * lkf->Q[i * n + j] * x_delayed[j];
        }

    /* Term 3: τ ẋᵀ R ẋ - ∫ ẋᵀ R ẋ ds
     * Using Jensen: -∫ ẋᵀRẋ ≤ -(1/τ)(x-x_d)ᵀR(x-x_d) */
    double xdot_R_xdot = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            xdot_R_xdot += xdot[i] * lkf->R[i * n + j] * xdot[j];
    dV += tau * xdot_R_xdot;

    /* Jensen bound for the integral term */
    double* diff = (double*)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) diff[i] = x_current[i] - x_delayed[i];
    double diff_R_diff = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            diff_R_diff += diff[i] * lkf->R[i * n + j] * diff[j];
    if (tau > LK_EPS) dV -= diff_R_diff / tau;

    /* Optional: use history to refine the integral term */
    if (xdot_history && n_hist > 1 && dt > 0) {
        double integral_term = 0.0;
        for (int k = 0; k < n_hist; k++) {
            const double* xdk = xdot_history + (size_t)k * n;
            double quad = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    quad += xdk[i] * lkf->R[i * n + j] * xdk[j];
            integral_term += quad * dt;
        }
        /* Blend Jensen bound with discretized integral */
        dV -= 0.5 * integral_term + 0.5 * diff_R_diff / tau;
    }

    free(xdot); free(diff);
    return dV;
}

/* ============================================================================
 * L4 — Positive Definiteness Check
 * ============================================================================ */

bool lkf_is_positive_definite(const LKFunctional* lkf) {
    if (!lkf) return false;
    int n = lkf->n;
    double* eig = (double*)malloc((size_t)n * sizeof(double));

    /* Check P > 0 */
    sym_eigenvalues(lkf->P, n, eig);
    for (int i = 0; i < n; i++)
        if (eig[i] <= LK_EPS) { free(eig); return false; }

    /* Check Q > 0 */
    sym_eigenvalues(lkf->Q, n, eig);
    for (int i = 0; i < n; i++)
        if (eig[i] <= LK_EPS) { free(eig); return false; }

    /* Check R > 0 */
    sym_eigenvalues(lkf->R, n, eig);
    for (int i = 0; i < n; i++)
        if (eig[i] <= LK_EPS) { free(eig); return false; }

    free(eig);
    return true;
}

/* ============================================================================
 * L4 — Decay Rate (negative-definite rate)
 * ============================================================================ */

double lkf_decay_rate(const LKFunctional* lkf,
                      const TimeDelaySystem* sys,
                      const double* x_current,
                      const double* x_delayed,
                      double tau) {
    if (!lkf || !sys || !x_current || !x_delayed) return 0.0;

    (void)lkf_evaluate(lkf, x_current, NULL, NULL, 0, tau, 0.0);
    double dV = lkf_derivative(lkf, sys, x_current, x_delayed,
                                NULL, NULL, 0, tau, 0.0);

    /* dV ≤ -ε ||x||² → ε = -dV / (xᵀx + ε0) */
    double x2 = 0.0;
    for (int i = 0; i < lkf->n; i++) x2 += x_current[i] * x_current[i];
    if (x2 < LK_EPS) return 0.0;
    return -dV / x2;
}

/* ============================================================================
 * Razumikhin Condition Check
 *
 * Theorem (Razumikhin): If there exists a Lyapunov FUNCTION V(x) and
 * functions α, β, γ ∈ 𝒦 such that:
 *   α(||x||) ≤ V(x) ≤ β(||x||)
 *   V̇(x(t)) ≤ -γ(||x(t)||)  whenever  V(x(t+θ)) ≤ p V(x(t)) ∀θ ∈ [-τ,0],
 *   for some p > 1.
 * Then the trivial solution is stable.
 *
 * Check: V(x(t+θ)) ≤ p V(x(t)) for all θ in the history buffer.
 * ============================================================================ */

bool razumikhin_condition_check(const double* x_history, int n_hist,
                                const double* x_current, int n,
                                double p_factor,
                                double (*V_func)(const double*, int)) {
    if (!x_history || !x_current || !V_func) return false;
    if (p_factor <= 1.0) return false;  /* p must be > 1 */

    double V_current = V_func(x_current, n);
    if (V_current < LK_EPS) return true;  /* Trivially holds at equilibrium */

    double threshold = p_factor * V_current;

    for (int k = 0; k < n_hist; k++) {
        double V_k = V_func(x_history + (size_t)k * n, n);
        if (V_k > threshold + LK_EPS) return false;
    }
    return true;
}

/* ============================================================================
 * LMI-Based Stability Checks
 *
 * Reference: Gu, Kharitonov, Chen (2003), Ch. 3
 *
 * Delay-Independent LMI:
 *   [ PA + AᵀP + Q    PA_d   ] < 0
 *   [   A_dᵀP         -Q    ]
 *
 * Delay-Dependent LMI:
 *   [ PA+AᵀP+Q  PA_d  τAᵀR ]   [ R  0   0 ]  [-I]           [-I  ]
 *   [  A_dᵀP    -Q   τA_dᵀR ] + [ 0  -R  0 ]  [ 0 ] R^(-1) [ 0  ]ᵀ < 0
 *   [  τRA     τRA_d  -τR  ]   [ 0   0   0 ]  [ I ]           [ I  ]
 *
 * Schur complement form of delay-dependent LMI:
 *   [ PA+AᵀP+Q+τAᵀRA  PA_d+τAᵀRA_d ]
 *   [ A_dᵀP+τA_dᵀRA  -Q+τA_dᵀRA_d ] < 0
 * ============================================================================ */

bool lmi_delay_independent_check(const TimeDelaySystem* sys,
                                 double* out_P, double* out_Q) {
    if (!sys || sys->n_states < 1) return false;
    /* Simplified: for scalar system, the condition reduces to checking
     * if P*a + a*P + Q < 0 and |a_d*P| < Q.
     * For 1D: 2pa + q < 0 and (p*a_d)^2 < q^2 → |a_d| < -a
     * This is equivalent to: a + |a_d| < 0.
     * So system is delay-independently stable iff a + |a_d| < 0. */

    int n = sys->n_states;
    if (n == 1) {
        double a = sys->A[0];
        double ad = sys->A_delayed[0];
        /* Choose P = 1, Q = |a_d| + ε */
        if (out_P) out_P[0] = 1.0;
        if (out_Q) out_Q[0] = fabs(ad) + 0.1;
        return (a + fabs(ad) < -LK_EPS);
    }

    /* For n>1: We need to verify existence.
     * Use a simple necessary condition: the delay-free system A+A_d
     * must be Hurwitz for any chance at delay-independent stability. */
    double* A_sum = (double*)malloc((size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n * n; i++)
        A_sum[i] = sys->A[i] + sys->A_delayed[i];

    /* Check if all eigenvalues of A+A_d have negative real parts.
     * This is a necessary (not sufficient!) condition. */
    double* eig = (double*)calloc((size_t)n, sizeof(double));
    sym_eigenvalues(A_sum, n, eig);

    bool all_neg = true;
    for (int i = 0; i < n; i++)
        if (eig[i] >= -LK_EPS) { all_neg = false; break; }

    /* If all_neg, set P = I, Q = some positive definite matrix */
    if (all_neg && out_P && out_Q) {
        for (int i = 0; i < n; i++) {
            out_P[i * n + i] = 1.0;
            out_Q[i * n + i] = 1.0;
        }
    }

    free(A_sum); free(eig);
    return all_neg;
}

bool lmi_delay_dependent_check(const TimeDelaySystem* sys,
                                double tau_max, double* out_P,
                                double* out_Q, double* out_R) {
    if (!sys || sys->n_states < 1 || tau_max <= 0) return false;
    int n = sys->n_states;

    if (n == 1) {
        double a = sys->A[0];
        double ad = sys->A_delayed[0];

        /* For scalar system, the characteristic equation is:
         * s + a + a_d e^{-τs} = 0
         *
         * Stability condition for τ = τ_max:
         * Find ω > 0 such that ω² = a_d² - a².
         * If a_d² ≤ a² → delay-independent stable.
         * If a_d² > a² → τ_max = (1/ω) [π - atan2(ω, -a)].
         *
         * Check if stable at given tau_max. */
        if (fabs(ad) <= fabs(a) && a < 0) {
            /* Delay-independently stable */
            if (out_P) out_P[0] = 1.0;
            if (out_Q) out_Q[0] = 1.0;
            if (out_R) out_R[0] = 1.0;
            return true;
        }

        /* Find crossing frequency */
        double omega2 = ad * ad - a * a;
        if (omega2 <= 0) return true;  /* No crossing */
        double omega = sqrt(omega2);

        /* Critical delay */
        double tau_crit = (M_PI - atan2(omega, -a)) / omega;
        if (tau_crit < 0) tau_crit = (2.0 * M_PI - atan2(omega, -a)) / omega;

        if (tau_max < tau_crit) {
            if (out_P) out_P[0] = 1.0;
            if (out_Q) out_Q[0] = tau_crit / (tau_max + LK_EPS);
            if (out_R) out_R[0] = 1.0 / tau_max;
            return true;
        }
        return false;
    }

    /* For n>1:
     * Use the delay-dependent LMI via Schur complement.
     * This is a simplified check — a full LMI solver would be needed
     * for guaranteed verification.
     *
     * Here we compute a sufficient condition by checking
     * ||(sI - A)^{-1} A_d|| < 1 for all s on the imaginary axis
     * (small-gain theorem approach). */

    /* Check at s = jω for a range of frequencies */
    bool feasible = true;
    for (int k = 0; k < 100 && feasible; k++) {
        double omega = 0.01 * (double)(k + 1);

        /* Compute (jωI - A)^{-1} */
        /* For simplicity, check if A is Hurwitz and A_d is "small enough" */

        /* Estimate via norm condition:
         * ||(jωI - A)^{-1} A_d|| ≤ ||(jωI - A)^{-1}|| · ||A_d|| */
        double Anorm = 0.0, Adnorm = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                Anorm += sys->A[i * n + j] * sys->A[i * n + j];
                Adnorm += sys->A_delayed[i * n + j] * sys->A_delayed[i * n + j];
            }
        Anorm = sqrt(Anorm); Adnorm = sqrt(Adnorm);

        /* Simple sufficient condition: if A is Hurwitz, check
         * Adnorm < sigma_min(jωI - A) for all ω. Approximate:
         * For ω small: sigma_min ≈ |Re(λ_min(A))|
         * For ω large: sigma_min ≈ ω */
        double A_sigma_min = 0.0;

        /* Compute |Re(λ_min(A))| */
        double* A_eig = (double*)calloc((size_t)n, sizeof(double));
        sym_eigenvalues(sys->A, n, A_eig);
        double max_re_eig = -INFINITY;
        for (int i = 0; i < n; i++)
            if (A_eig[i] > max_re_eig) max_re_eig = A_eig[i];
        if (max_re_eig >= 0) { free(A_eig); feasible = false; break; }

        A_sigma_min = -max_re_eig;
        if (omega > A_sigma_min) A_sigma_min = omega;
        free(A_eig);

        if (Adnorm >= A_sigma_min) { feasible = false; break; }
    }

    if (feasible && out_P && out_Q && out_R) {
        for (int i = 0; i < n; i++) {
            out_P[i * n + i] = 1.0;
            out_Q[i * n + i] = 1.0;
            out_R[i * n + i] = 1.0;
        }
    }

    return feasible;
}

/* ============================================================================
 * Augmented LK Functional
 * ============================================================================ */

AugmentedLKFunctional* alkf_create(int n, double tau_min,
                                    double tau_max, double mu) {
    AugmentedLKFunctional* alkf = (AugmentedLKFunctional*)
        calloc(1, sizeof(AugmentedLKFunctional));
    if (!alkf) return NULL;
    alkf->base.n = n;
    alkf->base.P = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->base.Q = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->base.R = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->base.S = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->Z1 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->Z2 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->N1 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->N2 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->M1 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->M2 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->M3 = (double*)calloc((size_t)(n * n), sizeof(double));
    alkf->tau_min = tau_min;
    alkf->tau_max = tau_max;
    alkf->mu = mu;
    return alkf;
}

void alkf_free(AugmentedLKFunctional* alkf) {
    if (!alkf) return;
    free(alkf->base.P); free(alkf->base.Q); free(alkf->base.R);
    free(alkf->base.S);
    free(alkf->Z1); free(alkf->Z2);
    free(alkf->N1); free(alkf->N2);
    free(alkf->M1); free(alkf->M2); free(alkf->M3);
    free(alkf);
}

double alkf_evaluate(const AugmentedLKFunctional* alkf,
                     const double* x_current,
                     const double* x_history,
                     const double* xdot_history,
                     int n_hist, double tau_t, double dt) {
    /* Base LK functional evaluation plus augmented terms */
    double V_base = lkf_evaluate(&alkf->base, x_current,
                                  x_history, xdot_history,
                                  n_hist, tau_t, dt);
    /* Additional interval-dependent terms for time-varying delay range.
     * These augment the base functional with extra slack terms. */
    return V_base;  /* Augmented terms would go here for full LMI-based V */
}

/* ============================================================================
 * Discretized LK Functional (Gu-Kharitonov-Chen)
 * ============================================================================ */

DiscretizedLKFunctional* dlkf_create(int n, int N, double tau) {
    DiscretizedLKFunctional* dlkf = (DiscretizedLKFunctional*)
        calloc(1, sizeof(DiscretizedLKFunctional));
    if (!dlkf) return NULL;
    dlkf->n = n;
    dlkf->N = N;
    dlkf->tau = tau;
    dlkf->P = (double*)calloc((size_t)(n * n), sizeof(double));
    dlkf->Q_mesh = (double*)calloc((size_t)((N + 1) * n * n), sizeof(double));
    dlkf->R_mesh = (double*)calloc((size_t)((N + 1) * (N + 1) * n * n),
                                    sizeof(double));
    dlkf->S_mesh = (double*)calloc((size_t)((N + 1) * n * n), sizeof(double));
    return dlkf;
}

void dlkf_free(DiscretizedLKFunctional* dlkf) {
    if (!dlkf) return;
    free(dlkf->P); free(dlkf->Q_mesh); free(dlkf->R_mesh);
    free(dlkf->S_mesh);
    free(dlkf);
}

double dlkf_evaluate(const DiscretizedLKFunctional* dlkf,
                     const double* x_current,
                     const double* x_history,
                     int n_hist, double dt) {
    if (!dlkf || !x_current) return 0.0;
    int n = dlkf->n;

    /* V = xᵀ P x + 2 xᵀ Σᵢ Q_i ∫ x dξ + Σᵢⱼ ∫∫ xᵀ R_{ij} x + Σᵢ ∫ xᵀ S_i x */
    double V = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            V += x_current[i] * dlkf->P[i * n + j] * x_current[j];

    /* Discretized integral terms using Q_mesh and S_mesh */
    if (x_history && n_hist > 0 && dt > 0.0) {
        for (int m = 0; m < n_hist; m++) {
            const double* xk = x_history + (size_t)m * n;
            /* S_mesh interpolation (nearest mesh point) */
            int mesh_idx = (int)((double)m / (double)n_hist * dlkf->N);
            if (mesh_idx > dlkf->N) mesh_idx = dlkf->N;

            double* Si = dlkf->S_mesh + (size_t)mesh_idx * n * n;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    V += xk[i] * Si[i * n + j] * xk[j] * dt;
        }
    }
    return V;
}
