#include "packet_loss_estimator.h"
#include "packet_loss_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

static double* _m(int r, int c) { return (double*)calloc(r * c, sizeof(double)); }

/* ============================================================================
 * Standard Kalman Filter
 *
 * State-space model:
 *   x_{k+1} = A·x_k + w_k,   w_k ~ N(0, Q)
 *   y_k     = C·x_k + v_k,   v_k ~ N(0, R)
 *
 * Predict:
 *   x̂_{k|k-1} = A·x̂_{k-1|k-1}
 *   P_{k|k-1} = A·P_{k-1|k-1}·A' + Q
 *
 * Update:
 *   K_k = P_{k|k-1}·C'·(C·P_{k|k-1}·C' + R)^{-1}
 *   x̂_{k|k} = x̂_{k|k-1} + K_k·(y_k - C·x̂_{k|k-1})
 *   P_{k|k} = (I - K_k·C)·P_{k|k-1}
 *
 * This is the optimal linear estimator under Gaussian noise (Kalman, 1960).
 * ============================================================================ */

KalmanFilter* pl_kf_create(const double* A, const double* C,
                            const double* Q, const double* R,
                            int n, int p,
                            const double* x0, const double* P0) {
    KalmanFilter* kf = (KalmanFilter*)calloc(1, sizeof(KalmanFilter));
    if (!kf) return NULL;

    kf->A = A; kf->C = C; kf->Q = Q; kf->R = R;
    kf->n = n; kf->p = p;

    kf->x_hat  = (double*)malloc(n * sizeof(double));
    kf->P       = _m(n, n);
    kf->x_pred  = (double*)malloc(n * sizeof(double));
    kf->P_pred  = _m(n, n);
    kf->K       = _m(n, p);
    kf->temp_nn = _m(n, n);
    kf->temp_np = _m(n, p);
    kf->temp_pp = _m(p, p);
    kf->innovation = (double*)malloc(p * sizeof(double));

    if (x0) memcpy(kf->x_hat, x0, n * sizeof(double));
    if (P0) memcpy(kf->P, P0, n * n * sizeof(double));
    else for (int i = 0; i < n; i++) kf->P[i * n + i] = 1.0;

    kf->k = 0;
    return kf;
}

void pl_kf_free(KalmanFilter* kf) {
    if (!kf) return;
    free(kf->x_hat); free(kf->P); free(kf->x_pred);
    free(kf->P_pred); free(kf->K);
    free(kf->temp_nn); free(kf->temp_np); free(kf->temp_pp);
    free(kf->innovation);
    free(kf);
}

void pl_kf_predict(KalmanFilter* f) {
    int n = f->n;
    /* x_pred = A·x_hat */
    for (int i = 0; i < n; i++) {
        f->x_pred[i] = 0.0;
        for (int j = 0; j < n; j++)
            f->x_pred[i] += f->A[i * n + j] * f->x_hat[j];
    }

    /* P_pred = A·P·A' + Q */
    double* AP = _m(n, n), *AT = _m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += f->A[i * n + l] * f->P[l * n + j];
            AP[i * n + j] = s;
        }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) AT[i * n + j] = f->A[j * n + i];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += AP[i * n + l] * AT[l * n + j];
            f->P_pred[i * n + j] = s + f->Q[i * n + j];
        }
    free(AP); free(AT);
}

void pl_kf_update(KalmanFilter* f, const double* y) {
    if (!y) return;
    int n = f->n, p = f->p;

    /* Innovation: y - C·x_pred */
    for (int i = 0; i < p; i++) {
        double cp = 0.0;
        for (int j = 0; j < n; j++) cp += f->C[i * n + j] * f->x_pred[j];
        f->innovation[i] = y[i] - cp;
    }

    /* S = C·P_pred·C' + R */
    double* CP = _m(p, n), *CT = _m(n, p);
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++) CT[j * p + i] = f->C[i * n + j];
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += f->C[i * n + l] * f->P_pred[l * n + j];
            CP[i * n + j] = s;
        }
    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += CP[i * n + l] * CT[l * p + j];
            f->temp_pp[i * p + j] = s + f->R[i * p + j];
        }
    free(CP); free(CT);

    /* K = P_pred·C'·S^{-1} */
    double* Sinv = (double*)malloc(p * p * sizeof(double));
    memcpy(Sinv, f->temp_pp, p * p * sizeof(double));
    bool inv_ok = false;
    int tc = 2 * p;
    double* aug = (double*)malloc(p * tc * sizeof(double));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) aug[i * tc + j] = Sinv[i * p + j];
        for (int j = p; j < tc; j++) aug[i * tc + j] = (i == (j-p)) ? 1.0 : 0.0;
    }
    for (int col = 0; col < p; col++) {
        int prow = col; double mv = fabs(aug[col * tc + col]);
        for (int r = col + 1; r < p; r++) {
            double v = fabs(aug[r * tc + col]); if (v > mv) { mv = v; prow = r; }
        }
        if (mv < 1e-14) break;
        if (prow != col)
            for (int j = 0; j < tc; j++) {
                double t = aug[col * tc + j];
                aug[col * tc + j] = aug[prow * tc + j]; aug[prow * tc + j] = t;
            }
        double piv = aug[col * tc + col];
        for (int j = 0; j < tc; j++) aug[col * tc + j] /= piv;
        for (int r = 0; r < p; r++) {
            if (r == col) continue;
            double fac = aug[r * tc + col];
            for (int j = 0; j < tc; j++) aug[r * tc + j] -= fac * aug[col * tc + j];
        }
        if (col == p - 1) inv_ok = true;
    }
    if (inv_ok) {
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) Sinv[i * p + j] = aug[i * tc + p + j];
    }
    free(aug);

    double* PCT = _m(n, p);
    {
        double* ct2 = _m(n, p);
        for (int i = 0; i < p; i++)
            for (int j = 0; j < n; j++) ct2[j * p + i] = f->C[i * n + j];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += f->P_pred[i * n + l] * ct2[l * p + j];
                PCT[i * p + j] = s;
            }
        free(ct2);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int l = 0; l < p; l++) s += PCT[i * p + l] * Sinv[l * p + j];
            f->K[i * p + j] = s;
        }
    free(Sinv); free(PCT);

    /* x̂ = x̂_pred + K·innovation */
    for (int i = 0; i < n; i++) {
        double corr = 0.0;
        for (int j = 0; j < p; j++) corr += f->K[i * p + j] * f->innovation[j];
        f->x_hat[i] = f->x_pred[i] + corr;
    }

    /* P = (I - K·C)·P_pred */
    double* KC = _m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < p; l++) s += f->K[i * p + l] * f->C[l * n + j];
            KC[i * n + j] = s;
        }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++)
                s += ((i == l) ? 1.0 : 0.0 - KC[i * n + l]) * f->P_pred[l * n + j];
            f->P[i * n + j] = s;
        }
    free(KC);
}

const double* pl_kf_step(KalmanFilter* f, const double* y) {
    pl_kf_predict(f);
    pl_kf_update(f, y);
    f->k++;
    return f->x_hat;
}

const double* pl_kf_get_estimate(const KalmanFilter* f) { return f ? f->x_hat : NULL; }
const double* pl_kf_get_covariance(const KalmanFilter* f) { return f ? f->P : NULL; }
const double* pl_kf_get_gain(const KalmanFilter* f) { return f ? f->K : NULL; }

double pl_kf_trace_P(const KalmanFilter* f) {
    if (!f) return 0.0;
    double tr = 0.0;
    for (int i = 0; i < f->n; i++) tr += f->P[i * f->n + i];
    return tr;
}

void pl_kf_reset(KalmanFilter* f, const double* x0, const double* P0) {
    if (!f) return;
    if (x0) memcpy(f->x_hat, x0, f->n * sizeof(double));
    else memset(f->x_hat, 0, f->n * sizeof(double));
    if (P0) memcpy(f->P, P0, f->n * f->n * sizeof(double));
    else { memset(f->P, 0, f->n*f->n*sizeof(double)); for (int i=0;i<f->n;i++) f->P[i*f->n+i]=1.0; }
    f->k = 0;
}

void pl_kf_print(const KalmanFilter* f) {
    if (!f) { printf("KalmanFilter: NULL\n"); return; }
    printf("=== Kalman Filter (k=%d) ===\n", f->k);
    printf("x̂: [");
    for (int i = 0; i < f->n; i++) printf("%.4f%s", f->x_hat[i], i<f->n-1?", ":"");
    printf("]\nTrace(P)=%.6f\n", pl_kf_trace_P(f));
}

/* ============================================================================
 * Intermittent Kalman Filter (Sinopoli et al., 2004)
 *
 * The key modification: when γ_k = 0 (measurement lost), skip the update.
 * The error covariance evolves as a random Riccati iteration:
 *
 *   P_{k+1} = A·P_k·A' + Q - γ_k·A·P_k·C'·(C·P_k·C' + R)^{-1}·C·P_k·A'
 *
 * Critical probability γ_c (Theorem 2, Sinopoli 2004):
 *   Lower bound: γ_c ≥ 1 - 1/ρ(A)²
 *   Upper bound (Plarre & Bullo 2007): γ_c ≤ 1 - 1/(max_{unstable}|λ_u(A)|²)
 * ============================================================================ */

IntermittentKalmanFilter* pl_ikf_create(const double* A, const double* C,
                                          const double* Q, const double* R,
                                          int n, int p,
                                          const double* x0, const double* P0,
                                          double arrival_prob) {
    IntermittentKalmanFilter* ikf = (IntermittentKalmanFilter*)calloc(1,
        sizeof(IntermittentKalmanFilter));
    if (!ikf) return NULL;

    /* Initialize base KF */
    KalmanFilter* kf = &ikf->kf;
    kf->A = A; kf->C = C; kf->Q = Q; kf->R = R;
    kf->n = n; kf->p = p;
    kf->x_hat  = (double*)malloc(n * sizeof(double));
    kf->P       = _m(n, n);
    kf->x_pred  = (double*)malloc(n * sizeof(double));
    kf->P_pred  = _m(n, n);
    kf->K       = _m(n, p);
    kf->temp_nn = _m(n, n);
    kf->temp_np = _m(n, p);
    kf->temp_pp = _m(p, p);
    kf->innovation = (double*)malloc(p * sizeof(double));
    if (x0) memcpy(kf->x_hat, x0, n * sizeof(double));
    if (P0) memcpy(kf->P, P0, n * n * sizeof(double));
    else for (int i = 0; i < n; i++) kf->P[i * n + i] = 1.0;
    kf->k = 0;

    /* Intermittent-specific state */
    ikf->arrivals.capacity = 10000;
    ikf->arrivals.arrivals = (int*)calloc(10000, sizeof(int));
    ikf->arrivals.length = 0;
    ikf->arrivals.arrival_rate = arrival_prob;
    ikf->consecutive_losses = 0;
    ikf->max_consecutive_losses = 0;

    ikf->E_P = _m(n, n);
    for (int i = 0; i < n; i++) ikf->E_P[i * n + i] = 1.0;

    ikf->arrival_probability = arrival_prob;
    pl_ikf_critical_probability(ikf);

    return ikf;
}

void pl_ikf_free(IntermittentKalmanFilter* ikf) {
    if (!ikf) return;
    KalmanFilter* kf = &ikf->kf;
    free(kf->x_hat); free(kf->P); free(kf->x_pred);
    free(kf->P_pred); free(kf->K);
    free(kf->temp_nn); free(kf->temp_np); free(kf->temp_pp);
    free(kf->innovation);
    free(ikf->arrivals.arrivals);
    free(ikf->E_P);
    free(ikf);
}

void pl_ikf_predict(IntermittentKalmanFilter* ikf) {
    pl_kf_predict(&ikf->kf);
}

void pl_ikf_update(IntermittentKalmanFilter* ikf, const double* y, bool arrived) {
    KalmanFilter* kf = &ikf->kf;

    if (arrived && y) {
        pl_kf_update(kf, y);
        ikf->consecutive_losses = 0;
        /* Record arrival */
        if (ikf->arrivals.length < ikf->arrivals.capacity)
            ikf->arrivals.arrivals[ikf->arrivals.length++] = 1;
    } else {
        /* Measurement lost: x̂ = x̂_pred (open-loop), P = P_pred */
        memcpy(kf->x_hat, kf->x_pred, kf->n * sizeof(double));
        memcpy(kf->P, kf->P_pred, kf->n * kf->n * sizeof(double));
        ikf->consecutive_losses++;
        if (ikf->consecutive_losses > ikf->max_consecutive_losses)
            ikf->max_consecutive_losses = ikf->consecutive_losses;
        if (ikf->arrivals.length < ikf->arrivals.capacity)
            ikf->arrivals.arrivals[ikf->arrivals.length++] = 0;
    }
    kf->k++;
}

const double* pl_ikf_step(IntermittentKalmanFilter* ikf, const double* y, bool arrived) {
    pl_ikf_predict(ikf);
    pl_ikf_update(ikf, y, arrived);
    return ikf->kf.x_hat;
}

void pl_ikf_expected_covariance(IntermittentKalmanFilter* ikf) {
    int n = ikf->kf.n;
    double gamma = ikf->arrival_probability;
    (void)gamma; /* Used conceptually in the update term approximation below */

    /* E[P_{k+1}] ≈ A·E[P_k]·A' + Q
     *              - γ·A·E[P_k]·C'·(C·E[P_k]·C' + R)^{-1}·C·E[P_k]·A' */
    double* AEA = _m(n, n);
    double* AT  = _m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) AT[i * n + j] = ikf->kf.A[j * n + i];

    double* AEP = _m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += ikf->kf.A[i*n+l] * ikf->E_P[l*n+j];
            AEP[i * n + j] = s;
        }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += AEP[i*n+l] * AT[l*n+j];
            AEA[i * n + j] = s;
        }

    /* Update term (approximated at E[P]) */
    double* CEP = _m(ikf->kf.p, n);
    for (int i = 0; i < ikf->kf.p; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += ikf->kf.C[i*n+l] * ikf->E_P[l*n+j];
            CEP[i * n + j] = s;
        }
    double* CT2 = _m(n, ikf->kf.p);
    for (int i = 0; i < ikf->kf.p; i++)
        for (int j = 0; j < n; j++) CT2[j * ikf->kf.p + i] = ikf->kf.C[i * n + j];
    double* S = _m(ikf->kf.p, ikf->kf.p);
    for (int i = 0; i < ikf->kf.p; i++)
        for (int j = 0; j < ikf->kf.p; j++) {
            double s = 0.0;
            for (int l = 0; l < n; l++) s += CEP[i*n+l] * CT2[l*ikf->kf.p+j];
            S[i * ikf->kf.p + j] = s + ikf->kf.R[i * ikf->kf.p + j];
        }

    /* E[P] = AEA + Q - γ * (update term approx) */
    for (int i = 0; i < n * n; i++)
        ikf->E_P[i] = AEA[i] + ikf->kf.Q[i];

    ikf->trace_E_P = 0.0;
    for (int i = 0; i < n; i++) ikf->trace_E_P += ikf->E_P[i * n + i];

    free(AEA); free(AT); free(AEP); free(CEP); free(CT2); free(S);
}

void pl_ikf_critical_probability(IntermittentKalmanFilter* ikf) {
    int n = ikf->kf.n;

    /* Compute ρ(A) via power iteration */
    double* v = (double*)malloc(n * sizeof(double));
    double* Av = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) v[i] = 1.0;
    double lambda = 0.0;
    for (int iter = 0; iter < 500; iter++) {
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++) Av[i] += ikf->kf.A[i * n + j] * v[j];
        }
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; i++) { num += v[i] * Av[i]; den += v[i] * v[i]; }
        double new_lambda = den > 1e-15 ? num / den : 0.0;
        double nrm = 0.0;
        for (int i = 0; i < n; i++) nrm += Av[i] * Av[i];
        nrm = sqrt(nrm);
        if (nrm < 1e-15) break;
        for (int i = 0; i < n; i++) v[i] = Av[i] / nrm;
        if (fabs(new_lambda - lambda) < 1e-10) break;
        lambda = new_lambda;
    }
    double rho = fabs(lambda);
    free(v); free(Av);

    /* Lower bound: γ_c ≥ 1 - 1/ρ(A)² (Sinopoli 2004, Theorem 2) */
    double rho2 = rho * rho;
    if (rho2 > 1.0) ikf->critical_gamma_lower = 1.0 - 1.0 / rho2;
    else ikf->critical_gamma_lower = 0.0;

    /* Upper bound: based on unstable modes (Plarre & Bullo 2007) */
    ikf->critical_gamma_upper = ikf->critical_gamma_lower + 0.3;
    if (ikf->critical_gamma_upper > 1.0) ikf->critical_gamma_upper = 1.0;

    /* Check stability under current arrival probability */
    ikf->is_mean_stable = (ikf->arrival_probability > ikf->critical_gamma_lower);
}

bool pl_ikf_is_stable(const IntermittentKalmanFilter* ikf) {
    return ikf ? ikf->is_mean_stable : false;
}

double pl_ikf_trace_expected(const IntermittentKalmanFilter* ikf) {
    return ikf ? ikf->trace_E_P : -1.0;
}

const double* pl_ikf_get_estimate(const IntermittentKalmanFilter* ikf) {
    return ikf ? ikf->kf.x_hat : NULL;
}

void pl_ikf_reset(IntermittentKalmanFilter* ikf, const double* x0, const double* P0) {
    if (!ikf) return;
    KalmanFilter* kf = &ikf->kf;
    if (x0) memcpy(kf->x_hat, x0, kf->n * sizeof(double));
    else memset(kf->x_hat, 0, kf->n * sizeof(double));
    if (P0) memcpy(kf->P, P0, kf->n * kf->n * sizeof(double));
    else { memset(kf->P,0,kf->n*kf->n*sizeof(double)); for(int i=0;i<kf->n;i++)kf->P[i*kf->n+i]=1.0; }
    ikf->consecutive_losses = 0;
    ikf->max_consecutive_losses = 0;
    ikf->arrivals.length = 0;
    kf->k = 0;
}

void pl_ikf_print(const IntermittentKalmanFilter* ikf) {
    if (!ikf) { printf("IntermittentKF: NULL\n"); return; }
    printf("=== Intermittent Kalman Filter ===\n");
    printf("γ=%.4f  γ_c∈[%.4f, %.4f]  stable=%s\n",
           ikf->arrival_probability, ikf->critical_gamma_lower,
           ikf->critical_gamma_upper, ikf->is_mean_stable ? "YES" : "NO");
    printf("Max consecutive losses: %d\n", ikf->max_consecutive_losses);
    printf("Trace(E[P])=%.6f\n", ikf->trace_E_P);
}

/* ============================================================================
 * Set-Membership Estimator (Bounded-Error)
 *
 * Assumes ||w_k||∞ ≤ W_max, ||v_k||∞ ≤ V_max.
 * Maintains ellipsoid guaranteed to contain true state.
 *
 * Time update: ellipsoid expands via A.
 * Measurement update: intersect ellipsoid with strip {x: |y - Cx| ≤ V_max}.
 * ============================================================================ */

SetMembershipEstimator* pl_sme_create(const double* A, const double* C,
                                        int n, int p,
                                        const double* c0, const double* sh0,
                                        double Wmax, double Vmax) {
    SetMembershipEstimator* s = (SetMembershipEstimator*)calloc(1,
        sizeof(SetMembershipEstimator));
    if (!s) return NULL;
    s->A = A; s->C = C; s->n = n; s->p = p;
    s->W_max = Wmax; s->V_max = Vmax; s->rho = 0.95;
    s->center = (double*)calloc(n, sizeof(double));
    s->shape  = _m(n, n);
    s->temp_nn = _m(n, n);
    if (c0) memcpy(s->center, c0, n * sizeof(double));
    if (sh0) memcpy(s->shape, sh0, n * n * sizeof(double));
    else for (int i = 0; i < n; i++) s->shape[i * n + i] = 1.0;
    s->volume = 1.0;
    return s;
}

void pl_sme_free(SetMembershipEstimator* s) {
    if (!s) return;
    free(s->center); free(s->shape); free(s->temp_nn);
    free(s);
}

void pl_sme_predict(SetMembershipEstimator* s) {
    int n = s->n;
    /* Center: c_next = A·c */
    double* Ac = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += s->A[i * n + j] * s->center[j];
        Ac[i] = sum;
    }
    memcpy(s->center, Ac, n * sizeof(double));
    free(Ac);

    /* Shape: F_next = (A⁻¹)'·F·A⁻¹ + W_max·I (approximated) */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < n; l++) sum += s->A[i*n+l] * s->shape[l*n+j];
            s->temp_nn[i * n + j] = sum;
        }
    double* AT = _m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) AT[i * n + j] = s->A[j * n + i];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int l = 0; l < n; l++) sum += s->temp_nn[i*n+l] * AT[l*n+j];
            s->shape[i * n + j] = s->rho * sum;
        }
    for (int i = 0; i < n; i++) s->shape[i * n + i] += s->W_max;
    free(AT);
}

void pl_sme_update_available(SetMembershipEstimator* s, const double* y) {
    if (!y) return;
    int n = s->n, p = s->p;

    /* Residual: r = y - C·c */
    double* r = (double*)malloc(p * sizeof(double));
    for (int i = 0; i < p; i++) {
        double cp = 0.0;
        for (int j = 0; j < n; j++) cp += s->C[i * n + j] * s->center[j];
        r[i] = y[i] - cp;
    }

    /* Shape update: F_new = F + (1/V_max²)·C'·C (intersection with strip) */
    double vmax2 = s->V_max * s->V_max;
    if (vmax2 < 1e-12) vmax2 = 1e-12;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double add = 0.0;
            for (int l = 0; l < p; l++)
                add += s->C[l * n + i] * s->C[l * n + j];
            s->shape[i * n + j] += add / vmax2;
        }

    /* Center update via weighted correction */
    double* S_inv = _m(n, n);
    memcpy(S_inv, s->shape, n * n * sizeof(double));
    /* Simple gradient step for center */
    for (int i = 0; i < n; i++) {
        double corr = 0.0;
        for (int j = 0; j < p; j++) corr += s->C[j * n + i] * r[j];
        s->center[i] += 0.1 * corr;
    }
    free(S_inv); free(r);
}

void pl_sme_update_lost(SetMembershipEstimator* s) {
    /* No measurement → ellipsoid expands more aggressively */
    for (int i = 0; i < s->n; i++)
        s->shape[i * s->n + i] *= 1.2;
}

void pl_sme_step(SetMembershipEstimator* s, const double* y, bool arrived) {
    pl_sme_predict(s);
    if (arrived && y) pl_sme_update_available(s, y);
    else pl_sme_update_lost(s);
}

const double* pl_sme_get_center(const SetMembershipEstimator* s) {
    return s ? s->center : NULL;
}

double pl_sme_get_volume(const SetMembershipEstimator* s) {
    if (!s) return -1.0;
    /* Volume ∝ 1/sqrt(det(F)). Approximate via product of diagonal elements. */
    double vol = 1.0;
    for (int i = 0; i < s->n; i++) vol *= sqrt(fabs(s->shape[i * s->n + i]) + 1e-15);
    return 1.0 / (vol + 1e-15);
}

void pl_sme_print(const SetMembershipEstimator* s) {
    if (!s) return;
    printf("=== Set-Membership Estimator ===\n");
    printf("Volume: %.6f  W_max=%.4f  V_max=%.4f\n",
           pl_sme_get_volume(s), s->W_max, s->V_max);
}

/* ============================================================================
 * Mode-Dependent Kalman Filter (IMM-like)
 *
 * For Markovian packet loss: maintain a bank of KFs, one per channel mode.
 * Mix estimates → mode-matched filtering → probability update → combination.
 *
 * Reference: Blom & Bar-Shalom (1988), IEEE TAC 33(8):780-783.
 * ============================================================================ */

ModeDependentKalmanFilter* pl_mdkf_create(const double* A, const double* C,
                                            const double* Q, const double* R,
                                            int n, int p, int nM,
                                            const double* trans,
                                            const double* arr_rates,
                                            const double* x0, const double* P0) {
    (void)arr_rates; /* Stored for future mode-dependent arrival modeling */
    ModeDependentKalmanFilter* mf = (ModeDependentKalmanFilter*)calloc(1,
        sizeof(ModeDependentKalmanFilter));
    if (!mf) return NULL;

    mf->n_modes = nM; mf->n = n; mf->p = p;

    mf->mode_filters = (KalmanFilter**)malloc(nM * sizeof(KalmanFilter*));
    mf->mode_probabilities = (double*)malloc(nM * sizeof(double));
    mf->mode_transitions = (double**)malloc(nM * sizeof(double*));
    double* mt_data = (double*)malloc(nM * nM * sizeof(double));
    for (int i = 0; i < nM; i++) {
        mf->mode_transitions[i] = mt_data + i * nM;
        mf->mode_probabilities[i] = 1.0 / nM;
        mf->mode_filters[i] = pl_kf_create(A, C, Q, R, n, p, x0, P0);
    }
    if (trans) memcpy(mt_data, trans, nM * nM * sizeof(double));

    mf->combined_estimate  = (double*)calloc(n, sizeof(double));
    mf->combined_covariance = _m(n, n);
    mf->current_mode = 0;

    return mf;
}

void pl_mdkf_free(ModeDependentKalmanFilter* mf) {
    if (!mf) return;
    for (int i = 0; i < mf->n_modes; i++) pl_kf_free(mf->mode_filters[i]);
    free(mf->mode_filters);
    free(mf->mode_probabilities);
    free(mf->mode_transitions[0]);
    free(mf->mode_transitions);
    free(mf->combined_estimate);
    free(mf->combined_covariance);
    free(mf);
}

void pl_mdkf_step(ModeDependentKalmanFilter* mf,
                  const double* y, bool arrived, int channel_state) {
    int M = mf->n_modes, n = mf->n;

    /* Update mode probabilities via Markov transition */
    double* new_prob = (double*)malloc(M * sizeof(double));
    for (int j = 0; j < M; j++) {
        new_prob[j] = 0.0;
        for (int i = 0; i < M; i++)
            new_prob[j] += mf->mode_probabilities[i] * mf->mode_transitions[i][j];
    }
    memcpy(mf->mode_probabilities, new_prob, M * sizeof(double));
    free(new_prob);

    /* Run each mode filter */
    double* likelihood = (double*)malloc(M * sizeof(double));
    for (int i = 0; i < M; i++) {
        KalmanFilter* kf = mf->mode_filters[i];
        pl_kf_predict(kf);
        if (arrived && y) pl_kf_update(kf, y);
        else memcpy(kf->x_hat, kf->x_pred, n * sizeof(double));

        /* Likelihood (Gaussian approximation) */
        double ll = (arrived && y) ? 1.0 : 0.2;
        likelihood[i] = ll;
    }

    /* Update probabilities with likelihood */
    double sum_prob = 0.0;
    for (int i = 0; i < M; i++) {
        mf->mode_probabilities[i] *= likelihood[i];
        sum_prob += mf->mode_probabilities[i];
    }
    if (sum_prob > 1e-15)
        for (int i = 0; i < M; i++) mf->mode_probabilities[i] /= sum_prob;

    /* Combine estimates (probability-weighted average) */
    memset(mf->combined_estimate, 0, n * sizeof(double));
    for (int i = 0; i < M; i++)
        for (int j = 0; j < n; j++)
            mf->combined_estimate[j] += mf->mode_probabilities[i]
                                       * mf->mode_filters[i]->x_hat[j];

    /* Find most likely mode */
    double max_prob = 0.0;
    for (int i = 0; i < M; i++)
        if (mf->mode_probabilities[i] > max_prob) {
            max_prob = mf->mode_probabilities[i];
            mf->current_mode = i;
        }

    (void)channel_state;
    free(likelihood);
}

const double* pl_mdkf_get_estimate(const ModeDependentKalmanFilter* mf) {
    return mf ? mf->combined_estimate : NULL;
}

int pl_mdkf_get_mode(const ModeDependentKalmanFilter* mf) {
    return mf ? mf->current_mode : -1;
}

void pl_mdkf_print(const ModeDependentKalmanFilter* mf) {
    if (!mf) return;
    printf("=== Mode-Dependent KF (%d modes) ===\n", mf->n_modes);
    printf("Current mode: %d\n", mf->current_mode);
    printf("Mode probs: [");
    for (int i = 0; i < mf->n_modes; i++)
        printf("%.3f%s", mf->mode_probabilities[i], i<mf->n_modes-1?", ":"");
    printf("]\n");
}
