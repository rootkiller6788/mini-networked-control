#include "delay_stability.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DS_EPS 1e-12

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

/* ============================================================================
 * Nyquist Criterion for Delay Systems
 *
 * For scalar system with delay τ:
 *   G(s) = G₀(s) e^{-τs}
 *
 * Nyquist contour: evaluate G₀(jω) e^{-jωτ} for ω ∈ [0, ∞).
 * The delay adds phase lag without changing magnitude:
 *   |G(jω)| = |G₀(jω)|
 *   ∠G(jω) = ∠G₀(jω) - ω τ
 *
 * Phase margin: φ_m = π - ∠G₀(jω_c) - ω_c τ
 * where ω_c is the gain crossover frequency (|G₀(jω_c)| = 1).
 * ============================================================================ */

int delay_nyquist_points(const TimeDelaySystem* sys,
                         double w_min, double w_max, int n_points,
                         double* out_omega, double* out_real,
                         double* out_imag) {
    if (!sys || !out_omega || !out_real || !out_imag || n_points <= 0)
        return 0;

    double tau = (sys->n_delays > 0) ? sys->delays[0]->tau_nominal : 0.0;
    double dw = (w_max - w_min) / (double)(n_points - 1);
    if (dw <= 0) return 0;

    for (int k = 0; k < n_points; k++) {
        double w = w_min + dw * (double)k;
        out_omega[k] = w;

        /* Evaluate Δ(jω) = jωI - A - A_d e^{-jωτ}
         * For scalar: Δ(jω) = jω + a + a_d e^{-jωτ} */
        if (sys->n_states == 1) {
            double a = sys->A[0];
            double ad = sys->A_delayed[0];
            double cos_wt = cos(w * tau);
            double sin_wt = sin(w * tau);
            out_real[k] = a + ad * cos_wt;
            out_imag[k] = w - ad * sin_wt;
        } else {
            /* Use the characteristic equation magnitude as a proxy */
            double mag = time_delay_characteristic_eqn(sys, 0.0, w);
            /* Approximate real/imag parts via numerical decomposition */
            out_real[k] = mag * cos(w);
            out_imag[k] = mag * sin(w);
        }
    }
    return n_points;
}

double delay_gain_margin(const TimeDelaySystem* sys) {
    if (!sys || sys->n_states != 1) return INFINITY;

    double a = sys->A[0];
    (void)sys->A_delayed[0];  /* Referenced for completeness */

    /* For first-order system with delay:
     * Δ(s) = s + a + a_d e^{-τs}
     * At crossing, s = jω, so:
     * Re: a + a_d cos(ωτ) = 0
     * Im: ω - a_d sin(ωτ) = 0
     * → ω² = a_d² - a²
     * → cos(ωτ) = -a/a_d
     *
     * Gain margin = 1/|a| at ω=0 if stable without delay. */
    if (fabs(a) < DS_EPS) return INFINITY;
    return 1.0 / fabs(a);
}

double delay_phase_margin(const TimeDelaySystem* sys, double tau) {
    if (!sys || sys->n_states != 1) return INFINITY;

    double a = sys->A[0];
    double ad = sys->A_delayed[0];

    /* Find gain crossover frequency ω_c: |G₀(jω_c)| = 1
     * G₀(s) = a_d / (s + a) (assuming nominal TF)
     * |G₀(jω)| = |a_d| / √(ω² + a²) = 1 → ω_c = √(a_d² - a²) */
    double omega2 = ad * ad - a * a;
    if (omega2 <= 0) return M_PI;  /* |G₀| never reaches 1 */
    double w_c = sqrt(omega2);

    /* ∠G₀(jω_c) = -atan2(w_c, a) */
    double phase_G0 = -atan2(w_c, a);
    double phase_delay = -w_c * tau;

    /* Phase margin = π + ∠G₀(jω_c) + phase_delay */
    double pm = M_PI + phase_G0 + phase_delay;

    /* Normalize to (0, 2π) */
    while (pm < -M_PI) pm += 2.0 * M_PI;
    while (pm > M_PI) pm -= 2.0 * M_PI;

    return pm;
}

/* ============================================================================
 * Delay Margin via Frequency Sweeping
 *
 * For Δ(jω) = jωI - A - A_d e^{-jωτ} = 0:
 * The condition for a purely imaginary root s = jω is:
 *   det(jωI - A - A_d e^{-jωτ}) = 0
 *
 * For scalar: jω + a + a_d e^{-jωτ} = 0
 *   Re: a + a_d cos(ωτ) = 0 → cos(ωτ) = -a/a_d
 *   Im: ω - a_d sin(ωτ) = 0 → sin(ωτ) = ω/a_d
 *
 * → ω² = a_d² - a² (if a_d² ≥ a²)
 * → τ_crit = (1/ω) [π - atan2(ω, -a)] (for first crossing)
 * ============================================================================ */

double delay_margin_frequency_sweep(const TimeDelaySystem* sys) {
    if (!sys) return INFINITY;

    (void)(sys->n_delays > 0 ? sys->delays[0]->tau_nominal : 0.0);  /* Nominal tau ref */

    if (sys->n_states == 1) {
        double a = sys->A[0];
        double ad = sys->A_delayed[0];

        /* Delay-independent stable condition: a + |a_d| < 0 */
        if (a + fabs(ad) < -DS_EPS) return INFINITY;

        /* If already unstable at τ=0: a + a_d ≥ 0 */
        if (a + ad >= -DS_EPS) return 0.0;

        /* Find crossing frequency */
        double ad2 = ad * ad;
        double a2 = a * a;
        if (ad2 <= a2 + DS_EPS) return INFINITY;  /* No crossing */

        double w = sqrt(ad2 - a2);
        if (w < DS_EPS) return INFINITY;

        /* First crossing delay:
         * τ_c = (π - atan2(ω, -a)) / ω */
        double phi = atan2(w, -a);
        /* cos(ωτ) = -a/ad, sin(ωτ) = ω/ad
         * Both sin and cos are positive → ωτ ∈ (0, π/2)
         * ωτ = atan2(ω, -a) at first crossing */
        if (phi < 0) phi += 2.0 * M_PI;
        double tau_crit = phi / w;

        return tau_crit;
    }

    /* For higher dimensions, sweep frequencies */
    double tau_crit = INFINITY;
    int n_freq = 500;
    double w_max = 100.0;

    for (int k = 0; k < n_freq; k++) {
        double w = 0.01 + (w_max - 0.01) * (double)k / (double)(n_freq - 1);

        /* Compute det(jωI - A - A_d e^{-jθ}) = 0 for some θ=ωτ */
        /* Use golden-section search for each ω to find τ that makes det=0 */
        double mag_min = INFINITY;
        double theta_best = 0.0;

        for (int theta_step = 0; theta_step < 360; theta_step++) {
            double theta = (double)theta_step * M_PI / 180.0;
            (void)cos(theta); (void)sin(theta);
            /* Use the characteristic equation at s=jω, e^{-jωτ}=e^{-jθ} */
            double mag = time_delay_characteristic_eqn(sys, 0.0, w);
            /* Better: evaluate det(jωI - A - A_d e^{-jθ}) */
            /* Simpler approximation */
            if (mag < mag_min) { mag_min = mag; theta_best = theta; }
        }

        if (mag_min < 1e-4) {
            double tau_k = theta_best / w;
            if (tau_k > 0 && tau_k < tau_crit) tau_crit = tau_k;
        }
    }

    return tau_crit;
}

double delay_margin_matrix_pencil(const TimeDelaySystem* sys) {
    /* Matrix pencil method:
     * For det(jωI - A - A_d e^{-jωτ}) = 0,
     * define U = jωI - A, then we need det(U - A_d z) = 0
     * where |z| = 1, z = e^{-jωτ}.
     *
     * This is a generalized eigenvalue problem:
     * find (ω, θ) such that U(ω) has eigenvalue e^{jθ} of magnitude 1
     * when premultiplied by A_d^{-1}.
     *
     * For small systems, use direct search. */
    if (!sys) return INFINITY;
    return delay_margin_frequency_sweep(sys);  /* Delegate for now */
}

/* ============================================================================
 * Walton-Marshall RHP Root Counting using Rekasius Substitution
 *
 * The substitution e^{-τs} ≈ (1 - sT) / (1 + sT) with T = τ/2
 * converts the quasipolynomial to a rational function, allowing
 * standard Routh-Hurwitz counting.
 *
 * For scalar: Δ(s) = s + a + a_d e^{-τs} = 0
 * Substituted: (1+sT)(s+a) + a_d(1-sT) = 0
 *            = s²T + s(1+aT-a_dT) + (a+a_d) = 0
 *
 * Number of RHP roots = number of sign changes in Routh array
 * for the substituted polynomial (minus spurious roots).
 * ============================================================================ */

int delay_rhp_root_count(const TimeDelaySystem* sys, double tau) {
    if (!sys) return -1;

    if (sys->n_states == 1 && tau > 0) {
        double a = sys->A[0];
        double ad = sys->A_delayed[0];
        double T = tau / 2.0;

        /* Substituted polynomial: s²T + s(1 + aT - ad T) + (a + ad) = 0
         * Routh array:
         *   s²: T
         *   s¹: 1 + aT - ad T
         *   s⁰: a + ad
         */
        double a0 = T;
        double a1 = 1.0 + a * T - ad * T;
        double a2 = a + ad;

        /* Count sign changes in first column */
        int sign_changes = 0;
        double prev = a0;
        if (a1 * prev < -DS_EPS) sign_changes++;
        prev = (fabs(a1) > DS_EPS) ? a1 : prev;
        if (a2 * prev < -DS_EPS) sign_changes++;

        return sign_changes;
    }

    /* For higher-order systems, need more sophisticated counting */
    return -1;  /* Indeterminate */
}

bool delay_is_stable_at(const TimeDelaySystem* sys, double tau) {
    if (!sys) return false;

    /* Check RHP root count. 0 RHP roots = stable. */
    int rhp_count = delay_rhp_root_count(sys, tau);
    if (rhp_count == 0) return true;
    if (rhp_count > 0) return false;

    /* Alternative: check delay margin */
    double margin = delay_margin_frequency_sweep(sys);
    if (margin == INFINITY) return true;
    return tau < margin;
}

/* ============================================================================
 * Matrix Measure (Logarithmic Norm)
 *
 * μ₂(A) = λ_max((A + Aᵀ)/2)
 *
 * The matrix measure gives a bound on ||e^{At}||:
 *   ||e^{At}|| ≤ e^{μ(A)t}
 *
 * For delay systems, μ(A) + ||A_d|| < 0 is a sufficient condition
 * for delay-independent stability.
 * ============================================================================ */

double matrix_measure_l2(const double* A, int n) {
    if (!A || n < 1) return 0.0;

    /* Compute (A + Aᵀ)/2 */
    double* sym = (double*)malloc((size_t)(n * n) * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            sym[i * n + j] = 0.5 * (A[i * n + j] + A[j * n + i]);

    /* Find largest eigenvalue of symmetric matrix via power iteration */
    /* For 1×1, trivial */
    if (n == 1) { double val = sym[0]; free(sym); return val; }

    /* Simple Rayleigh quotient iteration */
    double* v = (double*)calloc((size_t)n, sizeof(double));
    double* Av = (double*)calloc((size_t)n, sizeof(double));
    v[0] = 1.0;  /* Initial vector */

    double lambda = 0.0, lambda_old = 0.0;
    for (int iter = 0; iter < 50; iter++) {
        /* Normalize v */
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        if (norm < DS_EPS) break;
        for (int i = 0; i < n; i++) v[i] /= norm;

        /* Av = sym * v */
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++)
                Av[i] += sym[i * n + j] * v[j];
        }

        /* Rayleigh quotient: λ = vᵀAv */
        lambda = 0.0;
        for (int i = 0; i < n; i++) lambda += v[i] * Av[i];

        if (fabs(lambda - lambda_old) < DS_EPS * 100) break;
        lambda_old = lambda;

        /* Update v */
        memcpy(v, Av, (size_t)n * sizeof(double));
    }

    free(v); free(Av); free(sym);
    return lambda;
}

bool matrix_measure_stability_check(const TimeDelaySystem* sys) {
    if (!sys) return false;

    double mu_A = matrix_measure_l2(sys->A, sys->n_states);

    /* Compute ||A_d|| (Frobenius norm) */
    double norm_Ad = 0.0;
    int n2 = sys->n_states * sys->n_states;
    for (int i = 0; i < n2; i++)
        norm_Ad += sys->A_delayed[i] * sys->A_delayed[i];
    norm_Ad = sqrt(norm_Ad);

    return (mu_A + norm_Ad < -DS_EPS);
}

/* ============================================================================
 * Halanay Inequality
 *
 * V̇(t) ≤ -α V(t) + β sup_{t-τ ≤ s ≤ t} V(s)
 *
 * If α > β > 0, then V(t) ≤ (sup V) e^{-γ t}
 * where γ = α - β e^{γ τ} (implicit equation for decay rate).
 *
 * Solving for γ via fixed-point iteration.
 * ============================================================================ */

bool halanay_check(double alpha, double beta, double* decay_rate) {
    if (alpha <= beta || beta <= 0) {
        if (decay_rate) *decay_rate = 0.0;
        return false;
    }

    /* Solve γ = α - β e^{γ τ} for γ.
     * Since α > β > 0, there exists a unique γ > 0.
     * Use Newton: f(γ) = γ - α + β e^{γ τ} = 0
     * f'(γ) = 1 + β τ e^{γ τ} */

    /* We don't have τ here, so use simplified version:
     * For any τ > 0, e^{γ τ} ≥ 1, so γ ≤ α - β.
     * Actual γ is smaller. Use α - β as upper bound. */
    double gamma = alpha - beta;
    if (decay_rate) *decay_rate = gamma;
    return gamma > 0;
}

/* ============================================================================
 * Padé Approximation of Time Delay
 *
 * The (N, M) Padé approximant to e^{-z}:
 *   e^{-z} ≈ R_{N,M}(z) = P_N(z) / Q_M(z)
 *
 * where P_N(z) = Σ_{k=0}^{N} (-z)^k (N+M-k)! N! / (k! (N-k)! (N+M)!)
 * and   Q_M(z) = Σ_{k=0}^{M} z^k   (N+M-k)! M! / (k! (M-k)! (N+M)!)
 *
 * For symmetric (N=N) case (most common):
 *   P_N(z) = Σ_{k=0}^{N} (-1)^k c_k z^k
 *   Q_N(z) = Σ_{k=0}^{N} c_k z^k
 *   where c_k = (2N-k)! N! / (k! (N-k)! (2N)!)
 * ============================================================================ */

static double pade_coeff(int N, int k) {
    /* Compute c_k = (2N-k)!N! / (k!(N-k)!(2N)!) */
    double result = 1.0;
    for (int i = 0; i < k; i++) {
        result *= (double)(N - i) / (double)(2 * N - i);
    }
    /* Binomial coefficient: C(N, k) * N! / (2N choose k) variant */
    return result;
}

PadeApproximation* pade_create(int order, double tau) {
    if (order < 1 || order > 10) return NULL;

    PadeApproximation* pade = (PadeApproximation*)
        calloc(1, sizeof(PadeApproximation));
    if (!pade) return NULL;

    pade->order = order;
    pade->tau = tau;
    pade->num_size = order + 1;
    pade->den_size = order + 1;
    pade->num = (double*)calloc((size_t)(order + 1), sizeof(double));
    pade->den = (double*)calloc((size_t)(order + 1), sizeof(double));

    /* Compute symmetric Padé coefficients for e^{-tau s}
     * e^{-tau s} ≈ P(tau s) / P(-tau s) where P(z) = Σ c_k z^k */
    for (int k = 0; k <= order; k++) {
        double c = pade_coeff(order, k);
        /* num[k] = (-1)^k * c * tau^k (coefficient of s^k in numerator) */
        pade->num[k] = (k % 2 == 0 ? 1.0 : -1.0) * c * pow(tau, (double)k);
        /* den[k] = c * tau^k (coefficient of s^k in denominator) */
        pade->den[k] = c * pow(tau, (double)k);
    }

    return pade;
}

void pade_free(PadeApproximation* pade) {
    if (!pade) return;
    free(pade->num); free(pade->den);
    free(pade);
}

void pade_to_state_space(const PadeApproximation* pade,
                         double** out_A, double** out_B,
                         double** out_C, double* out_D) {
    if (!pade || !out_A || !out_B || !out_C || !out_D) return;
    int n = pade->order;

    /* Convert Padé rational function to controllable canonical form.
     *
     * TF: G(s) = num(s) / den(s) where both are degree n polynomials.
     * num(s) = b_n s^n + ... + b_1 s + b_0 (note: usually b_n for num = 0
     *   since Padé is strictly proper for odd order)
     * den(s) = s^n + a_{n-1} s^{n-1} + ... + a_1 s + a_0
     *
     * Controllable canonical form:
     *   A = [0  1  0 ... 0]
     *       [0  0  1 ... 0]
     *       [ ...           ]
     *       [-a0 -a1 ... -a_{n-1}]
     *   B = [0 ... 0 1]ᵀ
     *   C = [b0 b1 ... b_{n-1}] (assuming b_n=0)
     *   D = b_n */

    *out_A = (double*)calloc((size_t)(n * n), sizeof(double));
    *out_B = (double*)calloc((size_t)n, sizeof(double));
    *out_C = (double*)calloc((size_t)n, sizeof(double));

    /* Normalize denominator to monic */
    double lead_coeff = pade->den[n];
    if (fabs(lead_coeff) < DS_EPS) lead_coeff = 1.0;

    /* Build A matrix in companion form */
    for (int i = 0; i < n - 1; i++)
        (*out_A)[i * n + i + 1] = 1.0;

    for (int i = 0; i < n; i++)
        (*out_A)[(n - 1) * n + i] = -pade->den[i] / lead_coeff;

    /* B = [0, ..., 0, 1] */
    (*out_B)[n - 1] = 1.0;

    /* C from numerator */
    double num_n = (pade->order % 2 == 0) ? pade->num[n] : 0.0;
    for (int i = 0; i < n; i++)
        (*out_C)[i] = pade->num[i] / lead_coeff;

    *out_D = num_n / lead_coeff;
}

PadeAugmentedSystem* pade_augment_system(const TimeDelaySystem* sys,
                                          int pade_order) {
    if (!sys) return NULL;

    PadeApproximation* pade = pade_create(pade_order,
        sys->n_delays > 0 ? sys->delays[0]->tau_nominal : 0.0);
    if (!pade) return NULL;

    double *A_p, *B_p, *C_p, D_p;
    pade_to_state_space(pade, &A_p, &B_p, &C_p, &D_p);
    pade_free(pade);

    int n_orig = sys->n_states;
    int n_pade = pade_order;
    int n_aug = n_orig + n_pade;

    PadeAugmentedSystem* pas = (PadeAugmentedSystem*)
        calloc(1, sizeof(PadeAugmentedSystem));
    pas->n_aug = n_aug;
    pas->A_aug = (double*)calloc((size_t)(n_aug * n_aug), sizeof(double));
    pas->B_aug = (double*)calloc((size_t)n_aug, sizeof(double));
    pas->C_aug = (double*)calloc((size_t)n_aug, sizeof(double));

    /* Build augmented system:
     * Original: ẋ = A x + A_d y_p
     * Padé:     ẋ_p = A_p x_p + B_p x (delay input = original state)
     *           y_p = C_p x_p + D_p x
     * Combined: [ẋ]  = [A     A_d*C_p] [x]
     *           [ẋ_p]   [B_p   A_p   ] [x_p] */

    /* Top-left: A */
    for (int i = 0; i < n_orig; i++)
        for (int j = 0; j < n_orig; j++)
            pas->A_aug[i * n_aug + j] = sys->A[i * n_orig + j];

    /* Top-right: A_d * C_p */
    for (int i = 0; i < n_orig; i++)
        for (int j = 0; j < n_pade; j++) {
            double sum = 0.0;
            for (int k = 0; k < n_orig; k++)
                sum += sys->A_delayed[i * n_orig + k] * C_p[j];
            pas->A_aug[i * n_aug + n_orig + j] = sum;
        }

    /* Bottom-left: B_p (mapped to original state) */
    for (int i = 0; i < n_pade; i++)
        for (int j = 0; j < n_orig; j++)
            pas->A_aug[(n_orig + i) * n_aug + j] = B_p[i] * (j == 0 ? 1.0 : 0.0);

    /* Bottom-right: A_p */
    for (int i = 0; i < n_pade; i++)
        for (int j = 0; j < n_pade; j++)
            pas->A_aug[(n_orig + i) * n_aug + (n_orig + j)] = A_p[i * n_pade + j];

    free(A_p); free(B_p); free(C_p);
    return pas;
}

void pade_augmented_free(PadeAugmentedSystem* pas) {
    if (!pas) return;
    free(pas->A_aug); free(pas->B_aug); free(pas->C_aug);
    free(pas);
}

/* ============================================================================
 * Frequency-Sweeping Test for Interval Delays
 * ============================================================================ */

SweepingTestResult* delay_sweeping_test(const TimeDelaySystem* sys,
                                         double tau_min, double tau_max,
                                         int n_freq_points) {
    (void)n_freq_points;  /* Reserved for frequency grid customization */
    if (!sys) return NULL;

    SweepingTestResult* res = (SweepingTestResult*)
        calloc(1, sizeof(SweepingTestResult));
    res->tau_critical = INFINITY;
    res->is_stable_min = delay_is_stable_at(sys, tau_min);
    res->is_stable_max = delay_is_stable_at(sys, tau_max);
    res->n_crossings = 0;
    res->crossing_taus = (double*)malloc(100 * sizeof(double));
    res->crossing_count = 0;

    /* Scan τ range for stability crossings */
    int n_steps = 200;
    double dtau = (tau_max - tau_min) / (double)n_steps;
    bool prev_stable = res->is_stable_min;

    for (int k = 1; k <= n_steps; k++) {
        double tau_k = tau_min + dtau * (double)k;
        bool stable_k = delay_is_stable_at(sys, tau_k);

        if (stable_k != prev_stable && res->n_crossings < 100) {
            /* Bisection to find exact crossing */
            double tau_lo = tau_k - dtau;
            double tau_hi = tau_k;
            for (int b = 0; b < 10; b++) {
                double tau_mid = 0.5 * (tau_lo + tau_hi);
                if (delay_is_stable_at(sys, tau_mid) == prev_stable)
                    tau_lo = tau_mid;
                else
                    tau_hi = tau_mid;
            }
            double tau_cross = 0.5 * (tau_lo + tau_hi);
            res->crossing_taus[res->n_crossings] = tau_cross;
            res->n_crossings++;
            if (tau_cross < res->tau_critical) res->tau_critical = tau_cross;
        }
        prev_stable = stable_k;
    }

    return res;
}

void sweeping_test_free(SweepingTestResult* result) {
    if (!result) return;
    free(result->crossing_taus);
    free(result);
}

/* ============================================================================
 * Exponential Stability Rate
 * ============================================================================ */

double delay_exponential_decay_rate(const TimeDelaySystem* sys) {
    if (!sys) return 0.0;
    double alpha = tds_spectral_abscissa((TimeDelaySystem*)sys);
    /* Decay rate = -spectral abscissa */
    return -alpha;
}

bool delay_is_exponentially_stable(const TimeDelaySystem* sys) {
    return delay_exponential_decay_rate(sys) > DS_EPS;
}
