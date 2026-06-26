#include "time_delay_system.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * Internal Constants & Helpers
 * ============================================================================ */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TDS_EPS 1e-12
#define TDS_MAX_ROOTS_DEFAULT 64
#define TDS_HISTORY_POINTS_DEFAULT 512

static double tds_norm2(const double* v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}

__attribute__((unused))
static double tds_dot(const double* a, const double* b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

/* Complex exponential: e^{σ + jω} */
static void complex_exp(double sigma, double omega,
                        double* out_real, double* out_imag) {
    double mag = exp(sigma);
    *out_real = mag * cos(omega);
    *out_imag = mag * sin(omega);
}

/* Matrix-vector: y = A * x, A is n×n (row-major) */
__attribute__((unused))
static void mat_vec_mul(const double* A, const double* x, int n,
                         double* y) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i * n + j] * x[j];
        }
    }
}

/* ============================================================================
 * L1 — Core Definition: DelayDescriptor
 * ============================================================================ */

DelayDescriptor* delay_create_constant(double tau) {
    DelayDescriptor* d = (DelayDescriptor*)calloc(1, sizeof(DelayDescriptor));
    if (!d) return NULL;
    d->type = DELAY_CONSTANT;
    d->tau_nominal = tau;
    d->tau_min = tau;
    d->tau_max = tau;
    d->tau_variance = 0.0;
    d->derivative_bound = 0.0;   /* dτ/dt = 0 for constant */
    d->is_bounded = true;
    d->kernel = NULL;
    d->kernel_size = 0;
    return d;
}

DelayDescriptor* delay_create_time_varying(double tau_min, double tau_max,
                                            double deriv_bound) {
    DelayDescriptor* d = (DelayDescriptor*)calloc(1, sizeof(DelayDescriptor));
    if (!d) return NULL;
    d->type = DELAY_TIME_VARYING;
    d->tau_nominal = 0.5 * (tau_min + tau_max);
    d->tau_min = tau_min;
    d->tau_max = tau_max;
    d->tau_variance = (tau_max - tau_min) * (tau_max - tau_min) / 12.0;
    d->derivative_bound = deriv_bound;
    d->is_bounded = true;
    d->kernel = NULL;
    d->kernel_size = 0;
    return d;
}

DelayDescriptor* delay_create_stochastic(double tau_mean, double tau_var) {
    DelayDescriptor* d = (DelayDescriptor*)calloc(1, sizeof(DelayDescriptor));
    if (!d) return NULL;
    d->type = DELAY_STOCHASTIC;
    d->tau_nominal = tau_mean;
    d->tau_min = tau_mean - 3.0 * sqrt(tau_var);  /* 3-sigma bound */
    d->tau_max = tau_mean + 3.0 * sqrt(tau_var);
    if (d->tau_min < 0.0) d->tau_min = 0.0;
    d->tau_variance = tau_var;
    d->derivative_bound = 0.0;
    d->is_bounded = false;   /* Stochastic can exceed bounds */
    d->kernel = NULL;
    d->kernel_size = 0;
    return d;
}

void delay_free(DelayDescriptor* delay) {
    if (!delay) return;
    free(delay->kernel);
    free(delay);
}

void delay_print(const DelayDescriptor* delay) {
    if (!delay) { printf("DelayDescriptor: NULL\n"); return; }
    printf("DelayDescriptor: type=%d tau_nom=%.4f tau∈[%.4f, %.4f] "
           "var=%.6f |dτ/dt|≤%.4f bounded=%d\n",
           delay->type, delay->tau_nominal,
           delay->tau_min, delay->tau_max,
           delay->tau_variance, delay->derivative_bound,
           delay->is_bounded);
}

/* ============================================================================
 * TimeDelaySystem — Lifecycle
 * ============================================================================ */

TimeDelaySystem* tds_create(const char* name, int n_states,
                             int n_inputs, int n_outputs) {
    TimeDelaySystem* sys = (TimeDelaySystem*)calloc(1, sizeof(TimeDelaySystem));
    if (!sys) return NULL;

    sys->name = name ? strdup(name) : strdup("unnamed");
    sys->n_states = n_states;
    sys->n_inputs = n_inputs;
    sys->n_outputs = n_outputs;

    int n2 = n_states * n_states;
    sys->A = (double*)calloc(n2, sizeof(double));
    sys->A_delayed = (double*)calloc(n2, sizeof(double));
    sys->B = (double*)calloc(n_states * n_inputs, sizeof(double));
    sys->C = (double*)calloc(n_outputs * n_states, sizeof(double));

    sys->n_delays = 0;
    sys->delays = NULL;
    sys->f_rhs = NULL;
    sys->history = NULL;
    sys->rhs_params = NULL;
    sys->n_params = 0;
    sys->dde_type = DDE_RETARDED;

    sys->n_roots = 0;
    sys->roots_real = (double*)calloc(TDS_MAX_ROOTS_DEFAULT, sizeof(double));
    sys->roots_imag = (double*)calloc(TDS_MAX_ROOTS_DEFAULT, sizeof(double));
    sys->stability_class = DELAY_STABLE;
    sys->delay_margin = INFINITY;
    sys->spectral_abscissa = -INFINITY;

    sys->current_state = (double*)calloc(n_states, sizeof(double));
    sys->delayed_state = (double*)calloc(n_states, sizeof(double));
    sys->history_buffer = (double*)calloc(
        TDS_HISTORY_POINTS_DEFAULT * n_states, sizeof(double));
    sys->history_points = TDS_HISTORY_POINTS_DEFAULT;
    sys->t_current = 0.0;

    sys->settling_time = 0.0;
    sys->overshoot_percent = 0.0;
    sys->ise = 0.0;
    sys->iae = 0.0;

    return sys;
}

void tds_set_linear_model(TimeDelaySystem* sys,
                          const double* A, const double* A_d,
                          const double* B, const double* C) {
    if (!sys) return;
    int n2 = sys->n_states * sys->n_states;
    if (A)  memcpy(sys->A, A, n2 * sizeof(double));
    if (A_d) memcpy(sys->A_delayed, A_d, n2 * sizeof(double));
    if (B)  memcpy(sys->B, B, sys->n_states * sys->n_inputs * sizeof(double));
    if (C)  memcpy(sys->C, C, sys->n_outputs * sys->n_states * sizeof(double));
}

void tds_add_delay(TimeDelaySystem* sys, DelayType type,
                   double tau_nominal, double tau_min, double tau_max) {
    if (!sys) return;
    sys->n_delays++;
    sys->delays = (DelayDescriptor**)realloc(
        sys->delays, (size_t)sys->n_delays * sizeof(DelayDescriptor*));
    DelayDescriptor* dd = NULL;
    switch (type) {
        case DELAY_CONSTANT:
            dd = delay_create_constant(tau_nominal); break;
        case DELAY_TIME_VARYING:
            dd = delay_create_time_varying(tau_min, tau_max, 0.5); break;
        case DELAY_STOCHASTIC:
            dd = delay_create_stochastic(tau_nominal,
                (tau_max - tau_min) * (tau_max - tau_min) / 36.0); break;
        default:
            dd = delay_create_constant(tau_nominal); break;
    }
    sys->delays[sys->n_delays - 1] = dd;
}

void tds_set_history(TimeDelaySystem* sys, HistoryFunc phi) {
    if (!sys) return;
    sys->history = phi;
}

void tds_set_nonlinear_rhs(TimeDelaySystem* sys, DDERHSFunc f,
                            const double* params, int n_params) {
    if (!sys) return;
    sys->f_rhs = f;
    sys->n_params = n_params;
    free(sys->rhs_params);
    sys->rhs_params = NULL;
    if (n_params > 0 && params) {
        sys->rhs_params = (double*)malloc((size_t)n_params * sizeof(double));
        memcpy(sys->rhs_params, params, (size_t)n_params * sizeof(double));
    }
}

void tds_set_dde_type(TimeDelaySystem* sys, DDEType dtype) {
    if (sys) sys->dde_type = dtype;
}

void tds_free(TimeDelaySystem* sys) {
    if (!sys) return;
    free(sys->name);
    free(sys->A);
    free(sys->A_delayed);
    free(sys->B);
    free(sys->C);
    if (sys->delays) {
        for (int i = 0; i < sys->n_delays; i++) delay_free(sys->delays[i]);
        free(sys->delays);
    }
    free(sys->rhs_params);
    free(sys->roots_real);
    free(sys->roots_imag);
    free(sys->current_state);
    free(sys->delayed_state);
    free(sys->history_buffer);
    free(sys);
}

void tds_print(const TimeDelaySystem* sys) {
    if (!sys) { printf("TimeDelaySystem: NULL\n"); return; }
    printf("=== TimeDelaySystem: %s ===\n", sys->name);
    printf("  States: %d  Inputs: %d  Outputs: %d\n",
           sys->n_states, sys->n_inputs, sys->n_outputs);
    printf("  DDE Type: %d  Delays: %d\n", sys->dde_type, sys->n_delays);
    if (sys->n_delays > 0 && sys->delays) {
        printf("  Delay[0]: τ=%.4f type=%d\n",
               sys->delays[0]->tau_nominal, sys->delays[0]->type);
    }
    printf("  Stability class: %d  Delay margin: %.6f\n",
           sys->stability_class, sys->delay_margin);
    printf("  Spectral abscissa: %.6f\n", sys->spectral_abscissa);
    printf("  t=%.4f ||x||=%.6f\n",
           sys->t_current, tds_norm2(sys->current_state, sys->n_states));
}

/* ============================================================================
 * L3 — Characteristic Quasipolynomial
 *
 * For linear DDE: ẋ(t) = A x(t) + A_d x(t-τ)
 * Characteristic equation:
 *   Δ(s) = det(sI - A - A_d e^{-τs}) = 0
 *
 * For scalar case (n=1): Δ(s) = s + a + a_d e^{-τs}
 *   Real: σ + a + a_d e^{-στ} cos(ωτ)
 *   Imag: ω - a_d e^{-στ} sin(ωτ)
 * ============================================================================ */

double time_delay_characteristic_eqn(const TimeDelaySystem* sys,
                                     double sigma, double omega) {
    if (!sys || sys->n_states < 1) return 0.0;
    double tau = (sys->n_delays > 0 && sys->delays)
                 ? sys->delays[0]->tau_nominal : 0.0;

    /* For scalar system */
    if (sys->n_states == 1) {
        double a = sys->A[0];
        double ad = sys->A_delayed[0];
        double exp_real, exp_imag;
        complex_exp(-sigma * tau, -omega * tau, &exp_real, &exp_imag);
        /* Δ(s) = s + a + ad * e^{-τs}
         * s = sigma + j*omega, e^{-τs} = e^{-τσ} (cos(ωτ) - j sin(ωτ)) */
        double real_part = sigma + a + ad * exp_real * cos(omega * tau);
        double imag_part = omega - ad * exp_real * sin(omega * tau);
        return sqrt(real_part * real_part + imag_part * imag_part);
    }

    /* For n>1: Δ(s) = det(sI - A - A_d e^{-τs})
     * We approximate by computing the determinant numerically. */
    int n = sys->n_states;
    double exp_real, exp_imag;
    complex_exp(-sigma * tau, -omega * tau, &exp_real, &exp_imag);

    /* Build M = sI - A - A_d * e^{-τs} */
    /* M stores complex matrix as real parts then imag parts interleaved...
     * Actually let's use separate arrays for real and imag parts. */
    double* Mr = (double*)calloc((size_t)(n * n), sizeof(double));
    double* Mi = (double*)calloc((size_t)(n * n), sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a_ij = sys->A[i * n + j];
            double ad_ij = sys->A_delayed[i * n + j];
            double e_cos = exp_real * cos(omega * tau);
            double e_sin = -exp_real * sin(omega * tau);

            Mr[i * n + j] = -a_ij - ad_ij * e_cos;
            Mi[i * n + j] = -ad_ij * e_sin;

            if (i == j) {
                Mr[i * n + j] += sigma;
                Mi[i * n + j] += omega;
            }
        }
    }

    /* Determinant of complex matrix — use Gaussian elimination
     * for 2×2 or compute directly. For larger matrices use LU.
     * Here we handle 2×2 explicitly for common cases. */
    double det_r = 1.0, det_i = 0.0;
    if (n == 1) {
        det_r = Mr[0]; det_i = Mi[0];
    } else if (n == 2) {
        /* det = M[0][0]*M[1][1] - M[0][1]*M[1][0] */
        double a_r = Mr[0], a_i = Mi[0];
        double b_r = Mr[1], b_i = Mi[1];
        double c_r = Mr[2], c_i = Mi[2];
        double d_r = Mr[3], d_i = Mi[3];
        /* (a+bi)(d+di) - (b+bi)(c+ci) */
        det_r = a_r * d_r - a_i * d_i - (b_r * c_r - b_i * c_i);
        det_i = a_r * d_i + a_i * d_r - (b_r * c_i + b_i * c_r);
    } else {
        /* General n: use complex LU decomposition */
        /* For now, compute magnitude via QR-like approach */
        double* tmp_r = (double*)malloc((size_t)(n * n) * sizeof(double));
        double* tmp_i = (double*)malloc((size_t)(n * n) * sizeof(double));
        memcpy(tmp_r, Mr, (size_t)(n * n) * sizeof(double));
        memcpy(tmp_i, Mi, (size_t)(n * n) * sizeof(double));

        det_r = 1.0; det_i = 0.0;
        for (int k = 0; k < n; k++) {
            /* Find pivot */
            int pivot = k;
            double pmax = tmp_r[k * n + k] * tmp_r[k * n + k]
                        + tmp_i[k * n + k] * tmp_i[k * n + k];
            for (int i = k + 1; i < n; i++) {
                double v = tmp_r[i * n + k] * tmp_r[i * n + k]
                         + tmp_i[i * n + k] * tmp_i[i * n + k];
                if (v > pmax) { pmax = v; pivot = i; }
            }
            if (pmax < TDS_EPS) { det_r = 0.0; det_i = 0.0; break; }

            if (pivot != k) {
                for (int j = 0; j < n; j++) {
                    double tr = tmp_r[k * n + j];
                    double ti = tmp_i[k * n + j];
                    tmp_r[k * n + j] = tmp_r[pivot * n + j];
                    tmp_i[k * n + j] = tmp_i[pivot * n + j];
                    tmp_r[pivot * n + j] = tr;
                    tmp_i[pivot * n + j] = ti;
                }
                det_r = -det_r; det_i = -det_i;
            }

            double piv_r = tmp_r[k * n + k];
            double piv_i = tmp_i[k * n + k];
            double dnorm = piv_r * piv_r + piv_i * piv_i;
            /* Multiply determinant by pivot */
            double ndr = det_r * piv_r - det_i * piv_i;
            double ndi = det_r * piv_i + det_i * piv_r;
            det_r = ndr; det_i = ndi;

            for (int i = k + 1; i < n; i++) {
                double mik_r = (tmp_r[i * n + k] * piv_r + tmp_i[i * n + k] * piv_i) / dnorm;
                double mik_i = (tmp_i[i * n + k] * piv_r - tmp_r[i * n + k] * piv_i) / dnorm;
                for (int j = k; j < n; j++) {
                    tmp_r[i * n + j] -= mik_r * tmp_r[k * n + j] - mik_i * tmp_i[k * n + j];
                    tmp_i[i * n + j] -= mik_r * tmp_i[k * n + j] + mik_i * tmp_r[k * n + j];
                }
            }
        }
        free(tmp_r);
        free(tmp_i);
    }

    free(Mr); free(Mi);
    return sqrt(det_r * det_r + det_i * det_i);
}

void time_delay_char_eqn_parts(const TimeDelaySystem* sys,
                                double sigma, double omega,
                                double* out_real, double* out_imag) {
    if (!out_real || !out_imag) return;
    double tau = (sys && sys->n_delays > 0 && sys->delays)
                 ? sys->delays[0]->tau_nominal : 0.0;

    if (!sys || sys->n_states == 1) {
        double a = sys ? sys->A[0] : 0.0;
        double ad = sys ? sys->A_delayed[0] : 0.0;
        double e_r = exp(-sigma * tau) * cos(omega * tau);
        double e_i = -exp(-sigma * tau) * sin(omega * tau);
        *out_real = sigma + a + ad * e_r;
        *out_imag = omega + ad * e_i;
        return;
    }

    /* For larger systems, compute via numerical determinant */
    double mag = time_delay_characteristic_eqn(sys, sigma, omega);
    /* Phase approximation via finite differences */
    double eps_d = 1e-6;
    double m_plus_s = time_delay_characteristic_eqn(sys, sigma + eps_d, omega);
    double m_plus_w = time_delay_characteristic_eqn(sys, sigma, omega + eps_d);
    (void)((m_plus_s - mag) / eps_d);
    (void)((m_plus_w - mag) / eps_d);

    /* Not exact — decompose from characteristic equation structure */
    *out_real = mag;
    *out_imag = 0.0;
}

/* ============================================================================
 * L2 — State Norm for DDE: ||x_t|| = sup_{θ∈[-τ,0]} ||x(t+θ)||
 * ============================================================================ */

double time_delay_state_norm(const TimeDelaySystem* sys) {
    if (!sys) return 0.0;

    /* Current state norm */
    double max_norm = tds_norm2(sys->current_state, sys->n_states);

    /* Scan history buffer for maximum */
    if (sys->history_buffer && sys->history_points > 0) {
        for (int i = 0; i < sys->history_points; i++) {
            double norm = tds_norm2(
                sys->history_buffer + (size_t)i * sys->n_states,
                sys->n_states);
            if (norm > max_norm) max_norm = norm;
        }
    }
    return max_norm;
}

/* ============================================================================
 * L2 — Delay Rate Condition
 * For time-varying delays, |dτ/dt| ≤ d < 1 is often required
 * for Lyapunov-based stability criteria.
 * ============================================================================ */

bool time_delay_rate_check(const DelayDescriptor* delay) {
    if (!delay) return false;
    if (delay->type == DELAY_CONSTANT) return true;  /* dτ/dt = 0 < 1 */
    return delay->derivative_bound < 1.0;
}

/* ============================================================================
 * L3 — Characteristic Root Computation
 *
 * Method: Classical spectral discretization.
 * Approximate the solution operator using Chebyshev collocation
 * on an interval [-τ, 0] and compute eigenvalues of the
 * resulting finite-dimensional matrix.
 *
 * For ẋ = A x + A_d x(t-τ):
 * The infinitesimal generator of the semigroup has the form:
 *   𝒜 φ = φ'  on [-τ, 0]
 *   𝒟(𝒜) = {φ ∈ C¹ : φ'(0) = A φ(0) + A_d φ(-τ)}
 *
 * Discretize φ' at Chebyshev points and substitute.
 * ============================================================================ */

/* Build the discretized solution operator matrix
 * for eigenvalue computation using Chebyshev differentiation. */
static void build_operator_matrix(const TimeDelaySystem* sys,
                                   double* D, int N) {
    int n = sys->n_states;
    int M = n * (N + 1);
    double tau = (sys->n_delays > 0) ? sys->delays[0]->tau_nominal : 1.0;

    memset(D, 0, (size_t)(M * M) * sizeof(double));

    /* Chebyshev differentiation matrix (scaled to [-τ, 0])
     * D_ij ≈ d/dθ evaluated at Chebyshev points
     * Scale by 2/τ to map [-1,1] → [-τ,0] */
    double scale = 2.0 / tau;

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            if (i == j) {
                double s = 0.0;
                for (int k = 0; k <= N; k++) {
                    if (k != i) {
                        double xi = cos(M_PI * i / N);
                        double xk = cos(M_PI * k / N);
                        s += 1.0 / (xi - xk + 1e-15);
                    }
                }
                D[i * (N + 1) + j] = scale * s;
            } else {
                double xi = cos(M_PI * i / N);
                double xj = cos(M_PI * j / N);
                double num = 1.0;
                for (int k = 0; k <= N; k++) {
                    if (k != i && k != j) {
                        num *= (xi - cos(M_PI * k / N))
                              / (xj - cos(M_PI * k / N));
                    }
                }
                D[i * (N + 1) + j] = scale * num / (xj - xi + 1e-15);
            }
        }
    }

    /* Embed state-space matrices: build block matrix
     * [ D ⊗ I_n  ]
     * with the boundary condition at i=0 (corresponding to θ=0):
     * Replace row 0 with (A ⊗ e_0 + A_d ⊗ e_N) where e_0, e_N
     * select the boundary values. */

    /* Put A at (row=0, col=0) and A_d at (row=0, col=N) */
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            /* A at θ=0 */
            D[r * M + c] = sys->A[r * n + c];
            /* A_d at θ=-τ (last Chebyshev point) */
            D[r * M + N * n + c] = sys->A_delayed[r * n + c];
        }
    }

    /* For rows corresponding to θ ≠ 0, use differentiation matrix */
    for (int k = 1; k <= N; k++) {
        for (int j = 0; j <= N; j++) {
            double coeff = D[k * (N + 1) + j];
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) {
                    if (r == c) {
                        int row_idx = k * n + r;
                        int col_idx = j * n + c;
                        D[row_idx * M + col_idx] = coeff;
                    }
                }
            }
        }
    }
}

/* Simple power iteration to find the rightmost eigenvalue */
__attribute__((unused))
static double power_iteration_rightmost(const double* A_mat, int dim,
                                         int max_iter, double tol) {
    /* Use Arnoldi-like approach: shift and invert to get
     * rightmost eigenvalue. Simplified: use direct eigenvalue
     * search for small systems.
     *
     * For n=1 and N small, we directly compute the characteristic root
     * by solving s + a + a_d e^{-τs} = 0 numerically (Newton method). */
    (void)A_mat; (void)dim; (void)max_iter; (void)tol;
    return 0.0; /* Placeholder for full spectral discretization result */
}

int tds_compute_characteristic_roots(TimeDelaySystem* sys, int max_roots) {
    if (!sys || sys->n_states < 1) return 0;
    if (max_roots > TDS_MAX_ROOTS_DEFAULT) max_roots = TDS_MAX_ROOTS_DEFAULT;

    double tau = (sys->n_delays > 0) ? sys->delays[0]->tau_nominal : 1.0;

    /* For scalar systems: solve s + a + a_d e^{-τs} = 0 iteratively
     * using the Lambert W function branches or numerical continuation. */
    if (sys->n_states == 1) {
        double a = sys->A[0];
        double ad = sys->A_delayed[0];
        sys->n_roots = 0;

        /* Use Newton-Raphson on the complex plane for each frequency range.
         * Transform: s = σ + jω.
         * F(σ,ω) = [σ + a + a_d e^{-στ}cos(ωτ),
         *            ω - a_d e^{-στ}sin(ωτ)] */

        /* Start with the primary branch (rightmost roots).
         * For each interval of ω, solve the system. */
        int found = 0;
        (void)(2.0 * M_PI / tau * 0.5); /* Omega step reference */

        for (int k = 0; k < max_roots && k < 20; k++) {
            /* Initial guess: use asymptotic root locations
             * For large |s|, s ≈ (1/τ) [ln|ad/a| + j(π + 2πk)] */
            double omega_guess = (M_PI + 2.0 * M_PI * k) / tau;
            if (omega_guess < 0) omega_guess = -omega_guess;

            double sigma_guess = 0.0;
            if (fabs(ad) > TDS_EPS) {
                sigma_guess = log(fabs(ad / tau)) / tau;
            }
            if (sigma_guess > 10.0) sigma_guess = 10.0;
            if (sigma_guess < -10.0) sigma_guess = -10.0;

            /* Newton iteration */
            double s_r = sigma_guess, s_i = omega_guess;
            for (int iter = 0; iter < 50; iter++) {
                double e_r = exp(-s_r * tau);
                double e_i_cos = e_r * cos(s_i * tau);
                double e_i_sin = e_r * sin(s_i * tau);

                /* F = [s_r + a + ad*e^{-s_rτ}cos(s_iτ),
                 *      s_i - ad*e^{-s_rτ}sin(s_iτ)] */
                double F_r = s_r + a + ad * e_i_cos;
                double F_i = s_i - ad * e_i_sin;

                if (sqrt(F_r * F_r + F_i * F_i) < TDS_EPS * 1e2) break;

                /* Jacobian:
                 * ∂F_r/∂s_r = 1 - ad*τ*e^{-s_rτ}cos(s_iτ)
                 * ∂F_r/∂s_i = -ad*τ*e^{-s_rτ}sin(s_iτ)
                 * ∂F_i/∂s_r = ad*τ*e^{-s_rτ}sin(s_iτ)
                 * ∂F_i/∂s_i = 1 - ad*τ*e^{-s_rτ}cos(s_iτ) */
                double J11 = 1.0 - ad * tau * e_i_cos;
                double J12 = -ad * tau * e_i_sin;
                double J21 = ad * tau * e_i_sin;
                double J22 = 1.0 - ad * tau * e_i_cos;

                double detJ = J11 * J22 - J12 * J21;
                if (fabs(detJ) < TDS_EPS) break;

                /* J^{-1} * F */
                double ds_r = (J22 * F_r - J12 * F_i) / detJ;
                double ds_i = (-J21 * F_r + J11 * F_i) / detJ;

                s_r -= ds_r;
                s_i -= ds_i;

                if (sqrt(ds_r * ds_r + ds_i * ds_i) < TDS_EPS * 10.0) break;
            }

            /* Store the root */
            if (found < TDS_MAX_ROOTS_DEFAULT) {
                sys->roots_real[found] = s_r;
                sys->roots_imag[found] = s_i;
                found++;
            }

            /* Also try negative omega guess */
            if (k > 0 && found < max_roots) {
                double s_r2 = sigma_guess, s_i2 = -omega_guess;
                for (int iter = 0; iter < 50; iter++) {
                    double e_r2 = exp(-s_r2 * tau);
                    double e_c = e_r2 * cos(s_i2 * tau);
                    double e_s = e_r2 * sin(s_i2 * tau);
                    double F_r2 = s_r2 + a + ad * e_c;
                    double F_i2 = s_i2 - ad * e_s;
                    if (sqrt(F_r2 * F_r2 + F_i2 * F_i2) < TDS_EPS * 1e2) break;

                    double J11_2 = 1.0 - ad * tau * e_c;
                    double J12_2 = -ad * tau * e_s;
                    double J21_2 = ad * tau * e_s;
                    double J22_2 = 1.0 - ad * tau * e_c;
                    double detJ2 = J11_2 * J22_2 - J12_2 * J21_2;
                    if (fabs(detJ2) < TDS_EPS) break;

                    double ds_r2 = (J22_2 * F_r2 - J12_2 * F_i2) / detJ2;
                    double ds_i2 = (-J21_2 * F_r2 + J11_2 * F_i2) / detJ2;
                    s_r2 -= ds_r2; s_i2 -= ds_i2;
                    if (sqrt(ds_r2 * ds_r2 + ds_i2 * ds_i2) < TDS_EPS * 10.0) break;
                }
                if (found < TDS_MAX_ROOTS_DEFAULT) {
                    sys->roots_real[found] = s_r2;
                    sys->roots_imag[found] = s_i2;
                    found++;
                }
            }
        }
        sys->n_roots = found;
        return found;
    }

    /* For higher-dimensional systems, use spectral discretization */
    int N = 20; /* Chebyshev points */
    int M = sys->n_states * (N + 1);
    double* D = (double*)calloc((size_t)(M * M), sizeof(double));
    build_operator_matrix(sys, D, N);

    /* Compute eigenvalues via QR algorithm (simplified).
     * For this implementation, we use direct eigenvalue search
     * for the rightmost roots via Newton on det(sI - A - Ad e^{-τs}). */
    free(D);

    /* Direct search via Newton on the determinant */
    sys->n_roots = 0;
    int found = 0;

    /* Start from multiple initial guesses */
    for (int k = -3; k < 5 && found < max_roots; k++) {
        double omega_guess = (M_PI + 2.0 * M_PI * k) / tau;
        if (omega_guess < 0.5) omega_guess = 0.5;
        double sigma_guess = -1.0;  /* Assume initially stable */

        /* Newton on determinant */
        double s_r = sigma_guess, s_i = omega_guess;
        for (int iter = 0; iter < 30; iter++) {
            double f0 = time_delay_characteristic_eqn(sys, s_r, s_i);
            if (f0 < TDS_EPS * 100) break;

            double eps_d = 1e-7;
            double f_rp = time_delay_characteristic_eqn(sys, s_r + eps_d, s_i);
            double f_ip = time_delay_characteristic_eqn(sys, s_r, s_i + eps_d);
            double dF_dr = (f_rp - f0) / eps_d;
            double dF_di = (f_ip - f0) / eps_d;

            double grad_norm2 = dF_dr * dF_dr + dF_di * dF_di;
            if (grad_norm2 < TDS_EPS) break;

            double step = f0 / grad_norm2;
            s_r -= step * dF_dr;
            s_i -= step * dF_di;

            if (fabs(s_r) > 20.0) break;  /* Diverge guard */
        }

        /* Check if this is a new root (not duplicate) */
        bool is_new = true;
        for (int i = 0; i < found; i++) {
            double dr = s_r - sys->roots_real[i];
            double di = s_i - sys->roots_imag[i];
            if (sqrt(dr * dr + di * di) < 0.01) { is_new = false; break; }
        }
        if (is_new && found < TDS_MAX_ROOTS_DEFAULT) {
            sys->roots_real[found] = s_r;
            sys->roots_imag[found] = s_i;
            found++;
        }
    }
    sys->n_roots = found;
    return found;
}

bool tds_is_characteristic_root(const TimeDelaySystem* sys,
                                 double sigma, double omega, double eps) {
    double val = time_delay_characteristic_eqn(sys, sigma, omega);
    return val < eps;
}

double tds_spectral_abscissa(TimeDelaySystem* sys) {
    if (!sys || sys->n_roots == 0) {
        if (sys) tds_compute_characteristic_roots(sys, 8);
    }
    if (!sys || sys->n_roots == 0) return 0.0;

    double max_sigma = -INFINITY;
    for (int i = 0; i < sys->n_roots; i++) {
        if (sys->roots_real[i] > max_sigma)
            max_sigma = sys->roots_real[i];
    }
    sys->spectral_abscissa = max_sigma;

    /* Classify stability */
    if (max_sigma < -TDS_EPS) sys->stability_class = DELAY_STABLE;
    else if (fabs(max_sigma) < TDS_EPS * 1e3) sys->stability_class = DELAY_MARGINALLY_STABLE;
    else sys->stability_class = DELAY_UNSTABLE;

    return max_sigma;
}
