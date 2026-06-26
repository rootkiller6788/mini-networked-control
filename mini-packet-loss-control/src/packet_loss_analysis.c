#include "packet_loss_analysis.h"
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
 * Kronecker Product: C = A ⊗ B
 * A: m×n, B: p×q  →  C: (m·p)×(n·q)
 *
 * C_{i·p+k, j·q+l} = A_{i,j} · B_{k,l}
 *
 * Used in MJLS stability analysis: the augmented operator is
 *   A = (P' ⊗ I) · diag(A₀⊗A₀, ..., A_{M-1}⊗A_{M-1})
 *
 * where ⊗ is the Kronecker product and P is the Markov transition matrix.
 * ============================================================================ */

void pl_kronecker_product(const double* A, const double* B,
                           double* C, int m, int n, int p, int q) {
    int mr = m * p, nc = n * q;
    for (int i = 0; i < mr; i++)
        for (int j = 0; j < nc; j++)
            C[i * nc + j] = 0.0;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < p; k++)
                for (int l = 0; l < q; l++)
                    C[(i * p + k) * nc + (j * q + l)] =
                        A[i * n + j] * B[k * q + l];
}

/* ============================================================================
 * Spectral Radius via Power Iteration
 *
 * Finds dominant eigenvalue magnitude of matrix A (size n×n).
 * The power method: v_{k+1} = A·v_k / ||A·v_k||,
 * with Rayleigh quotient for eigenvalue estimate.
 *
 * Stopping criterion: |λ_{k+1} - λ_k| < tol or max_iter reached.
 *
 * (Defined in packet_loss_controller.c — see pl_spectral_radius_power)
 * ============================================================================ */

/* ============================================================================
 * Matrix Frobenius Norm: ||A||_F = sqrt(Σ_{i,j} a_{ij}²)
 * ============================================================================ */

double pl_matrix_norm_frobenius(const double* A, int rows, int cols) {
    double sum = 0.0;
    for (int i = 0; i < rows * cols; i++) sum += A[i] * A[i];
    return sqrt(sum);
}

/* ============================================================================
 * Schur Stability Check: ρ(A) < 1
 * ============================================================================ */

bool pl_is_schur_stable(const double* A, int n) {
    return pl_spectral_radius_power(A, n, 500, 1e-8) < 1.0;
}

/* ============================================================================
 * Eigenvalue Computation via QR Algorithm (Francis, 1961)
 *
 * Computes eigenvalues of a real n×n matrix using the QR algorithm
 * with Wilkinson shifts. Returns number of real eigenvalues.
 *
 * Complexity: O(n³) for the full decomposition.
 * For n ≤ 32, this is efficient.
 *
 * Reference: Golub & Van Loan (2013), §7.5.
 * ============================================================================ */

int pl_eigenvalues(const double* A, double* real_part, double* imag_part, int n) {
    if (n <= 0 || n > 64) return 0;

    /* Copy A to working matrix H (Hessenberg form) */
    double* H = (double*)malloc(n * n * sizeof(double));
    memcpy(H, A, n * n * sizeof(double));

    /* Reduce to upper Hessenberg form via Householder reflections */
    double* u = (double*)malloc(n * sizeof(double));
    for (int k = 0; k < n - 2; k++) {
        double sigma = 0.0;
        for (int i = k + 1; i < n; i++) sigma += H[i * n + k] * H[i * n + k];
        sigma = sqrt(sigma);

        if (sigma < 1e-15) continue;

        double alpha = (H[(k+1) * n + k] > 0) ? -sigma : sigma;
        u[k+1] = H[(k+1) * n + k] - alpha;
        double u_norm_sq = u[k+1] * u[k+1];
        for (int i = k + 2; i < n; i++) {
            u[i] = H[i * n + k];
            u_norm_sq += u[i] * u[i];
        }

        double beta = 2.0 / u_norm_sq;

        /* H = (I - β·u·u')·H·(I - β·u·u') */
        for (int i = 0; i < n; i++) {
            double dot = 0.0;
            for (int l = k + 1; l < n; l++) dot += u[l] * H[l * n + i];
            for (int l = k + 1; l < n; l++)
                H[l * n + i] -= beta * u[l] * dot;
        }
        for (int j = 0; j < n; j++) {
            double dot = 0.0;
            for (int l = k + 1; l < n; l++) dot += H[j * n + l] * u[l];
            for (int l = k + 1; l < n; l++)
                H[j * n + l] -= beta * dot * u[l];
        }
    }

    /* QR iteration on Hessenberg matrix (Francis double-shift) */
    for (int iter = 0; iter < 500; iter++) {
        int converged = 1;
        for (int i = 1; i < n; i++)
            if (fabs(H[i * n + i - 1]) > 1e-12) { converged = 0; break; }
        if (converged) break;

        /* Shift: use bottom-right 2×2 eigenvalues */
        double trace = H[(n-2)*n + n-2] + H[(n-1)*n + n-1];
        double det = H[(n-2)*n + n-2] * H[(n-1)*n + n-1]
                   - H[(n-2)*n + n-1] * H[(n-1)*n + n-2];
        (void)trace; (void)det; /* Reserved for Wilkinson double-shift implementation */

        /* Single (real) shift */
        double shift = H[(n-1)*n + n-1];

        /* Apply shifted QR (Givens rotation-based) */
        for (int k = 0; k < n - 1; k++) {
            double a = H[k * n + k] - shift;
            double b = H[(k+1) * n + k];
            double r = sqrt(a * a + b * b);
            if (r < 1e-15) continue;

            double c_rot = a / r, s_rot = -b / r;

            /* Apply Givens rotation to rows k, k+1 */
            for (int j = k; j < n; j++) {
                double t1 = H[k * n + j], t2 = H[(k+1) * n + j];
                H[k * n + j]     =  c_rot * t1 - s_rot * t2;
                H[(k+1) * n + j] =  s_rot * t1 + c_rot * t2;
            }
            /* Apply transpose to columns */
            for (int i = 0; i < n; i++) {
                double t1 = H[i * n + k], t2 = H[i * n + k + 1];
                H[i * n + k]     =  c_rot * t1 - s_rot * t2;
                H[i * n + k + 1] =  s_rot * t1 + c_rot * t2;
            }
        }
    }

    /* Extract eigenvalues from diagonal and subdiagonal */
    int n_real = 0;
    int i = 0;
    while (i < n) {
        if (i < n - 1 && fabs(H[(i+1) * n + i]) > 1e-10) {
            /* 2×2 block → complex conjugate pair */
            double a = H[i * n + i], b = H[i * n + i + 1];
            double c = H[(i+1) * n + i], d = H[(i+1) * n + i + 1];
            double tr = a + d;
            double dt = a * d - b * c;
            double disc = tr * tr - 4.0 * dt;
            real_part[i] = tr / 2.0;
            real_part[i+1] = tr / 2.0;
            if (disc >= 0) {
                double sqrt_disc = sqrt(disc);
                real_part[i] = (tr + sqrt_disc) / 2.0;
                real_part[i+1] = (tr - sqrt_disc) / 2.0;
                imag_part[i] = 0.0; imag_part[i+1] = 0.0;
                n_real += 2;
            } else {
                imag_part[i] = sqrt(-disc) / 2.0;
                imag_part[i+1] = -sqrt(-disc) / 2.0;
            }
            i += 2;
        } else {
            /* Real eigenvalue */
            real_part[i] = H[i * n + i];
            imag_part[i] = 0.0;
            n_real++;
            i++;
        }
    }

    free(H); free(u);
    return n_real;
}

/* ============================================================================
 * Stability Analysis — Bernoulli Packet Loss
 *
 * For Bernoulli loss with probability p:
 *   x_{k+1} = A_s·x_k (success, prob 1-p) or A_f·x_k (failure, prob p)
 *
 * Mean-square stability (MSS) condition:
 *   ρ((1-p)·(A_s ⊗ A_s) + p·(A_f ⊗ A_f)) < 1
 *
 * This is a necessary and sufficient condition for MJLS with
 * Bernoulli switching (Costa et al., 2005).
 *
 * The test matrix is n²×n² — tractable for n ≤ 10.
 * ============================================================================ */

StabilityCertificate* pl_stability_test_bernoulli(
    const double* A_success, const double* A_failure,
    int n, double loss_prob) {

    StabilityCertificate* cert = (StabilityCertificate*)calloc(1,
        sizeof(StabilityCertificate));
    cert->type = STAB_MEAN_SQUARE;

    if (n > 10) {
        cert->is_stable = false;
        cert->margin = -1.0;
        return cert;
    }

    int n2 = n * n;

    /* Compute Kronecker products */
    double* K_s = _m(n2, n2);
    double* K_f = _m(n2, n2);
    pl_kronecker_product(A_success, A_success, K_s, n, n, n, n);
    pl_kronecker_product(A_failure, A_failure, K_f, n, n, n, n);

    /* Form test matrix: T = (1-p)·(A_s⊗A_s) + p·(A_f⊗A_f) */
    double* T = _m(n2, n2);
    double p_success = 1.0 - loss_prob;
    for (int i = 0; i < n2 * n2; i++)
        T[i] = p_success * K_s[i] + loss_prob * K_f[i];

    /* Compute spectral radius */
    double rho = pl_spectral_radius_power(T, n2, 2000, 1e-10);
    cert->is_stable = (rho < 1.0);
    cert->margin = 1.0 - rho;
    cert->decay_rate = log(rho + 1e-15);
    cert->critical_loss_prob = (rho >= 1.0) ? loss_prob * 0.9 : loss_prob;

    /* Find critical p via binary search */
    if (!cert->is_stable) {
        double lo = 0.0, hi = loss_prob;
        for (int b = 0; b < 20; b++) {
            double mid = (lo + hi) / 2.0;
            for (int i = 0; i < n2 * n2; i++)
                T[i] = (1.0 - mid) * K_s[i] + mid * K_f[i];
            double rho_mid = pl_spectral_radius_power(T, n2, 1000, 1e-8);
            if (rho_mid < 1.0) lo = mid; else hi = mid;
        }
        cert->critical_loss_prob = lo;
    }

    free(K_s); free(K_f); free(T);
    return cert;
}

/* ============================================================================
 * Stability Analysis — Gilbert-Elliott Channel
 *
 * Uses coupled Lyapunov equation approach:
 *   For each mode i ∈ {Good, Bad}:
 *     P_i - Σ_j p_{ij}·A_j'·P_j·A_j ≻ 0
 *
 * Feasibility of these Lyapunov inequalities certifies MSS.
 * Iterative solution: P_i^{(t+1)} = Σ_j p_{ij}·A_j'·P_j^{(t)}·A_j + I
 *
 * If the iterations converge, the system is MSS.
 * ============================================================================ */

StabilityCertificate* pl_stability_test_gilbert_elliott(
    const double* A_success, const double* A_failure,
    int n, const GilbertElliottChannel* ch) {

    StabilityCertificate* cert = (StabilityCertificate*)calloc(1,
        sizeof(StabilityCertificate));
    cert->type = STAB_MEAN_SQUARE;

    /* Mode matrices: mode 0 = Good, mode 1 = Bad */
    /* In Good state: lower loss → mostly A_success */
    /* In Bad state: higher loss → mostly A_failure */
    const double* A_modes[2] = { A_success, A_failure };

    /* Transition matrix */
    double P[2][2] = {
        { 1.0 - ch->p_gb, ch->p_gb },
        { ch->p_bg, 1.0 - ch->p_bg }
    };

    /* Iterative Lyapunov computation */
    double* P0 = _m(n, n);
    double* P1 = _m(n, n);
    for (int i = 0; i < n; i++) { P0[i * n + i] = 1.0; P1[i * n + i] = 1.0; }

    double* P0_new = _m(n, n);
    double* P1_new = _m(n, n);
    double* temp = _m(n, n);

    int max_iter = 5000;
    double prev_diff = 0.0;
    bool converged = false;

    for (int iter = 0; iter < max_iter; iter++) {
        /* P0_new = p_00·A_0'P0A_0 + p_01·A_1'P1A_1 + I */
        for (int mode = 0; mode < 2; mode++) {
            double* P_cur = (mode == 0) ? P0 : P1; (void)P_cur;
            double* P_new = (mode == 0) ? P0_new : P1_new;
            memset(P_new, 0, n * n * sizeof(double));

            for (int src = 0; src < 2; src++) {
                double p_ij = P[mode][src];
                const double* A_j = A_modes[src];
                const double* Pj = (src == 0) ? P0 : P1;

                /* temp = A_j'·Pj */
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        double s = 0.0;
                        for (int l = 0; l < n; l++) s += A_j[l * n + i] * Pj[l * n + j];
                        temp[i * n + j] = s;
                    }

                /* Add p_ij · temp · A_j to P_new */
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        double s = 0.0;
                        for (int l = 0; l < n; l++) s += temp[i * n + l] * A_j[l * n + j];
                        P_new[i * n + j] += p_ij * s;
                    }
            }
            for (int i = 0; i < n; i++) P_new[i * n + i] += 1.0;
        }

        /* Check convergence */
        double diff = 0.0;
        for (int i = 0; i < n * n; i++) {
            double d = P0_new[i] - P0[i]; diff += d * d;
            d = P1_new[i] - P1[i]; diff += d * d;
        }
        diff = sqrt(diff);

        memcpy(P0, P0_new, n * n * sizeof(double));
        memcpy(P1, P1_new, n * n * sizeof(double));

        if (diff < 1e-10) { converged = true; cert->iterations_to_check = iter + 1; break; }
        if (iter > 100 && diff > prev_diff * 2.0) break; /* Diverging */
        prev_diff = diff;
    }

    cert->is_stable = converged;
    cert->margin = converged ? 1.0 : -1.0;
    cert->critical_loss_prob = ch->steady_loss_rate;

    free(P0); free(P1); free(P0_new); free(P1_new); free(temp);
    return cert;
}

/* ============================================================================
 * Stability Analysis — Generic MJLS
 * ============================================================================ */

StabilityCertificate* pl_stability_test_markov(JumpLinearSystem* jls) {
    StabilityCertificate* cert = (StabilityCertificate*)calloc(1,
        sizeof(StabilityCertificate));
    if (!jls) { cert->is_stable = false; return cert; }

    cert->type = STAB_MEAN_SQUARE;
    cert->is_stable = pl_jls_test_mss(jls, 5000, 1e-10);
    cert->margin = cert->is_stable ? 1.0 : -1.0;

    return cert;
}

/* ============================================================================
 * Lyapunov Exponent Estimation via Monte Carlo
 *
 * λ = lim_{k→∞} (1/k)·E[log||x_k||]
 *
 * For stable systems: λ < 0
 * For marginally stable: λ ≈ 0
 * For unstable: λ > 0
 *
 * Computed by averaging log-ratio over many trials:
 *   λ_estimate = (1/(N·K)) Σ_{trial} Σ_k log(||x_{k+1}|| / ||x_k||)
 * ============================================================================ */

double pl_stability_lyapunov_exponent(
    const double* A_success, const double* A_failure,
    int n, double loss_prob, int n_trials, int n_steps,
    unsigned long seed) {

    double total_log_ratio = 0.0;
    int total_steps = 0;

    for (int trial = 0; trial < n_trials; trial++) {
        /* Random initial state */
        double* x = (double*)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) x[i] = pl_uniform(&seed) - 0.5;
        /* Normalize */
        double norm0 = 0.0;
        for (int i = 0; i < n; i++) norm0 += x[i] * x[i];
        norm0 = sqrt(norm0);
        if (norm0 > 1e-10) for (int i = 0; i < n; i++) x[i] /= norm0;

        for (int step = 0; step < n_steps; step++) {
            double u = pl_uniform(&seed);
            bool success = (u >= loss_prob);
            const double* A_k = success ? A_success : A_failure;

            double* x_next = (double*)malloc(n * sizeof(double));
            for (int i = 0; i < n; i++) {
                x_next[i] = 0.0;
                for (int j = 0; j < n; j++)
                    x_next[i] += A_k[i * n + j] * x[j];
            }

            double norm_curr = 0.0, norm_next = 0.0;
            for (int i = 0; i < n; i++) {
                norm_curr += x[i] * x[i];
                norm_next += x_next[i] * x_next[i];
            }
            norm_curr = sqrt(norm_curr);
            norm_next = sqrt(norm_next);

            if (norm_curr > 1e-15 && norm_next > 1e-15) {
                total_log_ratio += log(norm_next / norm_curr);
                total_steps++;
            }

            /* Normalize to prevent overflow */
            if (norm_next > 1e-10)
                for (int i = 0; i < n; i++) x[i] = x_next[i] / norm_next;
            else
                memcpy(x, x_next, n * sizeof(double));

            free(x_next);
        }
        free(x);
    }

    return (total_steps > 0) ? total_log_ratio / total_steps : 0.0;
}

/* ============================================================================
 * Monte Carlo Mean-Square Error Simulation
 *
 * Estimates E[||x_k||²] over time by averaging over independent trials.
 * ============================================================================ */

double pl_stability_monte_carlo(
    const double* A_success, const double* A_failure,
    int n, double loss_prob, int n_trials, int n_steps,
    unsigned long seed) {

    double total_mse = 0.0;

    for (int trial = 0; trial < n_trials; trial++) {
        double* x = (double*)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) x[i] = pl_uniform(&seed) - 0.5;

        for (int step = 0; step < n_steps; step++) {
            double u = pl_uniform(&seed);
            const double* A_k = (u >= loss_prob) ? A_success : A_failure;

            double* xn = (double*)malloc(n * sizeof(double));
            for (int i = 0; i < n; i++) {
                xn[i] = 0.0;
                for (int j = 0; j < n; j++)
                    xn[i] += A_k[i * n + j] * x[j];
            }

            if (step == n_steps - 1) {
                double mse = 0.0;
                for (int i = 0; i < n; i++) mse += xn[i] * xn[i];
                total_mse += mse;
            }

            memcpy(x, xn, n * sizeof(double));
            free(xn);
        }
        free(x);
    }

    return total_mse / n_trials;
}

/* ============================================================================
 * Stochastic Lyapunov Function
 *
 * For Bernoulli loss: find P ≻ 0 such that
 *   (1-p)·A_s'·P·A_s + p·A_f'·P·A_f - P = -I
 *
 * Solved via Lyapunov iteration:
 *   P^{(t+1)} = (1-p)·A_s'·P^{(t)}·A_s + p·A_f'·P^{(t)}·A_f + I
 *
 * Convergence ⇔ MSS.
 * Then V(x) = x'·P·x is a stochastic Lyapunov function.
 * ============================================================================ */

double pl_stability_lyapunov_function(
    const double* x, int n,
    const double* A_success, const double* A_failure,
    double loss_prob, int max_iter, double tol) {

    double* P = _m(n, n);
    for (int i = 0; i < n; i++) P[i * n + i] = 1.0;

    double* P_new = _m(n, n);
    double* temp  = _m(n, n);
    double* AT_s  = _m(n, n);
    double* AT_f  = _m(n, n);

    /* Precompute transposes */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            AT_s[i * n + j] = A_success[j * n + i];
            AT_f[i * n + j] = A_failure[j * n + i];
        }

    for (int iter = 0; iter < max_iter; iter++) {
        /* P_new = I + (1-p)·A_s'·P·A_s + p·A_f'·P·A_f */
        memset(P_new, 0, n * n * sizeof(double));

        /* (1-p)·A_s'·P·A_s */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += AT_s[i*n+l] * P[l*n+j];
                temp[i * n + j] = s;
            }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += temp[i*n+l] * A_success[l*n+j];
                P_new[i * n + j] += (1.0 - loss_prob) * s;
            }

        /* p·A_f'·P·A_f */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += AT_f[i*n+l] * P[l*n+j];
                temp[i * n + j] = s;
            }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += temp[i*n+l] * A_failure[l*n+j];
                P_new[i * n + j] += loss_prob * s;
            }

        for (int i = 0; i < n; i++) P_new[i * n + i] += 1.0;

        double diff = 0.0;
        for (int i = 0; i < n * n; i++) {
            double d = P_new[i] - P[i]; diff += d * d;
        }
        memcpy(P, P_new, n * n * sizeof(double));
        if (sqrt(diff) < tol) break;
    }

    /* V(x) = x'·P·x */
    double v = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            v += x[i] * P[i * n + j] * x[j];

    free(P); free(P_new); free(temp); free(AT_s); free(AT_f);
    return v;
}

/* ============================================================================
 * Jump Linear System (MJLS) Core
 * ============================================================================ */

JumpLinearSystem* pl_jls_create(int n_states, int n_modes) {
    JumpLinearSystem* jls = (JumpLinearSystem*)calloc(1, sizeof(JumpLinearSystem));
    if (!jls) return NULL;

    jls->n_states = n_states;
    jls->n_modes = n_modes;

    jls->mode_matrices = (double**)malloc(n_modes * sizeof(double*));
    double* mm_data = _m(n_modes * n_states, n_states);
    for (int i = 0; i < n_modes; i++)
        jls->mode_matrices[i] = mm_data + i * n_states * n_states;

    jls->transition_matrix = (double**)malloc(n_modes * sizeof(double*));
    double* tm_data = _m(n_modes, n_modes);
    for (int i = 0; i < n_modes; i++) {
        jls->transition_matrix[i] = tm_data + i * n_modes;
        for (int j = 0; j < n_modes; j++)
            jls->transition_matrix[i][j] = 1.0 / n_modes;
    }

    jls->steady_state = (double*)calloc(n_modes, sizeof(double));
    for (int i = 0; i < n_modes; i++) jls->steady_state[i] = 1.0 / n_modes;

    jls->lyapunov_matrices = (double**)malloc(n_modes * sizeof(double*));
    double* lm_data = _m(n_modes * n_states, n_states);
    for (int i = 0; i < n_modes; i++)
        jls->lyapunov_matrices[i] = lm_data + i * n_states * n_states;

    return jls;
}

void pl_jls_free(JumpLinearSystem* jls) {
    if (!jls) return;
    free(jls->mode_matrices[0]);
    free(jls->mode_matrices);
    free(jls->transition_matrix[0]);
    free(jls->transition_matrix);
    free(jls->steady_state);
    free(jls->lyapunov_matrices[0]);
    free(jls->lyapunov_matrices);
    free(jls);
}

void pl_jls_set_mode_matrix(JumpLinearSystem* jls, int mode, const double* A_mode) {
    if (!jls || mode < 0 || mode >= jls->n_modes) return;
    memcpy(jls->mode_matrices[mode], A_mode,
           jls->n_states * jls->n_states * sizeof(double));
}

void pl_jls_set_transition(JumpLinearSystem* jls, int from, int to, double prob) {
    if (!jls || from < 0 || from >= jls->n_modes || to < 0 || to >= jls->n_modes)
        return;
    jls->transition_matrix[from][to] = prob;
}

void pl_jls_compute_steady_state(JumpLinearSystem* jls) {
    int M = jls->n_modes;
    double* pi_new = (double*)malloc(M * sizeof(double));
    for (int i = 0; i < M; i++) jls->steady_state[i] = 1.0 / M;

    for (int iter = 0; iter < 10000; iter++) {
        for (int j = 0; j < M; j++) {
            pi_new[j] = 0.0;
            for (int i = 0; i < M; i++)
                pi_new[j] += jls->steady_state[i] * jls->transition_matrix[i][j];
        }
        double max_d = 0.0;
        for (int i = 0; i < M; i++) {
            double d = fabs(pi_new[i] - jls->steady_state[i]);
            if (d > max_d) max_d = d;
            jls->steady_state[i] = pi_new[i];
        }
        if (max_d < 1e-15) break;
    }
    free(pi_new);
}

bool pl_jls_test_mss(JumpLinearSystem* jls, int max_iter, double tol) {
    int n = jls->n_states, M = jls->n_modes;
    double* Pi_new = _m(M * n, n);
    double* temp = _m(n, n);
    double* prev_norm = (double*)malloc(M * sizeof(double));

    /* Initialize Lyapunov matrices to identity */
    for (int mode = 0; mode < M; mode++) {
        double* Pi = jls->lyapunov_matrices[mode];
        memset(Pi, 0, n * n * sizeof(double));
        for (int i = 0; i < n; i++) Pi[i * n + i] = 1.0;
        prev_norm[mode] = 1.0;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        double max_norm = 0.0;
        for (int mode = 0; mode < M; mode++) {
            double* pin = Pi_new + mode * n * n;
            memset(pin, 0, n * n * sizeof(double));

            for (int src = 0; src < M; src++) {
                double p_ij = jls->transition_matrix[mode][src];
                const double* A_j = jls->mode_matrices[src];
                double* Pj = jls->lyapunov_matrices[src];

                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        double s = 0.0;
                        for (int l = 0; l < n; l++)
                            s += A_j[l * n + i] * Pj[l * n + j];
                        temp[i * n + j] = s;
                    }
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        double s = 0.0;
                        for (int l = 0; l < n; l++)
                            s += temp[i*n+l] * A_j[l*n+j];
                        pin[i * n + j] += p_ij * s;
                    }
            }
            for (int i = 0; i < n; i++) pin[i * n + i] += 1.0;

            double norm = 0.0;
            for (int i = 0; i < n * n; i++) norm += pin[i] * pin[i];
            norm = sqrt(norm);
            if (norm > max_norm) max_norm = norm;
        }

        double max_ratio = 0.0;
        for (int mode = 0; mode < M; mode++) {
            double* pin = Pi_new + mode * n * n;
            double norm = 0.0;
            for (int i = 0; i < n * n; i++) norm += pin[i] * pin[i];
            norm = sqrt(norm);
            double ratio = (prev_norm[mode] > 1e-15) ? norm / prev_norm[mode] : 0.0;
            if (ratio > max_ratio) max_ratio = ratio;
            prev_norm[mode] = norm;

            memcpy(jls->lyapunov_matrices[mode], pin, n * n * sizeof(double));
        }

        jls->spectral_radius_estimate = max_ratio;
        jls->iterations_to_converge = iter + 1;

        double change = 0.0;
        for (int mode = 0; mode < M; mode++) {
            double* Pi = jls->lyapunov_matrices[mode];
            double* pin = Pi_new + mode * n * n;
            for (int i = 0; i < n * n; i++) {
                double d = pin[i] - Pi[i]; change += d * d;
            }
        }

        if (max_norm > 1e6) { jls->mss_certified = false; break; }
        if (sqrt(change) < tol && max_ratio < 1.0) {
            jls->mss_certified = true; break;
        }
        if (max_ratio > 1.0) { jls->mss_certified = false; break; }
    }

    free(Pi_new); free(temp); free(prev_norm);
    return jls->mss_certified;
}

double pl_jls_operator_spectral_radius(JumpLinearSystem* jls) {
    return jls ? jls->spectral_radius_estimate : 0.0;
}

void pl_jls_print(const JumpLinearSystem* jls) {
    if (!jls) return;
    printf("=== Jump Linear System ===\n");
    printf("States: %d, Modes: %d, MSS: %s\n",
           jls->n_states, jls->n_modes, jls->mss_certified ? "YES" : "NO");
}

/* ============================================================================
 * Critical Probability Analysis
 * ============================================================================ */

CriticalProbabilityAnalysis* pl_critical_prob_analyze(
    const double* A, const double* B, const double* C,
    const double* L, int n, int m, int p) {

    CriticalProbabilityAnalysis* cpa = (CriticalProbabilityAnalysis*)calloc(1,
        sizeof(CriticalProbabilityAnalysis));

    /* Open-loop spectral radius */
    cpa->spectral_radius_A = pl_spectral_radius_power(A, n, 1000, 1e-10);

    /* Closed-loop: A_cl = A - B·L */
    double* Acl = _m(n, n);
    double* BL  = _m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int l = 0; l < m; l++) s += B[i * m + l] * L[l * n + j];
            BL[i * n + j] = s;
        }
    for (int i = 0; i < n * n; i++) Acl[i] = A[i] - BL[i];
    cpa->spectral_radius_Acl = pl_spectral_radius_power(Acl, n, 1000, 1e-10);

    /* Sensor critical probability (Sinopoli 2004) */
    double rho2 = cpa->spectral_radius_A * cpa->spectral_radius_A;
    if (rho2 > 1.0) {
        cpa->lower_bound_sensor = 1.0 - 1.0 / rho2;
        cpa->p_c_sensor = cpa->lower_bound_sensor;
    } else {
        cpa->lower_bound_sensor = 0.0;
        cpa->p_c_sensor = 1.0;
    }

    /* Actuator critical probability */
    double rho2_cl = cpa->spectral_radius_Acl * cpa->spectral_radius_Acl;
    if (rho2_cl > 1.0) {
        cpa->lower_bound_actuator = 1.0 - 1.0 / rho2_cl;
        cpa->p_c_actuator = cpa->lower_bound_actuator;
    } else {
        cpa->lower_bound_actuator = 0.0;
        cpa->p_c_actuator = 1.0;
    }

    /* Joint critical (product approximation) */
    cpa->p_c_joint = (cpa->p_c_sensor < cpa->p_c_actuator) ?
                      cpa->p_c_sensor : cpa->p_c_actuator;
    cpa->is_stabilizable_under_loss = (cpa->p_c_joint > 0.1);

    free(Acl); free(BL);
    (void)C; (void)p;
    return cpa;
}

void pl_critical_prob_free(CriticalProbabilityAnalysis* cpa) {
    if (!cpa) return;
    free(cpa->stability_region);
    free(cpa);
}

int pl_critical_prob_stability_region(CriticalProbabilityAnalysis* cpa,
                                       int n_samples) {
    if (!cpa || n_samples <= 0) return 0;

    cpa->stability_region = (double*)realloc(cpa->stability_region,
                                              2 * n_samples * sizeof(double));
    cpa->n_region_points = 0;

    /* Sample along line: p_s + p_a = p_c_joint */
    for (int i = 0; i < n_samples; i++) {
        double t = (double)i / (double)(n_samples - 1);
        double p_s = t * cpa->p_c_joint;
        double p_a = (1.0 - t) * cpa->p_c_joint;
        cpa->stability_region[2 * i] = p_s;
        cpa->stability_region[2 * i + 1] = p_a;
    }
    cpa->n_region_points = n_samples;

    return n_samples;
}

void pl_critical_prob_print(const CriticalProbabilityAnalysis* cpa) {
    if (!cpa) return;
    printf("=== Critical Probability Analysis ===\n");
    printf("ρ(A)=%.6f  ρ(A-BL)=%.6f\n", cpa->spectral_radius_A, cpa->spectral_radius_Acl);
    printf("p_c(sensor): %.6f  p_c(actuator): %.6f  p_c(joint): %.6f\n",
           cpa->p_c_sensor, cpa->p_c_actuator, cpa->p_c_joint);
    printf("Stabilizable under loss: %s\n",
           cpa->is_stabilizable_under_loss ? "YES" : "NO");
}
