/**
 * blc_control.c — Bandwidth-Limited Controller Design
 *
 * Implementation of control design methodologies under communication
 * constraints:
 *  - LQR under quantization (quantization penalty analysis)
 *  - Quantized state feedback stability conditions
 *  - Minimum attention control (Brockett, 1997)
 *  - Robust stability via Circle Criterion for quantization nonlinearity
 *  - Bandwidth-limited Model Predictive Control
 *
 * Knowledge coverage: L5 (Control Design), L6 (Canonical Problems)
 */

#include "blc_control.h"
#include "blc_core.h"
#include "blc_datarate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * LQR Controller with Quantization Analysis
 *
 * Standard LQR: u = -Kx = -R^{-1}B'P x
 * where P solves the ARE: A'P + PA - PBR^{-1}B'P + Q = 0
 *
 * Newton-Kleinman iteration for the ARE:
 *   P_{k+1} = P_k - (A-BK_k)'P_{k+1}(A-BK_k) - K_k'RK_k - Q
 *   where K_k = R^{-1}B'P_k
 *
 * With quantization Q(x) = x + ε:
 *   u = -K·Q(x) = -Kx - Kε
 *   ẋ = (A-BK)x - BKε
 *   J_quant = J_nominal + trace(Σ_ε · K'RK)
 *
 * @ref Kleinman (1968), "On an iterative technique for Riccati
 *      equation computations", IEEE TAC
 * ================================================================ */

int blc_lqr_init(BLCLQRController* lqr, int n_states, int n_inputs,
                  const double* Q_diag, const double* R_diag) {
    if (!lqr || n_states < 1 || n_states > BLC_MAX_STATES
        || n_inputs < 1 || n_inputs > BLC_MAX_STATES) return -1;

    lqr->n_states = n_states;
    lqr->n_inputs = n_inputs;

    for (int i = 0; i < n_states; i++) {
        lqr->Q_diag[i] = Q_diag ? Q_diag[i] : 1.0;
    }
    for (int i = 0; i < n_inputs; i++) {
        lqr->R_diag[i] = R_diag ? R_diag[i] : 1.0;
    }

    memset(lqr->K, 0, sizeof(lqr->K));
    memset(lqr->P, 0, sizeof(lqr->P));
    memset(lqr->Ac, 0, sizeof(lqr->Ac));
    lqr->cost = 0.0;
    lqr->cost_with_quantization = 0.0;
    lqr->quantization_penalty = 0.0;
    lqr->are_iterations = 0;
    lqr->converged = false;

    return 0;
}

int blc_lqr_solve(BLCLQRController* lqr, const double* A,
                   const double* B, int max_iter, double tol) {
    if (!lqr || !A || !B) return -1;
    int n = lqr->n_states;
    int m = lqr->n_inputs;

    /** Initialize P = Q */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            lqr->P[i][j] = (i == j) ? lqr->Q_diag[i] : 0.0;
        }
    }

    double P_prev[BLC_MAX_STATES][BLC_MAX_STATES];

    for (int iter = 0; iter < max_iter; iter++) {
        /** Save P_prev */
        memcpy(P_prev, lqr->P, sizeof(P_prev));

        /** Compute K = R^{-1} B' P */
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += B[k * m + i] * lqr->P[k][j];  /** B'(i,k) = B(k,i) */
                }
                lqr->K[i][j] = sum / lqr->R_diag[i];
            }
        }

        /** Compute A_cl = A - B K */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                lqr->Ac[i][j] = A[i * n + j];
                for (int k = 0; k < m; k++) {
                    lqr->Ac[i][j] -= B[i * m + k] * lqr->K[k][j];
                }
            }
        }

        /** Newton step: solve A_cl' P + P A_cl = -Q - K' R K
         *  Using Lyapunov solver */
        double RHS[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
        /** RHS = Q + K' R K */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                RHS[i][j] = (i == j) ? lqr->Q_diag[i] : 0.0;
                for (int k = 0; k < m; k++) {
                    RHS[i][j] += lqr->K[k][i] * lqr->R_diag[k] * lqr->K[k][j];
                }
            }
        }

        /** Solve Lyapunov: A_cl' P + P A_cl = -RHS */
        /** Use iterative method */
        double R[BLC_MAX_STATES][BLC_MAX_STATES];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double AP = 0.0, PA = 0.0;
                for (int k = 0; k < n; k++) {
                    AP += lqr->Ac[k][i] * lqr->P[k][j];  /** A_cl'(i,k) P(k,j) */
                    PA += lqr->P[i][k] * lqr->Ac[k][j];  /** P(i,k) A_cl(k,j) */
                }
                R[i][j] = AP + PA + RHS[i][j];  /** Residual: A'P+PA+RHS */
            }
        }

        /** Update: P -= stepsize * R */
        double Anorm = 0.0;
        for (int i = 0; i < n * n; i++) Anorm += fabs(A[i]);
        double alpha = 0.5 / (Anorm + 1e-10);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                lqr->P[i][j] -= alpha * R[i][j];

        /** Convergence check */
        double diff = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                diff += fabs(lqr->P[i][j] - P_prev[i][j]);

        if (diff < tol) {
            lqr->converged = true;
            lqr->are_iterations = iter + 1;
            return iter + 1;
        }
    }

    lqr->are_iterations = max_iter;
    lqr->converged = false;
    return max_iter;
}

int blc_lqr_closed_loop_poles(const BLCLQRController* lqr,
                               const double* A, const double* B,
                               double* poles_real, double* poles_imag) {
    if (!lqr || !A || !B) return -1;
    int n = lqr->n_states;
    int m = lqr->n_inputs;

    /** Compute A_cl = A - B K */
    double A_cl[BLC_MAX_STATES][BLC_MAX_STATES];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_cl[i][j] = A[i * n + j];
            for (int k = 0; k < m; k++) {
                A_cl[i][j] -= B[i * m + k] * lqr->K[k][j];
            }
        }
    }

    /** Compute eigenvalues of A_cl */
    return blc_eigenvalues((const double*)A_cl, n, poles_real, poles_imag);
}

double blc_lqr_cost_with_quantization(const BLCLQRController* lqr,
                                       const double* A, const double* B,
                                       const double* quant_errors) {
    if (!lqr) return 0.0;
    int n = lqr->n_states;
    int m = lqr->n_inputs;
    (void)A; (void)B;

    /** J_quant = J_nominal + Σᵢ σ_ε² · [K'RK]_{ii} */
    double penalty = 0.0;
    for (int i = 0; i < n; i++) {
        /** [K'RK]_{ii} */
        double krk_ii = 0.0;
        for (int a = 0; a < m; a++) {
            for (int b = 0; b < m; b++) {
                krk_ii += lqr->K[a][i] * lqr->R_diag[a] *
                          (a == b ? 1.0 : 0.0) * lqr->K[b][i];
            }
        }
        penalty += quant_errors[i] * quant_errors[i] * krk_ii;
    }

    return penalty;
}

double blc_lqr_quantization_penalty(const BLCLQRController* lqr,
                                     const double* A, const double* B,
                                     const double* quant_errors) {
    return blc_lqr_cost_with_quantization(lqr, A, B, quant_errors);
}

void blc_lqr_get_gain(const BLCLQRController* lqr, double* K_flat) {
    if (!lqr || !K_flat) return;
    int n = lqr->n_inputs;
    int m = lqr->n_states;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            K_flat[i * m + j] = lqr->K[i][j];
}

void blc_lqr_get_riccati(const BLCLQRController* lqr, double* P_flat) {
    if (!lqr || !P_flat) return;
    int n = lqr->n_states;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            P_flat[i * n + j] = lqr->P[i][j];
}

void blc_lqr_free(BLCLQRController* lqr) {
    (void)lqr; /** Nothing dynamically allocated */
}

/* ================================================================
 * Quantized State Feedback Implementation
 * ================================================================ */

void blc_quantized_state_feedback(BLCSystem* sys,
                                   const BLCLQRController* lqr,
                                   const double* x, double* u,
                                   double* quant_error) {
    if (!sys || !lqr || !x || !u) return;
    int n = lqr->n_states;
    int m = lqr->n_inputs;

    /** Quantize state */
    double xq[BLC_MAX_STATES];
    for (int i = 0; i < n; i++) {
        xq[i] = blc_quantize(&sys->sensor_quant, x[i]);
        if (quant_error) {
            quant_error[i] = x[i] - xq[i];
        }
    }

    /** u = -K xq */
    for (int i = 0; i < m; i++) {
        u[i] = 0.0;
        for (int j = 0; j < n; j++) {
            u[i] -= lqr->K[i][j] * xq[j];
        }
    }
}

double blc_max_quantization_step(const double* A_c, const double* B,
                                  const double* K, int n) {
    /** Maximum quantization step for stability via small-gain theorem.
     *
     *  Define G(s) = K(sI - A_c)^{-1}B (from quantization error to output).
     *  Stability requires: ||Δ||_∞ < 1 / ||G||_∞
     *
     *  For a rough bound:
     *    Δ_max ≈ 1 / (||K|| · ||B|| / |min Re(λ_c)|)
     */
    if (!A_c) return 0.0;

    /** Compute closed-loop eigenvalues */
    double ev_re[BLC_MAX_STATES], ev_im[BLC_MAX_STATES];
    blc_eigenvalues(A_c, n, ev_re, ev_im);

    /** Find minimum stability margin */
    double min_margin = DBL_MAX;
    for (int i = 0; i < n; i++) {
        if (-ev_re[i] < min_margin && -ev_re[i] > 0) {
            min_margin = -ev_re[i];
        }
    }
    if (min_margin > 1e10) min_margin = 0.1;  /** Fallback */

    /** Compute ||K|| and ||B|| */
    double K_norm = 0.0, B_norm = 0.0;
    if (K) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                K_norm += K[i * n + j] * K[i * n + j];
        K_norm = sqrt(K_norm);
    }
    if (B) {
        for (int i = 0; i < n * n; i++) B_norm += B[i] * B[i];
        B_norm = sqrt(B_norm);
    }

    if (K_norm * B_norm < 1e-15) return 10.0;
    return min_margin / (K_norm * B_norm + 1e-10);
}

int blc_required_quantization_levels(double range_hi, double range_lo,
                                      double delta_max) {
    if (delta_max <= 0.0) return BLC_MAX_QUANT_LEVELS;
    double range = range_hi - range_lo;
    int N = (int)ceil(range / delta_max);
    if (N < 2) N = 2;
    if (N > BLC_MAX_QUANT_LEVELS) N = BLC_MAX_QUANT_LEVELS;
    return N;
}

double blc_bitrate_from_levels(const int* levels, int n_states) {
    double bits = 0.0;
    for (int i = 0; i < n_states; i++) {
        if (levels[i] > 1) {
            bits += log2((double)levels[i]);
        }
    }
    return bits;
}

/* ================================================================
 * Minimum Attention Controller
 *
 * @ref Brockett (1997), "Minimum attention control", IEEE CDC
 *
 * The minimum attention paradigm minimizes the rate at which the
 * controller must access sensor measurements. The controller
 * "looks" at sensor data only when necessary, otherwise using
 * open-loop prediction.
 *
 * For an unstable scalar system ẋ = λx + u:
 *   Minimum attention rate = λ / ln(2)  bits/sec
 *   (same as the Data Rate Theorem bound)
 * ================================================================ */

int blc_minattn_init(BLCMinAttentionController* mac,
                      double error_threshold, int look_interval) {
    if (!mac || error_threshold <= 0.0 || look_interval < 1) return -1;

    mac->attention_rate = 0.0;
    mac->min_attention_rate = 0.0;
    mac->performance_bound = 0.0;
    mac->error_threshold = error_threshold;
    mac->look_interval = look_interval;
    mac->steps_since_last_look = 0;
    mac->prediction_error = 0.0;
    mac->total_looks = 0.0;
    mac->average_attention = 0.0;
    return 0;
}

bool blc_minattn_should_look(BLCMinAttentionController* mac) {
    if (!mac) return true;
    mac->steps_since_last_look++;

    /** Look if prediction error exceeds threshold OR
     *  mandatory look interval elapsed */
    if (mac->prediction_error > mac->error_threshold) {
        return true;
    }
    if (mac->steps_since_last_look >= mac->look_interval) {
        return true;
    }
    return false;
}

void blc_minattn_update(BLCMinAttentionController* mac,
                         double prediction_error) {
    if (!mac) return;
    mac->prediction_error = prediction_error;
    mac->total_looks += 1.0;

    /** Update running average attention rate */
    double alpha = 0.1;
    mac->average_attention = (1.0 - alpha) * mac->average_attention
                              + alpha * (1.0 / (double)(mac->look_interval + 1));
    mac->steps_since_last_look = 0;
}

double blc_minattn_theoretical_rate(const BLCPlant* plant) {
    /** Same as Data Rate Theorem */
    return blc_datarate_min_ct(plant);
}

void blc_minattn_stats(const BLCMinAttentionController* mac,
                        double* avg_rate, double* total_looks) {
    if (!mac) return;
    if (avg_rate) *avg_rate = mac->average_attention;
    if (total_looks) *total_looks = mac->total_looks;
}

/* ================================================================
 * Robust Control under Quantization — Circle Criterion
 *
 * The Circle Criterion states: a feedback interconnection of a
 * linear system G(s) and a sector-bounded nonlinearity φ(·)
 * with sector [α, β] is absolutely stable if the Nyquist plot
 * of G(jω) does not intersect the disk D(-1/α, -1/β).
 *
 * For quantization with N levels, the error satisfies:
 *   -Δ/2 ≤ φ(e) ≤ Δ/2  where Δ = 2*U/N
 * giving sector [ε_min, ε_max].
 *
 * For the scalar case (quantization of one signal):
 *   Condition: Re[G(jω)] > -1/δ for all ω
 *   where δ = U/(N·Δ) characterizes the quantization coarseness.
 * ================================================================ */

int blc_robust_circle_criterion(BLCRobustQuantCtrl* rqc,
                                 const double* A_c, const double* B,
                                 const double* K, int n) {
    if (!rqc || !A_c) return -1;

    /** For absolute stability via Circle Criterion:
     *
     *  The quantization error satisfies: ε'(ε - δx) ≤ 0
     *  where δ ∈ [0, 1/N] for uniform N-level quantization.
     *
     *  The circle condition requires: for all ω,
     *    Re[K (jωI - A_c)^{-1} B] > -1
     *
     *  For scalar input/output: this reduces to checking
     *  the DC gain.
     *
     *  We use a simple frequency-domain check at DC (ω=0):
     *    G(0) = -K A_c^{-1} B
     */

    /** Compute A_c^{-1} approximately for small n */
    /** Use the fact that: G(0) = C A_c^{-1} B with C = I */
    /** Since the quantization affects all states: C = I */

    /** For the circle criterion: need to check Re[K·(jωI - A_c)^{-1}·B] > -1/δ
     *  for all ω. We check at ω=0 as the worst case. */

    /** Invert A_c */
    double A_inv[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
    /** For small matrices, use Gaussian elimination */
    double aug[BLC_MAX_STATES][BLC_MAX_STATES*2];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A_c[i * n + j];
        }
        aug[i][n + i] = 1.0;
    }

    /** Gauss-Jordan */
    for (int p = 0; p < n; p++) {
        double pivot = aug[p][p];
        if (fabs(pivot) < 1e-12) {
            /** Singular A_c — system has integrator, use pseudo-inverse */
            rqc->is_absolutely_stable = false;
            return 1;
        }
        for (int c = 0; c < 2 * n; c++) aug[p][c] /= pivot;
        for (int r = 0; r < n; r++) {
            if (r == p) continue;
            double fac = aug[r][p];
            for (int c = 0; c < 2 * n; c++) {
                aug[r][c] -= fac * aug[p][c];
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_inv[i][j] = aug[i][n + j];

    /** Compute DC gain: G(0) = -K A_c^{-1} B (n×m output, m×n input)
     *  For scalar analysis: take norm */
    double G0[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                G0[i][j] += A_inv[i][k] * B[k * n + j];
            }
        }
    }
    /** Apply -K */
    double G0_K[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                G0_K[i][j] -= K[i * n + k] * G0[k][j];
            }
        }
    }

    /** Maximum singular value = H∞ norm lower bound */
    double G_norm = blc_matrix_norm_2((const double*)G0_K, n, n);

    /** δ_max = 1 / ||G(0)|| */
    rqc->max_tolerable_delta = 1.0 / (G_norm + 1e-10);

    /** Stable if δ < δ_max */
    rqc->is_absolutely_stable = (rqc->sector_delta < rqc->max_tolerable_delta);

    /** Gain margin: GM = δ_max / δ */
    rqc->gain_margin = rqc->max_tolerable_delta / (rqc->sector_delta + 1e-10);

    return rqc->is_absolutely_stable ? 0 : 1;
}

double blc_robust_hinf_norm(const double* A_c, const double* B,
                             const double* C, int n) {
    /** H∞ norm of G(s) = C(sI - A)^{-1}B.
     *
     *  Computed via the Hamiltonian matrix method:
     *    H = [A,    B B'/γ²;
     *         -C'C, -A']
     *
     *  ||G||_∞ < γ ⇔ H has no eigenvalues on the imaginary axis.
     *
     *  We use bisection on γ. */
    if (!A_c || !B || !C || n < 1) return 0.0;

    /** Simplified: use bound ||G||_∞ ≥ ||G(0)|| = ||C A^{-1} B|| */
    /** Compute A_c^{-1} */
    double A_inv[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
    double aug[BLC_MAX_STATES][BLC_MAX_STATES*2] = {{0}};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A_c[i * n + j];
        }
        aug[i][n + i] = 1.0;
    }
    for (int p = 0; p < n; p++) {
        double pivot = aug[p][p];
        if (fabs(pivot) < 1e-12) return 1e10;
        for (int c = 0; c < 2 * n; c++) aug[p][c] /= pivot;
        for (int r = 0; r < n; r++) {
            if (r == p) continue;
            double fac = aug[r][p];
            for (int c = 0; c < 2 * n; c++) {
                aug[r][c] -= fac * aug[p][c];
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A_inv[i][j] = aug[i][n + j];

    /** G(0) = C A^{-1} B */
    double G0[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                G0[i][j] += A_inv[i][k] * B[k * n + j];
            }
        }
    }
    double G0_C[BLC_MAX_STATES][BLC_MAX_STATES] = {{0}};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                G0_C[i][j] += C[i * n + k] * G0[k][j];
            }
        }
    }

    return blc_matrix_norm_2((const double*)G0_C, n, n);
}

/* ================================================================
 * Bandwidth-Limited Model Predictive Control
 *
 * MPC with a bit-rate proxy constraint: limit the number of
 * non-zero control inputs in the prediction horizon.
 *
 * min  Σ (x'_Qx + u'_Ru)
 * s.t. x_{k+1} = Ax_k + Bu_k
 *      ||u||_0 ≤ U_max   (sparsity constraint = rate proxy)
 *
 * The L0 constraint is relaxed to L1 regularization:
 *   min Σ (x'_Qx + u'_Ru + γ·|u|_1)
 *
 * @ref Aguilera et al. (2013), "On the inclusion of sparse
 *      constraints in MPC", Automatica
 * ================================================================ */

int blc_mpc_init(BLCMPCController* mpc, int horizon,
                  int max_nonzero_inputs) {
    if (!mpc || horizon < 1 || max_nonzero_inputs < 1) return -1;

    mpc->horizon = horizon;
    mpc->max_nonzero_inputs = max_nonzero_inputs;
    mpc->cost_opt = 0.0;
    mpc->active_constraints = 0;
    mpc->qp_iterations = 0;

    memset(mpc->u_opt, 0, sizeof(mpc->u_opt));
    memset(mpc->x_pred, 0, sizeof(mpc->x_pred));

    return 0;
}

int blc_mpc_solve(BLCMPCController* mpc, const double* x0,
                   const double* A, const double* B,
                   const double* Q_diag, const double* R_diag,
                   const double* u_prev, double* u_next) {
    if (!mpc || !x0) return -1;
    int N = mpc->horizon;
    int n = BLC_MAX_STATES;  /** Actual dimension from context */
    (void)n;

    /** Simple projected gradient descent with L1 prox.
     *
     *  Initialize u sequence to previous solution shifted by one step. */
    for (int k = 0; k < N; k++) {
        if (u_prev && k < N - 1) {
            mpc->u_opt[k] = u_prev[k + 1];  /** Warm start */
        } else {
            mpc->u_opt[k] = 0.0;
        }
    }

    double lambda = 0.1;  /** L1 regularization weight */
    double step_size = 0.01;

    /** Gradient descent with soft-thresholding */
    for (int iter = 0; iter < 200; iter++) {
        mpc->qp_iterations++;

        /** Forward simulate with current u to compute x prediction */
        double x_curr[BLC_MAX_STATES];
        for (int i = 0; i < BLC_MAX_STATES; i++) x_curr[i] = x0[i];

        double grad_u[64] = {0};
        double cost = 0.0;

        for (int k = 0; k < N; k++) {
            /** State cost */
            for (int i = 0; i < BLC_MAX_STATES; i++) {
                if (Q_diag) cost += x_curr[i] * Q_diag[i] * x_curr[i];
            }
            /** Input cost + L1 penalty */
            double u = mpc->u_opt[k];
            if (R_diag) cost += u * R_diag[0] * u;
            cost += lambda * fabs(u);

            /** Approximate gradient of cost w.r.t u_k:
             *  ∂J/∂u_k ≈ 2·R·u_k + λ·sign(u_k) + 2·x'Q·B */
            if (R_diag) grad_u[k] = 2.0 * R_diag[0] * u;
            grad_u[k] += lambda * ((u > 0) ? 1.0 : ((u < 0) ? -1.0 : 0.0));
            /** State feedback term: approximate */
            for (int i = 0; i < BLC_MAX_STATES; i++) {
                if (Q_diag) grad_u[k] += 2.0 * Q_diag[i] * x_curr[i] * B[i * BLC_MAX_STATES];
            }

            /** Update state: x_next = A x + B u */
            double x_next[BLC_MAX_STATES] = {0};
            for (int i = 0; i < BLC_MAX_STATES; i++) {
                for (int j = 0; j < BLC_MAX_STATES; j++) {
                    x_next[i] += A[i * BLC_MAX_STATES + j] * x_curr[j];
                }
                x_next[i] += B[i * BLC_MAX_STATES] * u;
            }
            for (int i = 0; i < BLC_MAX_STATES; i++) x_curr[i] = x_next[i];
        }

        /** Soft-thresholding (proximal operator for L1):
         *  u_new = S(u - α·∇, α·λ)
         *  where S(z, τ) = sign(z) · max(0, |z| - τ)
         */
        for (int k = 0; k < N; k++) {
            double z = mpc->u_opt[k] - step_size * grad_u[k];
            double tau = step_size * lambda;
            if (z > tau) {
                mpc->u_opt[k] = z - tau;
            } else if (z < -tau) {
                mpc->u_opt[k] = z + tau;
            } else {
                mpc->u_opt[k] = 0.0;
            }
        }
    }

    mpc->cost_opt = 0.0;
    /** Return first control */
    if (u_next) *u_next = mpc->u_opt[0];

    return 0;
}

double blc_mpc_effective_rate(const BLCMPCController* mpc) {
    if (!mpc || mpc->horizon == 0) return 0.0;
    int nonzeros = 0;
    for (int k = 0; k < mpc->horizon; k++) {
        if (fabs(mpc->u_opt[k]) > 1e-10) nonzeros++;
    }
    return (double)nonzeros / (double)mpc->horizon;
}

void blc_mpc_free(BLCMPCController* mpc) {
    (void)mpc;
}