#include "cps_detection.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CPS_EPS 1e-12

/* ============================================================================
 * Detector Lifecycle
 * ============================================================================ */

void cps_detector_init(CPSDetector* det, CPSDetectionMethod method,
                       double threshold) {
    if (!det) return;
    det->method = method;
    det->threshold = threshold;
    det->statistic = 0.0;
    det->alarm_time = 0.0;
    det->alarm_count = 0;
    det->alarm_active = 0;
    det->false_positive_rate = 0.0;
    det->detection_rate = 0.0;
    det->cusum_pos = 0.0;
    det->cusum_neg = 0.0;
    det->cusum_drift = 0.5;
    det->cusum_reset = 0.0;
    det->chi2_df = 1;
    det->chi2_pvalue = 1.0;
    det->watermark_energy = 1.0;
    det->watermark_expected = 0.0;
    det->history_length = 0;
    /* Always allocate fresh history buffer.
     * We use a sentinel capacity value to detect pre-allocated vs garbage.
     * cps_detector_free_internals sets capacity=0, so non-zero valid capacity
     * means the buffer was previously allocated through our API. */
    if (det->history_capacity > 0 && det->history_capacity <= 1000000) {
        /* Valid previously-allocated buffer — safe to free */
        free(det->residual_history);
    }
    det->history_capacity = 256;
    det->residual_history = (double*)malloc(
        (size_t)det->history_capacity * sizeof(double));
}

void cps_detector_reset(CPSDetector* det) {
    if (!det) return;
    det->statistic = 0.0;
    det->alarm_active = 0;
    det->cusum_pos = 0.0;
    det->cusum_neg = 0.0;
    det->chi2_pvalue = 1.0;
    det->history_length = 0;
}

void cps_detector_free_internals(CPSDetector* det) {
    if (!det) return;
    free(det->residual_history);
    det->residual_history = NULL;
    det->history_capacity = 0;
    det->history_length = 0;
}

/* ============================================================================
 * Chi-Squared Detector (L5: Algorithms)
 *
 * Theory: Under H0 (no attack), the residual r[k] = y[k] - C*x_hat[k]
 * is zero-mean Gaussian with covariance Sigma. The quadratic form
 * g[k] = r[k]' * Sigma^{-1} * r[k] follows a chi-squared distribution
 * with p degrees of freedom.
 *
 * Detection rule: alarm if g[k] > chi2_{1-alpha}(p)
 * where alpha is the desired false-positive rate (typically 0.05).
 *
 * For p=1: chi2_{0.95}(1) = 3.841
 * For p=2: chi2_{0.95}(2) = 5.991
 * For p=3: chi2_{0.95}(3) = 7.815
 *
 * Reference: Mo & Sinopoli (2010), "False Data Injection Attacks
 *            in Control Systems", Section III-B
 * ============================================================================ */

double cps_chi2_test(CPSDetector* det, const double* residual,
                     const double* covariance_diag, int p) {
    if (!det || !residual || p <= 0) return 0.0;
    double g = 0.0;
    for (int i = 0; i < p; i++) {
        double sigma_sq = covariance_diag
            ? fabs(covariance_diag[i]) + CPS_EPS
            : 1.0;
        g += (residual[i] * residual[i]) / sigma_sq;
    }
    det->statistic = g;
    det->chi2_df = p;

    /* Approximate p-value using Wilson-Hilferty transformation
     * (chi2 -> normal approximation for p >= 3) */
    if (p == 1) {
        det->chi2_pvalue = 1.0 - erf(sqrt(g) / sqrt(2.0));
    } else if (p == 2) {
        det->chi2_pvalue = exp(-g / 2.0);
    } else {
        double z = (pow(g / p, 1.0/3.0) - (1.0 - 2.0/(9.0*p)))
                   / sqrt(2.0/(9.0*p));
        det->chi2_pvalue = 0.5 * (1.0 - erf(z / sqrt(2.0)));
    }

    /* Store residual history */
    if (det->history_length < det->history_capacity) {
        det->residual_history[det->history_length] = g;
        det->history_length++;
    }

    cps_chi2_check_alarm(det);
    return g;
}

bool cps_chi2_check_alarm(CPSDetector* det) {
    if (!det) return false;
    if (det->statistic > det->threshold) {
        det->alarm_active = 1;
        det->alarm_count++;
        return true;
    }
    det->alarm_active = 0;
    return false;
}

/* ============================================================================
 * CUSUM (Cumulative Sum) Detector (L5: Algorithms)
 *
 * CUSUM is optimal for detecting a change in the mean of a sequence
 * from mu_0 to mu_1 (Page, 1954).
 *
 * Recursion:
 *   S_pos[k] = max(0, S_pos[k-1] + y[k] - mu_0 - drift/2)
 *   S_neg[k] = max(0, S_neg[k-1] - y[k] + mu_0 - drift/2)
 *
 * Alarm when S_pos > threshold or S_neg > threshold.
 *
 * The CUSUM drift parameter (delta = mu_1 - mu_0) controls
 * sensitivity to the expected attack magnitude. Smaller drift
 * detects smaller attacks but increases false alarms.
 *
 * Reference: Basseville & Nikiforov (1993),
 *            "Detection of Abrupt Changes: Theory and Application"
 * ============================================================================ */

double cps_cusum_update(CPSDetector* det, double log_likelihood_ratio) {
    if (!det) return 0.0;

    /* Standard CUSUM on log-likelihood ratio */
    det->cusum_pos += log_likelihood_ratio - det->cusum_drift;
    if (det->cusum_pos < det->cusum_reset)
        det->cusum_pos = det->cusum_reset;

    det->cusum_neg -= log_likelihood_ratio + det->cusum_drift;
    if (det->cusum_neg < det->cusum_reset)
        det->cusum_neg = det->cusum_reset;

    det->statistic = (det->cusum_pos > det->cusum_neg)
                     ? det->cusum_pos : det->cusum_neg;

    if (det->history_length < det->history_capacity) {
        det->residual_history[det->history_length] = det->statistic;
        det->history_length++;
    }

    cps_cusum_check_alarm(det);
    return det->statistic;
}

bool cps_cusum_check_alarm(CPSDetector* det) {
    if (!det) return false;
    if (det->statistic > det->threshold) {
        det->alarm_active = 1;
        det->alarm_count++;
        return true;
    }
    det->alarm_active = 0;
    return false;
}

void cps_cusum_reset_accumulators(CPSDetector* det) {
    if (!det) return;
    det->cusum_pos = det->cusum_reset;
    det->cusum_neg = det->cusum_reset;
    det->statistic = 0.0;
}

/* ============================================================================
 * Kalman Filter Residual Computation (L5: Algorithms)
 *
 * The Kalman filter provides optimal state estimation for linear
 * Gaussian systems. The innovation (residual) sequence under normal
 * operation is zero-mean white Gaussian, making it the standard
 * basis for anomaly detection in CPS.
 *
 * Predict step:
 *   x_pred[k] = A * x_est[k-1] + B * u[k-1]
 *   P_pred[k] = A * P_est[k-1] * A' + Q
 *
 * Update step:
 *   K[k] = P_pred * C' * (C * P_pred * C' + R)^{-1}
 *   x_est[k] = x_pred[k] + K[k] * (y[k] - C * x_pred[k])
 *   P_est[k] = (I - K[k] * C) * P_pred[k]
 *
 * The residual r[k] = y[k] - C * x_pred[k] is the "innovation"
 * ============================================================================ */

void cps_kalman_predict(double* x_pred, const double* x_est,
                        const double* A, const double* B,
                        const double* u, int n, int m) {
    /* x_pred = A * x_est + B * u */
    double* Ax = (double*)malloc((size_t)n * sizeof(double));
    double* Bu = (double*)calloc((size_t)n, sizeof(double));

    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int k = 0; k < n; k++)
            s += A[i * n + k] * x_est[k];
        Ax[i] = s;
    }
    if (u && m > 0) {
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int k = 0; k < m; k++)
                s += B[i * m + k] * u[k];
            Bu[i] = s;
        }
    }
    for (int i = 0; i < n; i++)
        x_pred[i] = Ax[i] + Bu[i];

    free(Ax); free(Bu);
}

void cps_kalman_update(double* x_est, const double* x_pred,
                        const double* y, const double* C,
                        const double* K, int n, int p) {
    /* innovation = y - C * x_pred */
    double* innov = (double*)malloc((size_t)p * sizeof(double));
    double* Cx = (double*)malloc((size_t)p * sizeof(double));

    for (int i = 0; i < p; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++)
            s += C[i * n + j] * x_pred[j];
        Cx[i] = s;
        innov[i] = y[i] - Cx[i];
    }

    /* x_est = x_pred + K * innovation */
    for (int i = 0; i < n; i++) {
        double correction = 0.0;
        for (int j = 0; j < p; j++)
            correction += K[i * p + j] * innov[j];
        x_est[i] = x_pred[i] + correction;
    }

    free(innov); free(Cx);
}

void cps_compute_residual(double* residual, const double* y,
                          const double* C, const double* x_pred,
                          int p, int n) {
    for (int i = 0; i < p; i++) {
        double y_hat = 0.0;
        for (int j = 0; j < n; j++)
            y_hat += C[i * n + j] * x_pred[j];
        residual[i] = y[i] - y_hat;
    }
}

/* ============================================================================
 * Kalman Gain via DARE Iteration (L5: Algorithms)
 *
 * Solves the Discrete Algebraic Riccati Equation (DARE):
 *   P = A*P*A' - A*P*C'*(C*P*C' + R)^{-1}*C*P*A' + Q
 *
 * Uses iterative method (pushing P forward in time) until convergence.
 * This is the standard method when dlqr/dare from control toolboxes
 * is not available. Convergence is guaranteed for detectable systems.
 *
 * The steady-state Kalman gain is:
 *   K = P * C' * (C * P * C' + R)^{-1}
 *
 * Returns number of iterations until convergence (or max_iter if not).
 * ============================================================================ */

int cps_kalman_gain_dare(double* K, const double* A, const double* C,
                          double Q_scale, double R_scale,
                          int n, int p, int max_iter) {
    if (n <= 0 || p <= 0 || !K || !A || !C) return -1;

    double* P = (double*)calloc((size_t)(n * n), sizeof(double));
    double* P_next = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* CT = (double*)malloc((size_t)(n * p) * sizeof(double));
    double* CPC_R = (double*)malloc((size_t)(p * p) * sizeof(double));
    double* APC = (double*)malloc((size_t)(n * p) * sizeof(double));
    double* temp_nn = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp_pn = (double*)malloc((size_t)(p * n) * sizeof(double));

    /* Initialize P = Q */
    for (int i = 0; i < n; i++) P[i * n + i] = Q_scale;

    /* Transpose C */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            CT[i * p + j] = C[j * n + i];

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* C * P * C' */
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    double c_p = 0.0;
                    for (int l = 0; l < n; l++)
                        c_p += C[i * n + l] * P[l * n + k];
                    s += c_p * CT[k * p + j];
                }
                CPC_R[i * p + j] = s;
            }
            CPC_R[i * p + i] += R_scale;
        }

        /* A * P * C' */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    double p_c = 0.0;
                    for (int l = 0; l < n; l++)
                        p_c += P[i * n + l] * CT[l * p + j];
                    s += A[i * n + k] * p_c;
                }
                APC[i * p + j] = s;
            }
        }

        /* Kalman gain: K = APC * (CPC_R)^{-1} */
        if (p == 1) {
            double inv = 1.0 / (CPC_R[0] + CPS_EPS);
            for (int i = 0; i < n; i++) K[i] = APC[i] * inv;
        } else if (p == 2) {
            double det = CPC_R[0]*CPC_R[3] - CPC_R[1]*CPC_R[2];
            double inv_det = 1.0 / (fabs(det) + CPS_EPS);
            for (int i = 0; i < n; i++) {
                K[i*2]   = (APC[i*2]*CPC_R[3] - APC[i*2+1]*CPC_R[2]) * inv_det;
                K[i*2+1] = (APC[i*2+1]*CPC_R[0] - APC[i*2]*CPC_R[1]) * inv_det;
            }
        } else {
            /* Diagonal approximation */
            for (int i = 0; i < n; i++)
                for (int j = 0; j < p; j++)
                    K[i * p + j] = APC[i * p + j]
                        / (CPC_R[j * p + j] + CPS_EPS);
        }

        /* P_next = A*P*A' - K*(C*P*A') + Q */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) {
                    double pa = 0.0;
                    for (int l = 0; l < n; l++)
                        pa += P[i * n + l] * A[j * n + l]; /* A' at (l,j) */
                    s += A[i * n + k] * pa;
                }
                temp_nn[i * n + j] = s;
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double kcpa = 0.0;
                for (int l = 0; l < p; l++) {
                    double cpa = 0.0;
                    for (int m = 0; m < n; m++)
                        cpa += C[l * n + m] * temp_nn[m * n + j];
                    kcpa += K[i * p + l] * cpa;
                }
                P_next[i * n + j] = temp_nn[i * n + j] - kcpa;
                if (i == j) P_next[i * n + j] += Q_scale;
            }
        }

        /* Check convergence: ||P_next - P||_F / ||P||_F */
        double norm_diff = 0.0, norm_P = 0.0;
        for (int i = 0; i < n * n; i++) {
            double d = P_next[i] - P[i];
            norm_diff += d * d;
            norm_P += P[i] * P[i];
        }
        if (sqrt(norm_diff) < 1e-8 * (sqrt(norm_P) + 1.0))
            break;

        memcpy(P, P_next, (size_t)(n * n) * sizeof(double));
    }

    free(P); free(P_next); free(CT); free(CPC_R);
    free(APC); free(temp_nn); free(temp_pn);
    return iter < max_iter ? iter + 1 : max_iter;
}

/* ============================================================================
 * Sequential Probability Ratio Test (SPRT) (L5: Algorithms)
 *
 * Wald's SPRT (1945) is optimal for sequential hypothesis testing
 * in the sense that it minimizes the expected sample size for
 * given Type I (alpha) and Type II (beta) error probabilities.
 *
 * Log-likelihood ratio: Lambda[k] = sum_{i=1}^{k} log(p1(y_i)/p0(y_i))
 *
 * Thresholds:
 *   A = log(beta/(1-alpha))  -- accept H0 if Lambda <= A
 *   B = log((1-beta)/alpha)  -- accept H1 if Lambda >= B
 *
 * For CPS security:
 *   H0: no attack (normal operation)
 *   H1: attack present
 * ============================================================================ */

void cps_sprt_init(CPSSPRT* sprt, double alpha, double beta) {
    if (!sprt) return;
    sprt->alpha = alpha;
    sprt->beta = beta;
    sprt->threshold_A = log(beta / (1.0 - alpha + CPS_EPS));
    sprt->threshold_B = log((1.0 - beta + CPS_EPS) / alpha);
    sprt->log_lambda = 0.0;
    sprt->decision = 0;
    sprt->sample_count = 0;
}

int cps_sprt_update(CPSSPRT* sprt, double log_lr) {
    if (!sprt) return 0;
    sprt->log_lambda += log_lr;
    sprt->sample_count++;

    if (sprt->log_lambda <= sprt->threshold_A) {
        sprt->decision = 1;  /* Accept H0: no attack */
        return 1;
    }
    if (sprt->log_lambda >= sprt->threshold_B) {
        sprt->decision = 2;  /* Accept H1: attack detected */
        return 2;
    }
    sprt->decision = 0;  /* Continue sampling */
    return 0;
}

void cps_sprt_reset(CPSSPRT* sprt) {
    if (!sprt) return;
    sprt->log_lambda = 0.0;
    sprt->decision = 0;
    sprt->sample_count = 0;
}

/* ============================================================================
 * Residual Evaluation Functions (L5: Algorithms)
 *
 * Moving average and EWMA provide smoothed estimates of the residual
 * trend, reducing noise variance at the cost of detection delay.
 * SNR quantifies how detectable the attack is relative to noise.
 * ============================================================================ */

double cps_residual_moving_average(CPSDetector* det, int window) {
    if (!det || window <= 0 || det->history_length < window)
        return 0.0;

    double sum = 0.0;
    int start = det->history_length - window;
    for (int i = start; i < det->history_length; i++)
        sum += det->residual_history[i];
    return sum / window;
}

double cps_residual_ewma(CPSDetector* det, double lambda) {
    if (!det || det->history_length == 0) return 0.0;
    if (lambda < 0.0) lambda = 0.0;
    if (lambda > 1.0) lambda = 1.0;

    double ewma = det->residual_history[0];
    for (int i = 1; i < det->history_length; i++)
        ewma = lambda * det->residual_history[i] + (1.0 - lambda) * ewma;
    return ewma;
}

double cps_residual_snr(const double* residual, int p, double noise_std) {
    if (!residual || p <= 0) return 0.0;
    double signal_power = 0.0;
    for (int i = 0; i < p; i++)
        signal_power += residual[i] * residual[i];
    signal_power /= p;
    double noise_power = noise_std * noise_std;
    if (noise_power < CPS_EPS) return signal_power / CPS_EPS;
    return signal_power / noise_power;
}

/* ============================================================================
 * Multiple Detector Fusion (L8: Advanced)
 * ============================================================================ */

double cps_detector_fusion_vote(CPSDetector* detectors, int n_detectors,
                                 const double* weights) {
    if (!detectors || n_detectors <= 0) return 0.0;
    double weighted_sum = 0.0;
    double weight_total = 0.0;

    for (int i = 0; i < n_detectors; i++) {
        double w = weights ? weights[i] : 1.0;
        double vote = detectors[i].alarm_active ? 1.0 : 0.0;
        weighted_sum += w * vote;
        weight_total += w;
    }
    return (weight_total > CPS_EPS)
           ? weighted_sum / weight_total : 0.0;
}

double cps_bayesian_attack_probability(CPSDetector* detectors,
                                        int n_detectors,
                                        double prior_attack_prob) {
    if (!detectors || n_detectors <= 0) return prior_attack_prob;

    /* P(attack | detections) =
     *   prior * P(detections|attack) / P(detections) */
    double log_posterior = log(prior_attack_prob + CPS_EPS);
    double log_prior_normal = log(1.0 - prior_attack_prob + CPS_EPS);

    for (int i = 0; i < n_detectors; i++) {
        double tpr_i = detectors[i].detection_rate;
        double fpr_i = detectors[i].false_positive_rate;
        if (tpr_i < CPS_EPS) tpr_i = 0.5;
        if (fpr_i < CPS_EPS) fpr_i = 0.01;

        if (detectors[i].alarm_active) {
            log_posterior += log(tpr_i);
            log_prior_normal += log(fpr_i);
        } else {
            log_posterior += log(1.0 - tpr_i);
            log_prior_normal += log(1.0 - fpr_i);
        }
    }

    double posterior_odds = exp(log_posterior - log_prior_normal);
    return posterior_odds / (1.0 + posterior_odds);
}

/* ============================================================================
 * Performance Metrics (L2: Core Concepts)
 * ============================================================================ */

double cps_fpr_estimate(CPSDetector* det, int total_normal_steps) {
    if (!det || total_normal_steps <= 0) return 0.0;
    /* Estimate FPR from alarm count during known-normal periods */
    return (double)det->alarm_count / (double)total_normal_steps;
}

double cps_tpr_estimate(CPSDetector* det, int total_attack_steps) {
    if (!det || total_attack_steps <= 0) return 0.0;
    return (double)det->alarm_count / (double)total_attack_steps;
}

double cps_average_detection_delay(CPSDetector* det,
                                    double* attack_times,
                                    int n_attacks) {
    if (!det || !attack_times || n_attacks <= 0) return 0.0;
    /* Simple estimate: assuming detection time is alarm_time */
    double total_delay = 0.0;
    for (int i = 0; i < n_attacks; i++) {
        double delay = det->alarm_time - attack_times[i];
        if (delay > 0) total_delay += delay;
    }
    return total_delay / n_attacks;
}

void cps_roc_point(CPSDetector* det,
                    const double* normal_residuals, int n_normal,
                    const double* attack_residuals, int n_attack,
                    double* fpr_out, double* tpr_out) {
    if (!det || !fpr_out || !tpr_out) return;

    int fp = 0, tp = 0;
    double thresh = det->threshold;

    if (normal_residuals) {
        for (int i = 0; i < n_normal; i++)
            if (normal_residuals[i] > thresh) fp++;
    }
    if (attack_residuals) {
        for (int i = 0; i < n_attack; i++)
            if (attack_residuals[i] > thresh) tp++;
    }

    *fpr_out = (n_normal > 0) ? (double)fp / n_normal : 0.0;
    *tpr_out = (n_attack > 0) ? (double)tp / n_attack : 0.0;
}

/* ============================================================================
 * Generalized Likelihood Ratio Test (GLRT) (L8: Advanced)
 *
 * GLRT handles the case where attack parameters are unknown:
 *   lambda = max_{theta} p(y | H1, theta) / p(y | H0)
 *
 * For Gaussian observations with unknown mean change:
 *   theta_hat = mean(residuals)
 *   GLRT statistic = n * (theta_hat^2 / sigma^2)
 *
 * Reference: Kay (1998), "Fundamentals of Statistical Signal Processing,
 *            Vol. II: Detection Theory"
 * ============================================================================ */

double cps_glrt_statistic(const double* residuals, int n,
                           const double* covariance_diag,
                           double* theta_hat) {
    if (!residuals || n <= 0) return 0.0;

    /* Estimate unknown mean shift */
    double mean_residual = 0.0;
    for (int i = 0; i < n; i++)
        mean_residual += residuals[i];
    mean_residual /= n;

    if (theta_hat) *theta_hat = mean_residual;

    /* Average variance */
    double avg_var = 1.0;
    if (covariance_diag) {
        avg_var = 0.0;
        for (int i = 0; i < n; i++)
            avg_var += covariance_diag[i % n];
        avg_var /= n;
    }
    if (avg_var < CPS_EPS) avg_var = CPS_EPS;

    /* GLRT statistic: n * (theta_hat^2 / sigma^2) ~ chi2(1) under H0 */
    return n * (mean_residual * mean_residual) / avg_var;
}
