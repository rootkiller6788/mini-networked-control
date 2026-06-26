#ifndef CPS_DETECTION_H
#define CPS_DETECTION_H

#include "cps_security_core.h"

/* ============================================================================
 * CPS Attack Detection Subsystem (L2: Core Concepts, L5: Algorithms)
 *
 * Implements classical and state-of-the-art attack detectors:
 *   - Chi-squared detector (standard in FDI detection)
 *   - CUSUM detector (Page 1954, optimal for change detection)
 *   - Kalman filter residual analysis
 *   - Sequential probability ratio test
 *   - Moves detection theory from text to code
 * ============================================================================ */

/* --- Detector Lifecycle --- */

void cps_detector_init(CPSDetector* det, CPSDetectionMethod method,
                       double threshold);
void cps_detector_reset(CPSDetector* det);
void cps_detector_free_internals(CPSDetector* det);

/* --- Core Detection Algorithms (L5: Algorithms) --- */

/* Chi-squared detector: g[k] = r[k]' * Sigma^{-1} * r[k] ~ chi2(p)
 * Theorem: Under H0 (no attack), g[k] follows chi-squared(p)
 * If g[k] > chi2_{alpha}(p), reject H0 (attack detected)
 * Reference: Mo & Sinopoli (2010) */
double cps_chi2_test(CPSDetector* det, const double* residual,
                     const double* covariance_diag, int p);
bool cps_chi2_check_alarm(CPSDetector* det);

/* CUSUM (Cumulative Sum) detector
 * S[k] = max(0, S[k-1] + log(p1(y[k])/p0(y[k])))
 * Reference: Page (1954), Basseville & Nikiforov (1993) */
double cps_cusum_update(CPSDetector* det, double log_likelihood_ratio);
bool cps_cusum_check_alarm(CPSDetector* det);
void cps_cusum_reset_accumulators(CPSDetector* det);

/* Kalman filter residual computation
 * residual = y[k] - C * x_hat_predicted[k]
 * Standard for attack detection in linear Gaussian systems */
void cps_kalman_predict(double* x_pred, const double* x_est,
                        const double* A, const double* B,
                        const double* u, int n, int m);
void cps_kalman_update(double* x_est, const double* x_pred,
                        const double* y, const double* C,
                        const double* K, int n, int p);
void cps_compute_residual(double* residual, const double* y,
                          const double* C, const double* x_pred,
                          int p, int n);
/* Compute Kalman gain offline via DARE (discrete algebraic Riccati) */
int cps_kalman_gain_dare(double* K, const double* A, const double* C,
                          double Q_scale, double R_scale,
                          int n, int p, int max_iter);

/* --- Sequential Probability Ratio Test (SPRT) (L5) --- */

/* Wald's SPRT for attack detection
 * Accept H0 if lambda[k] <= A, Accept H1 if lambda[k] >= B
 * Otherwise continue sampling */
typedef struct {
    double log_lambda;          /* Log-likelihood ratio accumulator */
    double threshold_A;         /* Lower threshold (accept H0) */
    double threshold_B;         /* Upper threshold (accept H1) */
    double alpha;               /* Type I error probability */
    double beta;                /* Type II error probability */
    int decision;               /* 0=undecided, 1=H0, 2=H1 */
    int sample_count;
} CPSSPRT;

void cps_sprt_init(CPSSPRT* sprt, double alpha, double beta);
int cps_sprt_update(CPSSPRT* sprt, double log_lr);
void cps_sprt_reset(CPSSPRT* sprt);

/* --- Residual Evaluation --- */

/* Moving average of residuals */
double cps_residual_moving_average(CPSDetector* det, int window);
/* Exponentially weighted moving average */
double cps_residual_ewma(CPSDetector* det, double lambda);
/* Signal-to-noise ratio of residual */
double cps_residual_snr(const double* residual, int p, double noise_std);

/* --- Multiple Detector Fusion (L8: Advanced) --- */

/* Weighted voting across detector types */
double cps_detector_fusion_vote(CPSDetector* detectors, int n_detectors,
                                 const double* weights);
/* Bayesian posterior of attack given detections */
double cps_bayesian_attack_probability(CPSDetector* detectors,
                                        int n_detectors,
                                        double prior_attack_prob);

/* --- Performance Metrics --- */

/* Compute empirical false positive rate */
double cps_fpr_estimate(CPSDetector* det, int total_normal_steps);
/* Compute empirical detection rate */
double cps_tpr_estimate(CPSDetector* det, int total_attack_steps);
/* Average detection delay (SADD — Steady-state Average Detection Delay) */
double cps_average_detection_delay(CPSDetector* det, double* attack_times,
                                    int n_attacks);
/* ROC curve point: (fpr, tpr) for given threshold */
void cps_roc_point(CPSDetector* det, const double* normal_residuals,
                    int n_normal, const double* attack_residuals,
                    int n_attack, double* fpr_out, double* tpr_out);

/* --- Advanced: Generalized Likelihood Ratio Test (GLRT) --- */

/* GLRT for unknown attack magnitude
 * lambda = max_theta p(y|H1,theta) / p(y|H0)
 * Used when attack parameters are unknown */
double cps_glrt_statistic(const double* residuals, int n,
                           const double* covariance_diag,
                           double* theta_hat);

#endif /* CPS_DETECTION_H */
