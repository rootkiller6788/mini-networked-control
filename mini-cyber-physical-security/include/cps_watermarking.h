#ifndef CPS_WATERMARKING_H
#define CPS_WATERMARKING_H

#include "cps_security_core.h"

/* ============================================================================
 * CPS Physical Watermarking for Attack Detection (L5: Algorithms, L8: Advanced)
 *
 * Reference: Mo & Sinopoli (2010) — "Physical Watermarking for Replay Attack
 *            Detection in Cyber-Physical Systems"
 * Reference: Ferrari & Teixeira (2021) — "Detection of Cyber Attacks"
 *
 * Physical watermarking injects a known, random signal into the control
 * input. The response of the system to this watermark is then checked
 * against the expected response — an attacker cannot replicate the
 * watermark effect without knowledge of the secret signal.
 * ============================================================================ */

/* --- Watermark Signal Types --- */

typedef enum {
    CPS_WATERMARK_GAUSSIAN = 0,  /* Gaussian i.i.d. watermark */
    CPS_WATERMARK_BINARY = 1,    /* Binary {+Delta, -Delta} watermark */
    CPS_WATERMARK_SINUSOIDAL = 2,/* Sinusoidal probing signal */
    CPS_WATERMARK_PN_SEQUENCE = 3,/* Pseudo-noise (max-length sequence) */
    CPS_WATERMARK_CHAOTIC = 4,   /* Chaotic oscillator watermark */
    CPS_WATERMARK_ADAPTIVE = 5   /* Adaptive watermark (energy based on state) */
} CPSWatermarkType;

typedef struct {
    CPSWatermarkType type;
    double* signal;              /* The watermark signal values */
    double* expected_response;   /* Expected system output response */
    int length;
    int capacity;
    int position;                /* Current position in the watermark cycle */
    double energy;               /* Watermark energy = var(signal) */
    double amplitude;            /* Amplitude of watermark */
    double frequency;            /* For sinusoidal */
    unsigned long seed;          /* RNG seed for reproducibility */
    double* secret_key;          /* Secret key for PN sequence generation */
    int key_length;

    /* Watermark injection dynamics */
    double* B_w;                 /* Watermark input matrix (n × 1) */

    /* Detection statistics */
    double correlation_expected; /* Expected correlation wm→response */
    double correlation_current;  /* Current correlation estimate */
    double threshold_upper;      /* Upper detection threshold */
    double threshold_lower;      /* Lower detection threshold */

    /* Covariance of watermark response (for chi2 test) */
    double* response_covariance; /* Covariance matrix of Δy[wmk] */
} CPSWatermark;

/* --- Watermark Lifecycle --- */

void cps_watermark_init(CPSWatermark* wm, CPSWatermarkType type,
                         double amplitude, unsigned long seed);
void cps_watermark_free(CPSWatermark* wm);

/* Generate next watermark value to inject */
double cps_watermark_next(CPSWatermark* wm);

/* Generate full watermark sequence of given length */
void cps_watermark_generate_sequence(CPSWatermark* wm, int n_samples);

/* --- Watermark Injection and Response Analysis --- */

/* Inject watermark: u[k] = u_nominal[k] + Delta_u[k]
 * where Delta_u[k] = B_w * wm[k] (the physical watermark) */
void cps_watermark_inject(const CPSWatermark* wm,
                           const double* u_nominal,
                           double* u_watermarked, int dim);

/* Compute expected system response to watermark
 * y_expected[k] = C * sum_{i=0}^{k} A^{k-i} * B_w * wm[i]
 * This is what the legitimate system should output */
void cps_watermark_expected_response(const CPSWatermark* wm,
                                      const double* A, const double* C,
                                      int n, int p, int horizon,
                                      double* y_expected);

/* --- Watermark-Based Detection (L5: Algorithms) --- */

/* Chi-squared test on watermark-induced measurement difference
 * g_wm = (y - y_expected)' * Sigma^{-1} * (y - y_expected)
 * Under no-attack: g_wm ~ chi2(p). Under replay: g_wm follows
 * non-central chi2 with non-centrality = watermark energy / noise */
double cps_watermark_chi2_test(const CPSWatermark* wm,
                                const double* y_measured,
                                const double* y_expected,
                                int p);

/* Correlation-based watermark check
 * rho = corr(Delta_y, Delta_y_expected)
 * Under no attack: rho ≈ 1. Under attack: rho ≈ 0 (replay) or biased */
double cps_watermark_correlation_test(CPSWatermark* wm,
                                       const double* y_residual,
                                       int horizon);

/* --- Moving Target Defense via Time-Varying Watermark (L8: Advanced) --- */

/* Change watermark parameters to prevent attacker adaptation
 * This is the core of Moving Target Defense (MTD) for CPS */
void cps_watermark_mtd_rotate(CPSWatermark* wm, double new_amplitude,
                               unsigned long new_seed);

/* Frequency hopping: change sinusoidal watermark frequency */
void cps_watermark_frequency_hop(CPSWatermark* wm, double new_freq);

/* --- Advanced: Encrypted Watermark (L8: Advanced) --- */

/* Generate encrypted watermark using stream cipher (simple XOR-based)
 * The attacker cannot generate valid encrypted watermark without the key */
void cps_watermark_encrypted_generate(CPSWatermark* wm,
                                       const unsigned char* key,
                                       int key_len,
                                       double* encrypted_output);

/* Homomorphic verification: check watermark response without
 * decrypting measurements (for privacy-preserving detection) */
double cps_watermark_homomorphic_check(const double* encrypted_y,
                                        const double* encrypted_expected,
                                        int p,
                                        const double* homomorphic_key);

/* --- Watermark Performance Analysis --- */

/* Kullback-Leibler divergence between distributions with/without attack
 * D_KL(P_detection|no_attack || P_detection|attack)
 * Higher KL → better detectability */
double cps_watermark_kl_detectability(const CPSWatermark* wm,
                                       double attack_snr);

/* Trade-off: control performance degradation vs detection improvement
 * Returns Pareto-optimal watermark energy for given weight w */
double cps_watermark_pareto_energy(double control_cost_weight,
                                    double detection_weight,
                                    double system_sensitivity);

/* --- Adaptive Watermarking (L8: Advanced) --- */

/* Adjust watermark energy based on current security state
 * Higher security risk → larger watermark energy */
void cps_watermark_adapt_energy(CPSWatermark* wm,
                                 CPSSecurityState security_state,
                                 double base_energy);

#endif /* CPS_WATERMARKING_H */
