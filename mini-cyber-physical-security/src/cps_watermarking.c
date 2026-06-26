#include "cps_watermarking.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CPS_EPS 1e-12
#define PI 3.14159265358979323846

/* ============================================================================
 * Physical Watermarking for Attack Detection (L5: Algorithms, L8: Advanced)
 *
 * Physical watermarking is a proactive defense mechanism where the
 * defender injects a known random signal (watermark) into the control
 * input. The system's response is then checked against the expected
 * response. An attacker cannot replicate the correct watermark response
 * without knowing the secret signal.
 *
 * Theorem (Mo & Sinopoli, 2010):
 *   - Without watermark: replay attack is undetectable (identical
 *     residual distribution under normal and attack conditions)
 *   - With watermark: the Kullback-Leibler divergence between the
 *     normal and attacked residual distributions is proportional
 *     to the watermark energy, enabling detection
 *
 * The trade-off: larger watermark energy improves detection but
 * degrades control performance (adds noise to the system).
 * ============================================================================ */

/* ============================================================================
 * Pseudo-random Number Generator (Linear Congruential)
 * Used for reproducible watermark sequence generation.
 * ============================================================================ */

static unsigned long lcg_state = 12345;

static void lcg_seed(unsigned long seed) {
    lcg_state = seed;
}

static double lcg_rand(void) {
    lcg_state = lcg_state * 1103515245UL + 12345UL;
    return (double)(lcg_state % 2147483648UL) / 2147483647.0;
}

static double lcg_gaussian(double mu, double sigma) {
    double u1 = lcg_rand();
    double u2 = lcg_rand();
    if (u1 < 1e-10) u1 = 1e-10;
    return mu + sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ============================================================================
 * Watermark Lifecycle
 * ============================================================================ */

void cps_watermark_init(CPSWatermark* wm, CPSWatermarkType type,
                         double amplitude, unsigned long seed) {
    if (!wm) return;
    wm->type = type;
    wm->amplitude = amplitude;
    wm->seed = seed;
    wm->capacity = 256;
    wm->length = 0;
    wm->position = 0;
    wm->energy = amplitude * amplitude;
    wm->frequency = 1.0;
    wm->signal = (double*)malloc((size_t)wm->capacity * sizeof(double));
    wm->expected_response = (double*)malloc(
        (size_t)wm->capacity * sizeof(double));
    wm->correlation_expected = 1.0;
    wm->correlation_current = 0.0;
    wm->threshold_upper = 1.5;
    wm->threshold_lower = 0.5;
    wm->response_covariance = NULL;
    wm->B_w = NULL;
    wm->secret_key = NULL;
    wm->key_length = 0;

    lcg_seed(seed);
}

void cps_watermark_free(CPSWatermark* wm) {
    if (!wm) return;
    free(wm->signal);
    free(wm->expected_response);
    free(wm->response_covariance);
    free(wm->B_w);
    free(wm->secret_key);
}

/* ============================================================================
 * Watermark Signal Generation (L5: Algorithms)
 * ============================================================================ */

double cps_watermark_next(CPSWatermark* wm) {
    if (!wm) return 0.0;

    if (wm->length >= wm->capacity) {
        int new_cap = wm->capacity * 2;
        wm->signal = (double*)realloc(wm->signal,
            (size_t)new_cap * sizeof(double));
        wm->expected_response = (double*)realloc(
            wm->expected_response,
            (size_t)new_cap * sizeof(double));
        wm->capacity = new_cap;
    }

    double val;
    switch (wm->type) {
        case CPS_WATERMARK_GAUSSIAN:
            val = lcg_gaussian(0.0, wm->amplitude);
            break;
        case CPS_WATERMARK_BINARY:
            val = (lcg_rand() > 0.5) ? wm->amplitude : -wm->amplitude;
            break;
        case CPS_WATERMARK_SINUSOIDAL:
            val = wm->amplitude * sin(2.0 * PI * wm->frequency
                                       * (double)wm->length);
            break;
        case CPS_WATERMARK_PN_SEQUENCE:
            /* Maximum-length sequence approximation */
            val = ((wm->length & (1 << (wm->length % 7))) != 0)
                  ? wm->amplitude : -wm->amplitude;
            break;
        case CPS_WATERMARK_CHAOTIC:
            /* Logistic map: x_{n+1} = r * x_n * (1 - x_n), r=3.99 */
            {
                static double chaotic_state = 0.5;
                chaotic_state = 3.99 * chaotic_state
                    * (1.0 - chaotic_state);
                val = wm->amplitude * (2.0 * chaotic_state - 1.0);
            }
            break;
        case CPS_WATERMARK_ADAPTIVE:
            val = lcg_gaussian(0.0, wm->amplitude
                * (1.0 + 0.1 * wm->length));
            break;
        default:
            val = 0.0;
            break;
    }

    wm->signal[wm->length] = val;
    wm->length++;
    return val;
}

void cps_watermark_generate_sequence(CPSWatermark* wm, int n_samples) {
    if (!wm) return;
    for (int i = 0; i < n_samples; i++)
        cps_watermark_next(wm);
}

/* ============================================================================
 * Watermark Injection and Response Analysis
 * ============================================================================ */

void cps_watermark_inject(const CPSWatermark* wm,
                           const double* u_nominal,
                           double* u_watermarked, int dim) {
    if (!wm || !u_watermarked || dim <= 0) return;
    double wm_val = (wm->length > 0)
        ? wm->signal[wm->position % wm->length] : 0.0;

    for (int i = 0; i < dim; i++) {
        u_watermarked[i] = (u_nominal ? u_nominal[i] : 0.0);
        if (wm->B_w) {
            u_watermarked[i] += wm_val * wm->B_w[i];
        } else {
            /* Default: inject equally into all inputs */
            u_watermarked[i] += wm_val / dim;
        }
    }
}

void cps_watermark_expected_response(const CPSWatermark* wm,
                                      const double* A, const double* C,
                                      int n, int p, int horizon,
                                      double* y_expected) {
    if (!wm || !A || !C || !y_expected || n <= 0 || p <= 0)
        return;

    /* y_expected = C * sum_{i=0}^{k} A^{k-i} * b_w * wm[i]
     * b_w = B_w if defined, else ones(n,1)/n */
    double* b_w = (double*)malloc((size_t)n * sizeof(double));
    if (wm->B_w) {
        for (int i = 0; i < n; i++) b_w[i] = wm->B_w[i];
    } else {
        for (int i = 0; i < n; i++) b_w[i] = 1.0 / n;
    }

    double* A_pow = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* temp = (double*)malloc((size_t)(n * n) * sizeof(double));
    double* state_contrib = (double*)calloc((size_t)n, sizeof(double));

    for (int i = 0; i < n; i++) A_pow[i * n + i] = 1.0;

    for (int k = 0; k < horizon && k < wm->length; k++) {
        /* Contribution = A^{k} * b_w * wm[position - k - 1] */
        int idx = (wm->position - k - 1 + wm->length) % wm->length;
        double wm_k = wm->signal[idx];

        /* b_w_contrib = A_pow * b_w * wm_k */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++)
                s += A_pow[i * n + j] * b_w[j];
            state_contrib[i] += s * wm_k;
        }

        /* A_pow = A_pow * A */
        cps_matrix_multiply(temp, A_pow, A, n, n, n);
        memcpy(A_pow, temp, (size_t)(n * n) * sizeof(double));
    }

    /* y_expected = C * state_contrib */
    for (int i = 0; i < p; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++)
            s += C[i * n + j] * state_contrib[j];
        y_expected[i] = s;
    }

    free(b_w); free(A_pow); free(temp); free(state_contrib);
}

/* ============================================================================
 * Watermark-Based Detection (L5: Algorithms)
 * ============================================================================ */

double cps_watermark_chi2_test(const CPSWatermark* wm,
                                const double* y_measured,
                                const double* y_expected,
                                int p) {
    if (!wm || !y_measured || !y_expected || p <= 0) return 0.0;

    double g = 0.0;
    for (int i = 0; i < p; i++) {
        double diff = y_measured[i] - y_expected[i];
        g += diff * diff;
    }
    /* Normalize by watermark energy */
    return g / (wm->energy + CPS_EPS);
}

double cps_watermark_correlation_test(CPSWatermark* wm,
                                       const double* y_residual,
                                       int horizon) {
    if (!wm || !y_residual || horizon <= 0) return 0.0;
    if (wm->length < horizon) horizon = wm->length;

    double mean_sig = 0.0, mean_res = 0.0;
    for (int i = 0; i < horizon; i++) {
        mean_sig += wm->signal[i];
        mean_res += y_residual[i];
    }
    mean_sig /= horizon;
    mean_res /= horizon;

    double cov = 0.0, var_sig = 0.0, var_res = 0.0;
    for (int i = 0; i < horizon; i++) {
        double ds = wm->signal[i] - mean_sig;
        double dr = y_residual[i] - mean_res;
        cov += ds * dr;
        var_sig += ds * ds;
        var_res += dr * dr;
    }

    if (var_sig < CPS_EPS || var_res < CPS_EPS) return 0.0;
    double rho = cov / sqrt(var_sig * var_res);

    wm->correlation_current = rho;
    return rho;
}

/* ============================================================================
 * Moving Target Defense (L8: Advanced)
 * ============================================================================ */

void cps_watermark_mtd_rotate(CPSWatermark* wm, double new_amplitude,
                               unsigned long new_seed) {
    if (!wm) return;
    wm->amplitude = new_amplitude;
    wm->energy = new_amplitude * new_amplitude;
    wm->seed = new_seed;
    wm->position = 0;
    wm->length = 0;
    lcg_seed(new_seed);
}

void cps_watermark_frequency_hop(CPSWatermark* wm, double new_freq) {
    if (!wm) return;
    wm->frequency = new_freq;
}

/* ============================================================================
 * Encrypted Watermark (L8: Advanced)
 *
 * Uses XOR-based stream cipher for watermark encryption.
 * The attacker cannot generate valid watermark without the secret key.
 *
 * encrypted_wm[i] = wm[i] XOR PRNG(key, i)
 * where PRNG is a deterministic stream derived from the key.
 * ============================================================================ */

void cps_watermark_encrypted_generate(CPSWatermark* wm,
                                       const unsigned char* key,
                                       int key_len,
                                       double* encrypted_output) {
    if (!wm || !key || key_len <= 0 || !encrypted_output) return;

    /* Generate keystream from key using simple stream cipher */
    unsigned long ks_state = 0;
    for (int i = 0; i < key_len; i++)
        ks_state = ks_state * 31 + key[i];

    for (int i = 0; i < wm->length; i++) {
        ks_state = ks_state * 1103515245UL + 12345UL;
        /* XOR watermark with keystream in mantissa bits */
        unsigned long wm_bits = *(unsigned long*)&wm->signal[i];
        wm_bits ^= ks_state;
        encrypted_output[i] = *(double*)&wm_bits;
    }
}

double cps_watermark_homomorphic_check(const double* encrypted_y,
                                        const double* encrypted_expected,
                                        int p,
                                        const double* homomorphic_key) {
    if (!encrypted_y || !encrypted_expected || p <= 0) return 0.0;

    /* Simple homomorphic comparison: correlate the encrypted values
     * (only works if encryption is additively homomorphic) */
    double corr = 0.0;
    for (int i = 0; i < p; i++)
        corr += encrypted_y[i] * encrypted_expected[i];
    return corr / (p + CPS_EPS);
}

/* ============================================================================
 * Watermark Performance Analysis (L8: Advanced)
 * ============================================================================ */

double cps_watermark_kl_detectability(const CPSWatermark* wm,
                                       double attack_snr) {
    if (!wm) return 0.0;

    /* KL divergence between N(0, sigma^2) and N(attack_snr, sigma^2)
     * D_KL = (attack_snr)^2 / (2 * sigma^2)
     * With watermark: effective sigma^2 = noise + ||C*B_w||^2 * energy */
    double signal_power = wm->energy;
    double noise_power = 1.0;
    return (attack_snr * attack_snr + signal_power)
           / (2.0 * noise_power + CPS_EPS);
}

double cps_watermark_pareto_energy(double control_cost_weight,
                                    double detection_weight,
                                    double system_sensitivity) {
    /* Pareto optimal watermark energy: balance between
     * control degradation (proportional to energy) and
     * detection improvement (proportional to sqrt(energy))
     *
     * minimize: control_cost_weight * E + detection_weight / sqrt(E)
     * optimal E* = (detection_weight / (2 * control_cost_weight))^{2/3}
     */
    if (control_cost_weight < CPS_EPS) return 1.0;
    double ratio = detection_weight
        / (2.0 * control_cost_weight * system_sensitivity);
    if (ratio < CPS_EPS) ratio = CPS_EPS;
    return pow(ratio, 2.0 / 3.0);
}

/* ============================================================================
 * Adaptive Watermarking (L8: Advanced)
 * ============================================================================ */

void cps_watermark_adapt_energy(CPSWatermark* wm,
                                 CPSSecurityState security_state,
                                 double base_energy) {
    if (!wm) return;

    double factor;
    switch (security_state) {
        case CPS_SECURE_NORMAL:
            factor = 0.5;  /* Low watermark during normal operation */
            break;
        case CPS_SECURE_SUSPICIOUS:
            factor = 1.5;  /* Increase watermark when suspicious */
            break;
        case CPS_SECURE_ATTACKED:
            factor = 2.0;  /* Maximum watermark during attack */
            break;
        case CPS_SECURE_DEGRADED:
            factor = 1.0;
            break;
        case CPS_SECURE_RECOVERING:
            factor = 1.2;
            break;
        default:
            factor = 1.0;
            break;
    }

    wm->energy = base_energy * factor;
    wm->amplitude = sqrt(wm->energy);
}
