#include "packet_loss_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Statistical PRNG Utilities — xorshift32 (Marsaglia 2003)
 * ============================================================================ */

unsigned long pl_xorshift32(unsigned long* state) {
    unsigned long x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return x;
}

double pl_uniform(unsigned long* state) {
    return (double)(pl_xorshift32(state) & 0x7FFFFFFF) / (double)0x80000000;
}

double pl_exponential(unsigned long* state, double lambda) {
    if (lambda <= 0.0) return 0.0;
    double u = pl_uniform(state);
    if (u < 1e-15) u = 1e-15;
    return -log(u) / lambda;
}

int pl_geometric(unsigned long* state, double p) {
    if (p <= 0.0) return 1;
    if (p >= 1.0) return 0;
    double u = pl_uniform(state);
    if (u < 1e-15) u = 1e-15;
    return (int)floor(log(u) / log(1.0 - p));
}

/* ============================================================================
 * Bernoulli Channel (i.i.d. loss)
 *
 * Model: Each packet is independently lost with probability p.
 * P(loss_k | history) = p for all k.
 *
 * This is the baseline model — most tractable analytically.
 * Key property: memoryless, so the loss sequence is a Bernoulli process.
 * ============================================================================ */

BernoulliChannel* pl_bernoulli_create(double loss_prob, unsigned long seed) {
    BernoulliChannel* ch = (BernoulliChannel*)calloc(1, sizeof(BernoulliChannel));
    if (!ch) return NULL;
    ch->loss_probability = (loss_prob < 0.0) ? 0.0 :
                           (loss_prob > 1.0) ? 1.0 : loss_prob;
    ch->seed = (seed == 0) ? 123456789UL : seed;
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;
    return ch;
}

void pl_bernoulli_free(BernoulliChannel* ch) {
    free(ch);
}

PacketStatus pl_bernoulli_transmit(BernoulliChannel* ch) {
    double u = pl_uniform(&ch->seed);
    ch->total_sent++;
    if (u < ch->loss_probability) {
        ch->total_lost++;
        ch->observed_loss_rate = (double)ch->total_lost / (double)ch->total_sent;
        return PACKET_LOST;
    }
    ch->observed_loss_rate = (double)ch->total_lost / (double)ch->total_sent;
    return PACKET_RECEIVED;
}

void pl_bernoulli_reset(BernoulliChannel* ch) {
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;
}

double pl_bernoulli_get_loss_rate(const BernoulliChannel* ch) {
    return ch->observed_loss_rate;
}

bool pl_bernoulli_is_stable(double loss_prob, double spectral_radius) {
    /* Critical loss probability: p_c = 1 - 1/ρ(A)²
     * System is mean-square stable iff p < p_c (i.e., 1-p > 1/ρ(A)²).
     *
     * Reference: Sinopoli et al. (2004), Theorem 2.
     */
    if (spectral_radius <= 1.0) {
        /* Open-loop stable system — stable for any loss rate
         * (though performance degrades with higher loss). */
        return true;
    }
    double rho2 = spectral_radius * spectral_radius;
    double p_critical = 1.0 - 1.0 / rho2;
    return loss_prob < p_critical;
}

/* ============================================================================
 * Gilbert-Elliott Channel (2-state Markov)
 *
 * Captures temporal correlations: losses tend to come in bursts.
 * The channel state evolves as a 2-state Markov chain:
 *
 *   P = | 1-p_gb    p_gb   |
 *       |  p_bg    1-p_bg  |
 *
 * Steady state: π_G = p_bg/(p_gb+p_bg), π_B = p_gb/(p_gb+p_bg)
 * ============================================================================ */

GilbertElliottChannel* pl_gilbert_elliott_create(
    double p_gb, double p_bg, double loss_good, double loss_bad,
    unsigned long seed)
{
    GilbertElliottChannel* ch = (GilbertElliottChannel*)calloc(1,
        sizeof(GilbertElliottChannel));
    if (!ch) return NULL;

    ch->p_gb = (p_gb < 0.0) ? 0.0 : (p_gb > 1.0) ? 1.0 : p_gb;
    ch->p_bg = (p_bg < 0.0) ? 0.0 : (p_bg > 1.0) ? 1.0 : p_bg;
    ch->loss_rate_good = (loss_good < 0.0) ? 0.0 : (loss_good > 1.0) ? 1.0 : loss_good;
    ch->loss_rate_bad  = (loss_bad  < 0.0) ? 0.0 : (loss_bad  > 1.0) ? 1.0 : loss_bad;
    ch->seed = (seed == 0) ? 987654321UL : seed;
    ch->current_state = 0;  /* Start in Good state */
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;

    /* Compute steady-state distribution */
    pl_gilbert_elliott_compute_steady_state(ch);

    return ch;
}

void pl_gilbert_elliott_free(GilbertElliottChannel* ch) {
    free(ch);
}

void pl_gilbert_elliott_compute_steady_state(GilbertElliottChannel* ch) {
    double denom = ch->p_gb + ch->p_bg;
    if (denom < 1e-15) {
        /* Degenerate: no transitions */
        ch->steady_p_good = (ch->current_state == 0) ? 1.0 : 0.0;
        ch->steady_p_bad  = (ch->current_state == 1) ? 1.0 : 0.0;
    } else {
        ch->steady_p_good = ch->p_bg / denom;
        ch->steady_p_bad  = ch->p_gb / denom;
    }
    ch->steady_loss_rate = ch->steady_p_good * ch->loss_rate_good
                         + ch->steady_p_bad  * ch->loss_rate_bad;
}

PacketStatus pl_gilbert_elliott_transmit(GilbertElliottChannel* ch) {
    double u;

    /* First: determine if this packet is lost in the current state */
    double loss_rate = (ch->current_state == 0) ?
                        ch->loss_rate_good : ch->loss_rate_bad;

    u = pl_uniform(&ch->seed);
    ch->total_sent++;

    PacketStatus result;
    if (u < loss_rate) {
        ch->total_lost++;
        result = PACKET_LOST;
    } else {
        result = PACKET_RECEIVED;
    }

    /* Then: update the Markov state for next transmission */
    u = pl_uniform(&ch->seed);
    double trans_prob = (ch->current_state == 0) ? ch->p_gb : ch->p_bg;
    if (u < trans_prob) {
        ch->current_state = 1 - ch->current_state;  /* Toggle state */
    }

    ch->observed_loss_rate = (double)ch->total_lost / (double)ch->total_sent;
    return result;
}

double pl_gilbert_elliott_burstiness(const GilbertElliottChannel* ch) {
    /* Burstiness = P(loss_k | loss_{k-1}) / P(loss_k)
     *
     * For GE model:
     *   P(loss_k | loss_{k-1}) ≈ weighted combination of
     *   (1-p_gb)*(1-p_bg) + transitions + loss rates in each state.
     *
     * Simplified computation using steady-state and transition structure.
     * If p_gb is small and loss_bad is large → high burstiness.
     */
    double steady_loss = ch->steady_loss_rate;
    if (steady_loss < 1e-12) return 0.0;

    /* Conditional probability: P(loss now | loss before)
     * Approximated via joint state distribution. */
    double PL_given_G = ch->loss_rate_good * (1.0 - ch->p_gb)
                       + ch->loss_rate_bad * ch->p_gb;
    double PL_given_B = ch->loss_rate_bad * (1.0 - ch->p_bg)
                       + ch->loss_rate_good * ch->p_bg;

    double conditional_PL = ch->steady_p_good * PL_given_G * ch->loss_rate_good
                           + ch->steady_p_bad  * PL_given_B * ch->loss_rate_bad;
    conditional_PL /= steady_loss;

    return conditional_PL / steady_loss;
}

void pl_gilbert_elliott_reset(GilbertElliottChannel* ch) {
    ch->current_state = 0;
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;
}

/* ============================================================================
 * K-State Markov Channel
 *
 * Generalizes Gilbert-Elliott to K ≥ 2 states.
 * Used for finer-grained channel modeling (e.g., signal quality tiers).
 *
 * Transition matrix P is K×K, row-stochastic (each row sums to 1).
 * Steady-state π solved via power iteration: π^{(t+1)} = π^{(t)} · P.
 * ============================================================================ */

MarkovChannel* pl_markov_create(int n_states, unsigned long seed) {
    if (n_states < 2) n_states = 2;

    MarkovChannel* ch = (MarkovChannel*)calloc(1, sizeof(MarkovChannel));
    if (!ch) return NULL;

    ch->n_states = n_states;
    ch->seed = (seed == 0) ? 5551212UL : seed;
    ch->current_state = 0;

    /* Allocate transition matrix (row-major) */
    ch->transition_matrix = (double**)malloc(n_states * sizeof(double*));
    double* tm_data = (double*)calloc(n_states * n_states, sizeof(double));
    for (int i = 0; i < n_states; i++) {
        ch->transition_matrix[i] = tm_data + i * n_states;
    }

    /* Initialize with uniform transitions */
    double uniform = 1.0 / (double)n_states;
    for (int i = 0; i < n_states; i++) {
        for (int j = 0; j < n_states; j++) {
            ch->transition_matrix[i][j] = uniform;
        }
    }

    ch->loss_rates = (double*)calloc(n_states, sizeof(double));
    ch->steady_state = (double*)calloc(n_states, sizeof(double));
    for (int i = 0; i < n_states; i++) {
        ch->steady_state[i] = uniform;
    }

    return ch;
}

void pl_markov_free(MarkovChannel* ch) {
    if (!ch) return;
    if (ch->transition_matrix) {
        free(ch->transition_matrix[0]);  /* contiguous allocation */
        free(ch->transition_matrix);
    }
    free(ch->loss_rates);
    free(ch->steady_state);
    free(ch);
}

void pl_markov_set_transition(MarkovChannel* ch, int from, int to, double prob) {
    if (from < 0 || from >= ch->n_states || to < 0 || to >= ch->n_states) return;
    ch->transition_matrix[from][to] = prob;
}

void pl_markov_set_loss_rate(MarkovChannel* ch, int state, double rate) {
    if (state < 0 || state >= ch->n_states) return;
    ch->loss_rates[state] = (rate < 0.0) ? 0.0 : (rate > 1.0) ? 1.0 : rate;
}

void pl_markov_compute_steady_state(MarkovChannel* ch) {
    int K = ch->n_states;
    double* pi = ch->steady_state;
    double* pi_new = (double*)malloc(K * sizeof(double));

    /* Initialize with uniform distribution */
    for (int i = 0; i < K; i++) pi[i] = 1.0 / (double)K;

    /* Power iteration: π = π · P */
    int max_iter = 10000;
    double tol = 1e-12;
    for (int iter = 0; iter < max_iter; iter++) {
        for (int j = 0; j < K; j++) {
            pi_new[j] = 0.0;
            for (int i = 0; i < K; i++) {
                pi_new[j] += pi[i] * ch->transition_matrix[i][j];
            }
        }

        double max_diff = 0.0;
        for (int i = 0; i < K; i++) {
            double diff = fabs(pi_new[i] - pi[i]);
            if (diff > max_diff) max_diff = diff;
            pi[i] = pi_new[i];
        }
        if (max_diff < tol) break;
    }

    free(pi_new);
}

PacketStatus pl_markov_transmit(MarkovChannel* ch) {
    int K = ch->n_states;

    /* Determine loss in current state */
    double loss_rate = ch->loss_rates[ch->current_state];
    double u = pl_uniform(&ch->seed);
    ch->total_sent++;

    PacketStatus result;
    if (u < loss_rate) {
        ch->total_lost++;
        result = PACKET_LOST;
    } else {
        result = PACKET_RECEIVED;
    }

    /* Transition to next state */
    u = pl_uniform(&ch->seed);
    double cumulative = 0.0;
    for (int j = 0; j < K; j++) {
        cumulative += ch->transition_matrix[ch->current_state][j];
        if (u < cumulative) {
            ch->current_state = j;
            break;
        }
    }

    ch->observed_loss_rate = (double)ch->total_lost / (double)ch->total_sent;
    return result;
}

void pl_markov_reset(MarkovChannel* ch) {
    ch->current_state = 0;
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;
}

/* ============================================================================
 * Burst Channel
 *
 * Alternating burst/gap process:
 *   In burst: each packet is lost
 *   In gap: each packet succeeds
 *   Burst length ~ Geometric(1/μ_B), Gap length ~ Geometric(1/μ_G)
 *
 * Long-run loss rate: p = μ_B / (μ_B + μ_G)
 * ============================================================================ */

BurstChannel* pl_burst_create(double mean_burst, double mean_gap, unsigned long seed) {
    BurstChannel* ch = (BurstChannel*)calloc(1, sizeof(BurstChannel));
    if (!ch) return NULL;

    ch->mean_burst_length = (mean_burst < 1.0) ? 1.0 : mean_burst;
    ch->mean_gap_length   = (mean_gap   < 1.0) ? 1.0 : mean_gap;
    ch->loss_probability = ch->mean_burst_length /
                          (ch->mean_burst_length + ch->mean_gap_length);
    ch->seed = (seed == 0) ? 314159265UL : seed;
    ch->current_burst_length = 0;
    ch->current_gap_length = 0;
    ch->in_burst = false;
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;

    return ch;
}

void pl_burst_free(BurstChannel* ch) {
    free(ch);
}

PacketStatus pl_burst_transmit(BurstChannel* ch) {
    ch->total_sent++;

    if (ch->in_burst) {
        /* In burst: this packet is lost */
        ch->current_burst_length++;
        ch->total_lost++;
        ch->observed_loss_rate = (double)ch->total_lost / (double)ch->total_sent;

        /* Check if burst should end */
        double p_end_burst = 1.0 / ch->mean_burst_length;
        if (pl_uniform(&ch->seed) < p_end_burst ||
            ch->current_burst_length >= (int)(3.0 * ch->mean_burst_length)) {
            ch->in_burst = false;
            ch->current_burst_length = 0;
            ch->current_gap_length = 0;
        }
        return PACKET_LOST;
    } else {
        /* In gap: this packet succeeds */
        ch->current_gap_length++;

        /* Check if gap should end (burst should start) */
        double p_start_burst = 1.0 / ch->mean_gap_length;
        if (pl_uniform(&ch->seed) < p_start_burst ||
            ch->current_gap_length >= (int)(3.0 * ch->mean_gap_length)) {
            ch->in_burst = true;
            ch->current_gap_length = 0;
            ch->current_burst_length = 0;
        }

        ch->observed_loss_rate = (double)ch->total_lost / (double)ch->total_sent;
        return PACKET_RECEIVED;
    }
}

void pl_burst_reset(BurstChannel* ch) {
    ch->in_burst = false;
    ch->current_burst_length = 0;
    ch->current_gap_length = 0;
    ch->total_sent = 0;
    ch->total_lost = 0;
    ch->observed_loss_rate = 0.0;
}

/* ============================================================================
 * Network Packet Lifecycle
 * ============================================================================ */

NetworkPacket* pl_packet_create(int payload_size) {
    NetworkPacket* pkt = (NetworkPacket*)calloc(1, sizeof(NetworkPacket));
    if (!pkt) return NULL;

    pkt->payload_size = (payload_size > 0) ? payload_size : 1;
    pkt->payload = (double*)calloc(pkt->payload_size, sizeof(double));
    pkt->sequence_number = 0;
    pkt->timestamp = 0.0;
    pkt->original_timestamp = 0.0;
    pkt->status = PACKET_RECEIVED;
    pkt->retransmission_count = 0;

    return pkt;
}

void pl_packet_free(NetworkPacket* pkt) {
    if (!pkt) return;
    free(pkt->payload);
    free(pkt);
}

void pl_packet_set_payload(NetworkPacket* pkt, const double* data, int size) {
    if (!pkt || !data || size <= 0) return;
    int copy_size = (size > pkt->payload_size) ? pkt->payload_size : size;
    memcpy(pkt->payload, data, copy_size * sizeof(double));
}

const double* pl_packet_get_payload(const NetworkPacket* pkt) {
    if (!pkt) return NULL;
    return pkt->payload;
}

void pl_packet_increment_retry(NetworkPacket* pkt) {
    if (!pkt) return;
    pkt->retransmission_count++;
}

/* ============================================================================
 * Packet History
 * ============================================================================ */

PacketHistory* pl_history_create(unsigned long capacity) {
    PacketHistory* hist = (PacketHistory*)calloc(1, sizeof(PacketHistory));
    if (!hist) return NULL;

    hist->capacity = (capacity > 0) ? capacity : 1024;
    hist->entries = (PacketLogEntry*)calloc(hist->capacity, sizeof(PacketLogEntry));
    hist->count = 0;

    return hist;
}

void pl_history_free(PacketHistory* hist) {
    if (!hist) return;
    free(hist->entries);
    free(hist);
}

void pl_history_record(PacketHistory* hist, unsigned long seq,
                        double time, PacketStatus st) {
    if (!hist) return;

    /* Ring-buffer behavior */
    unsigned long idx = hist->count % hist->capacity;
    hist->entries[idx].seq = seq;
    hist->entries[idx].timestamp = time;
    hist->entries[idx].status = st;
    hist->count++;
}

double pl_history_loss_rate(const PacketHistory* hist, double window) {
    if (!hist || hist->count == 0) return 0.0;

    unsigned long lookback = (unsigned long)window;
    if (lookback == 0 || lookback > hist->count) lookback = hist->count;
    if (lookback > hist->capacity) lookback = hist->capacity;

    unsigned long losses = 0;
    unsigned long start = (hist->count > lookback) ? (hist->count - lookback) : 0;

    for (unsigned long i = start; i < hist->count; i++) {
        unsigned long idx = i % hist->capacity;
        if (hist->entries[idx].status != PACKET_RECEIVED) {
            losses++;
        }
    }

    return (double)losses / (double)lookback;
}

double pl_history_burst_index(const PacketHistory* hist) {
    if (!hist || hist->count < 2) return 0.0;

    unsigned long loss_after_loss = 0;
    unsigned long total_losses = 0;
    unsigned long lookback = hist->count;
    if (lookback > hist->capacity) lookback = hist->capacity;

    unsigned long start = (hist->count > lookback) ? (hist->count - lookback + 1) : 1;
    for (unsigned long i = start; i < hist->count; i++) {
        unsigned long prev_idx = (i - 1) % hist->capacity;
        unsigned long curr_idx = i % hist->capacity;

        bool prev_lost = (hist->entries[prev_idx].status != PACKET_RECEIVED);
        bool curr_lost = (hist->entries[curr_idx].status != PACKET_RECEIVED);

        if (prev_lost) {
            total_losses++;
            if (curr_lost) loss_after_loss++;
        }
    }

    if (total_losses == 0) return 0.0;
    return (double)loss_after_loss / (double)total_losses;
}

void pl_history_print(const PacketHistory* hist) {
    if (!hist) { printf("History: NULL\n"); return; }
    printf("=== Packet History ===\n");
    printf("Entries: %lu / %lu (capacity)\n", hist->count, hist->capacity);
    printf("Recent loss rate: %.4f\n", pl_history_loss_rate(hist, 100));
    printf("Burst index: %.4f\n", pl_history_burst_index(hist));

    /* Show last 10 entries */
    unsigned long show = (hist->count < 10) ? hist->count : 10;
    unsigned long start = (hist->count > show) ? (hist->count - show) : 0;
    for (unsigned long i = start; i < hist->count; i++) {
        unsigned long idx = i % hist->capacity;
        const char* st_str = "UNKNOWN";
        switch (hist->entries[idx].status) {
            case PACKET_RECEIVED:  st_str = "OK"; break;
            case PACKET_LOST:      st_str = "LOST"; break;
            case PACKET_DELAYED:   st_str = "DELAYED"; break;
            case PACKET_CORRUPTED: st_str = "CORRUPT"; break;
        }
        printf("  seq=%6lu t=%.4f %s\n",
               hist->entries[idx].seq,
               hist->entries[idx].timestamp, st_str);
    }
}

/* ============================================================================
 * Packet Channel Generic Interface
 * ============================================================================ */

PacketChannel* pl_channel_create(ChannelModelType type) {
    PacketChannel* ch = (PacketChannel*)calloc(1, sizeof(PacketChannel));
    if (!ch) return NULL;
    ch->type = type;
    switch (type) {
        case CHANNEL_BERNOULLI:
            ch->model.bernoulli = pl_bernoulli_create(0.1, 42);
            break;
        case CHANNEL_GILBERT_ELLIOTT:
            ch->model.gilbert_elliott = pl_gilbert_elliott_create(0.1, 0.5, 0.01, 0.3, 42);
            break;
        case CHANNEL_MARKOV_K:
            ch->model.markov = pl_markov_create(3, 42);
            break;
        case CHANNEL_BURST:
            ch->model.burst = pl_burst_create(5.0, 50.0, 42);
            break;
        case CHANNEL_FADING:
        default:
            /* Create Bernoulli as fallback for unimplemented channel types */
            free(ch);
            ch = (PacketChannel*)calloc(1, sizeof(PacketChannel));
            ch->type = CHANNEL_BERNOULLI;
            ch->model.bernoulli = pl_bernoulli_create(0.1, 42);
            break;
    }
    return ch;
}

void pl_channel_free(PacketChannel* ch) {
    if (!ch) return;
    switch (ch->type) {
        case CHANNEL_BERNOULLI:
            pl_bernoulli_free(ch->model.bernoulli);
            break;
        case CHANNEL_GILBERT_ELLIOTT:
            pl_gilbert_elliott_free(ch->model.gilbert_elliott);
            break;
        case CHANNEL_MARKOV_K:
            pl_markov_free(ch->model.markov);
            break;
        case CHANNEL_BURST:
            pl_burst_free(ch->model.burst);
            break;
        case CHANNEL_FADING:
        default:
            if (ch->model.bernoulli) pl_bernoulli_free(ch->model.bernoulli);
            break;
    }
    free(ch);
}

PacketStatus pl_channel_transmit(PacketChannel* ch) {
    if (!ch) return PACKET_LOST;
    switch (ch->type) {
        case CHANNEL_BERNOULLI:
            return pl_bernoulli_transmit(ch->model.bernoulli);
        case CHANNEL_GILBERT_ELLIOTT:
            return pl_gilbert_elliott_transmit(ch->model.gilbert_elliott);
        case CHANNEL_MARKOV_K:
            return pl_markov_transmit(ch->model.markov);
        case CHANNEL_BURST:
            return pl_burst_transmit(ch->model.burst);
        case CHANNEL_FADING:
        default:
            return pl_bernoulli_transmit(ch->model.bernoulli);
    }
}

const char* pl_channel_type_name(ChannelModelType type) {
    switch (type) {
        case CHANNEL_BERNOULLI:       return "Bernoulli (i.i.d.)";
        case CHANNEL_GILBERT_ELLIOTT:  return "Gilbert-Elliott (2-Markov)";
        case CHANNEL_MARKOV_K:         return "K-State Markov";
        case CHANNEL_BURST:            return "Burst";
        case CHANNEL_FADING:           return "Fading";
        default:                       return "Unknown";
    }
}

double pl_channel_empirical_loss(const PacketChannel* ch) {
    if (!ch) return 0.0;
    switch (ch->type) {
        case CHANNEL_BERNOULLI:
            return pl_bernoulli_get_loss_rate(ch->model.bernoulli);
        case CHANNEL_GILBERT_ELLIOTT:
            return ch->model.gilbert_elliott->observed_loss_rate;
        case CHANNEL_MARKOV_K:
            return ch->model.markov->observed_loss_rate;
        case CHANNEL_BURST:
            return ch->model.burst->observed_loss_rate;
        case CHANNEL_FADING:
        default:
            return pl_bernoulli_get_loss_rate(ch->model.bernoulli);
    }
}

/* ============================================================================
 * Fading Channel (Wireless) — Rayleigh Fading Model
 *
 * Models SNR variations in wireless channels.
 * If instantaneous SNR < threshold → packet loss.
 *
 * Rayleigh fading: signal amplitude follows Rayleigh distribution.
 * SNR (dB) = SNR_mean + X where X ~ N(0, σ²) (log-normal shadowing)
 * + Rayleigh fast fading component.
 * ============================================================================ */

/* Generate correlated Rayleigh fading sample (L8: Fading Channels) */
double pl_rayleigh_fade_sample(double prev_fade, double coherence_factor,
                                unsigned long* seed) {
    /* Generate two independent Gaussian samples (Box-Muller) */
    double u1 = pl_uniform(seed);
    double u2 = pl_uniform(seed);
    if (u1 < 1e-15) u1 = 1e-15;
    double g1 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    double g2 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

    /* Rayleigh amplitude = sqrt(g1² + g2²) */
    double amplitude = sqrt(g1 * g1 + g2 * g2);

    /* Correlate with previous sample (exponential autocorrelation) */
    double new_fade = 20.0 * log10(amplitude + 1e-10);
    double correlated = coherence_factor * prev_fade
                      + (1.0 - coherence_factor) * new_fade;

    return correlated;
}

/* Note: Fading channel functions are partially implemented in this file.
 * Full Rayleigh/Rician fading with Doppler is an advanced topic (L8). */