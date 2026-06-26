#ifndef PACKET_LOSS_CORE_H
#define PACKET_LOSS_CORE_H

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

/* ============================================================================
 * Packet Loss Core — Networked Control Systems (NCS)
 *
 * Foundational types and operations for modeling and analyzing packet loss
 * in networked control systems. Based on the canonical works of:
 *   - Schenato, Sinopoli, Franceschetti, Poolla, Sastry (2007)
 *     "Foundations of Control and Estimation Over Lossy Networks"
 *   - Sinopoli, Schenato, Franceschetti, Poolla, Jordan, Sastry (2004)
 *     "Kalman Filtering with Intermittent Observations"
 *   - Hespanha, Naghshtabrizi, Xu (2007)
 *     "A Survey of Recent Results in Networked Control Systems"
 *   - Zhang, Branicky, Phillips (2001)
 *     "Stability of Networked Control Systems"
 *   - Gupta, Hassibi, Murray (2007)
 *     "Optimal LQG Control Across Packet-Dropping Links"
 * ============================================================================ */

/* --- Packet Status Enumeration --- */

/**
 * State of a single packet transmission attempt.
 * - PACKET_RECEIVED: Successfully delivered within deadline
 * - PACKET_LOST: Dropped by network (congestion, collision, fading)
 * - PACKET_DELAYED: Arrived after deadline (logically equivalent to loss)
 * - PACKET_CORRUPTED: Frame check sequence (FCS) failed at receiver
 */
typedef enum {
    PACKET_RECEIVED = 0,
    PACKET_LOST = 1,
    PACKET_DELAYED = 2,
    PACKET_CORRUPTED = 3
} PacketStatus;

/**
 * Network transport protocol abstraction.
 * - PROTO_TCP_LIKE: Acknowledged delivery with automatic retransmission.
 *   The controller knows which packets were delivered (Schenato 2007).
 * - PROTO_UDP_LIKE: Fire-and-forget, no ACKs, no retransmission.
 *   The controller does NOT know if a packet was received (Gupta 2007).
 */
typedef enum {
    PROTO_TCP_LIKE = 0,
    PROTO_UDP_LIKE = 1
} TransportProtocol;

/**
 * Hold strategy when control/sensor packet is lost.
 * - HOLD_ZERO_INPUT: Apply zero control (most conservative).
 * - HOLD_ZERO_ORDER: Continue applying last successfully transmitted input.
 * - HOLD_PREDICTIVE: Use a predictive model to estimate missing control.
 * - HOLD_LQG_OPTIMAL: LQG-optimal compensation (separation holds for TCP-like).
 */
typedef enum {
    HOLD_ZERO_INPUT = 0,
    HOLD_ZERO_ORDER = 1,
    HOLD_PREDICTIVE = 2,
    HOLD_LQG_OPTIMAL = 3
} HoldStrategy;

/**
 * Channel model type for the packet loss process.
 * - CHANNEL_BERNOULLI: i.i.d. loss probability p ∈ [0,1]
 * - CHANNEL_GILBERT_ELLIOTT: 2-state Markov chain (Good/Bad)
 * - CHANNEL_MARKOV_K: K-state finite Markov chain
 * - CHANNEL_BURST: Bursty loss with explicit run-length model
 * - CHANNEL_FADING: Wireless fading channel (Rayleigh/Rician)
 */
typedef enum {
    CHANNEL_BERNOULLI = 0,
    CHANNEL_GILBERT_ELLIOTT = 1,
    CHANNEL_MARKOV_K = 2,
    CHANNEL_BURST = 3,
    CHANNEL_FADING = 4
} ChannelModelType;

/* --- Core Data Structures --- */

/**
 * Bernoulli (i.i.d.) packet loss model.
 *
 * Each packet independently lost with probability p. This is the baseline
 * model in most theoretical work:
 *
 *   Theorem (Sinopoli et al., 2004): For Kalman filtering with intermittent
 *   observations over a Bernoulli channel with arrival probability γ = 1-p,
 *   there exists a critical value γ_c such that the expected estimation error
 *   covariance is bounded iff γ > γ_c.
 *
 * Key property: Memoryless — P(loss_n | loss_{n-1}) = P(loss_n) = p.
 */
typedef struct {
    double loss_probability;       /* p ∈ [0, 1] — i.i.d. loss probability */
    unsigned long seed;            /* PRNG seed for reproducibility */
    unsigned long total_sent;      /* Running count of transmission attempts */
    unsigned long total_lost;      /* Running count of lost packets */
    double observed_loss_rate;     /* Empirical loss rate = total_lost/total_sent */
} BernoulliChannel;

/**
 * Gilbert-Elliott (GE) channel: 2-state Markov model.
 *
 * The GE model captures temporal correlations in packet loss:
 *   State 0 (Good): Low loss probability
 *   State 1 (Bad): High loss probability (bursty)
 *
 * Transition probabilities:
 *   P(Good → Bad) = p_gb    (enter burst)
 *   P(Bad → Good)  = p_bg   (exit burst)
 *
 * Steady-state distribution:
 *   π_G = p_bg / (p_gb + p_bg),  π_B = p_gb / (p_gb + p_bg)
 *   Average loss rate = π_G·loss_good + π_B·loss_bad
 *
 * Memory: 1-step Markov — P(state_n | state_{n-1}, ..., state_0) = P(state_n | state_{n-1})
 *
 * Reference: Gilbert (1960), "Capacity of a Burst-Noise Channel"
 *            Elliott (1963), "Estimates of Error Rates for Codes on Burst-Noise Channels"
 */
typedef struct {
    double p_gb;                   /* Good→Bad transition probability */
    double p_bg;                   /* Bad→Good transition probability */
    double loss_rate_good;         /* Loss probability when in Good state */
    double loss_rate_bad;          /* Loss probability when in Bad state */
    int current_state;             /* 0 = Good, 1 = Bad */
    unsigned long seed;
    unsigned long total_sent;
    unsigned long total_lost;
    double observed_loss_rate;
    double steady_p_good;          /* Stationary probability π_G */
    double steady_p_bad;           /* Stationary probability π_B */
    double steady_loss_rate;       /* Long-run average loss rate */
} GilbertElliottChannel;

/**
 * K-state Markov chain for general packet loss processes.
 *
 * Generalizes Gilbert-Elliott to K ≥ 2 states, allowing finer modeling
 * of channel conditions (e.g., signal quality tiers, multiple interference
 * levels in wireless networks).
 *
 * The channel evolves according to a K×K stochastic matrix P
 * where P[i][j] = P(state_n+1 = j | state_n = i).
 * Each state i has an associated loss rate λ_i ∈ [0,1].
 *
 * Steady-state: π = π·P, solved via linear equations.
 */
typedef struct {
    int n_states;
    double** transition_matrix;    /* P[from][to], row-stochastic */
    double* loss_rates;            /* Loss probability in each state */
    double* steady_state;          /* Stationary distribution π */
    int current_state;
    unsigned long seed;
    unsigned long total_sent;
    unsigned long total_lost;
    double observed_loss_rate;
} MarkovChannel;

/**
 * Burst loss model: losses occur in contiguous runs.
 *
 * Characterized by mean burst length μ_B and mean gap length μ_G.
 * The process alternates between "burst" (series of losses) and
 * "gap" (series of successes).
 *
 * Long-run loss rate: p_loss = μ_B / (μ_B + μ_G)
 *
 * This model captures the empirical observation that packet losses
 * in wireless networks tend to be bursty (Arauz, 2004).
 */
typedef struct {
    double mean_burst_length;      /* Expected consecutive losses */
    double mean_gap_length;        /* Expected consecutive successful deliveries */
    double loss_probability;       /* Long-run average: μ_B/(μ_B+μ_G) */
    int current_burst_length;
    int current_gap_length;
    bool in_burst;
    unsigned long seed;
    unsigned long total_sent;
    unsigned long total_lost;
    double observed_loss_rate;
} BurstChannel;

/** Fading channel model for wireless communications. */
typedef struct {
    double snr_threshold_db;       /* SNR threshold for successful decode */
    double path_loss_exponent;     /* Path loss exponent (2=free space, 3-5=urban) */
    double doppler_frequency;      /* Maximum Doppler shift (Hz) */
    double coherence_time_ms;      /* Channel coherence time */
    double fade_margin_db;         /* Fade margin for reliability */
    unsigned long seed;
    unsigned long total_sent;
    unsigned long total_lost;
    double observed_loss_rate;
    /* Fading envelope state */
    double current_fade_db;
    double previous_fade_db;
} FadingChannel;

/** Generic channel tagged union. */
typedef struct {
    ChannelModelType type;
    union {
        BernoulliChannel* bernoulli;
        GilbertElliottChannel* gilbert_elliott;
        MarkovChannel* markov;
        BurstChannel* burst;
        FadingChannel* fading;
    } model;
} PacketChannel;

/* --- Packet Structure --- */

/**
 * A single packet in the networked control loop.
 *
 * Carries either a sensor measurement (y_k) from sensor to controller,
 * or a control signal (u_k) from controller to actuator.
 *
 * The sequence_number enables detection of out-of-order delivery
 * and duplicate detection.
 */
typedef struct {
    unsigned long sequence_number;
    double* payload;               /* Sensor reading or control signal values */
    int payload_size;
    double timestamp;              /* Transmission time */
    PacketStatus status;
    int retransmission_count;      /* TCP-like retransmission count */
    double original_timestamp;     /* First transmission time (for retransmitted) */
} NetworkPacket;

/** Single entry in the packet transmission history log. */
typedef struct {
    unsigned long seq;
    double timestamp;
    PacketStatus status;
} PacketLogEntry;

/** Complete transmission history for offline analysis. */
typedef struct {
    PacketLogEntry* entries;
    unsigned long capacity;
    unsigned long count;
} PacketHistory;

/* --- Core API --- */

/* ===== Bernoulli Channel ===== */

/**
 * Create a Bernoulli channel with i.i.d. loss probability p.
 * Complexity: O(1). Seed allows reproducible simulation.
 */
BernoulliChannel* pl_bernoulli_create(double loss_prob, unsigned long seed);

/** Free Bernoulli channel. Complexity: O(1). */
void pl_bernoulli_free(BernoulliChannel* ch);

/**
 * Simulate a packet transmission. Returns PACKET_RECEIVED or PACKET_LOST.
 * Updates internal counters. Complexity: O(1).
 */
PacketStatus pl_bernoulli_transmit(BernoulliChannel* ch);

/** Reset all counters to initial state. */
void pl_bernoulli_reset(BernoulliChannel* ch);

/** Get empirical loss rate observed so far. */
double pl_bernoulli_get_loss_rate(const BernoulliChannel* ch);

/**
 * Check if the system is mean-square stable under Bernoulli loss.
 * For a system with spectral radius ρ(A) (open-loop),
 * the critical loss probability is p_c = 1 - 1/ρ(A)^2.
 * Returns true iff p < p_c (i.e., stability is possible).
 *
 * Reference: Sinopoli et al. (2004), Theorem 2.
 */
bool pl_bernoulli_is_stable(double loss_prob, double spectral_radius);

/* ===== Gilbert-Elliott Channel ===== */

/**
 * Create a Gilbert-Elliott channel.
 * @param p_gb: Good→Bad transition probability
 * @param p_bg: Bad→Good transition probability
 * @param loss_good: Loss probability in Good state
 * @param loss_bad: Loss probability in Bad state
 * @param seed: PRNG seed
 */
GilbertElliottChannel* pl_gilbert_elliott_create(
    double p_gb, double p_bg, double loss_good, double loss_bad,
    unsigned long seed);

void pl_gilbert_elliott_free(GilbertElliottChannel* ch);

PacketStatus pl_gilbert_elliott_transmit(GilbertElliottChannel* ch);

/**
 * Compute steady-state distribution and long-run loss rate.
 * π_G = p_bg/(p_gb+p_bg), π_B = p_gb/(p_gb+p_bg).
 * Sets internal fields steady_p_good, steady_p_bad, steady_loss_rate.
 */
void pl_gilbert_elliott_compute_steady_state(GilbertElliottChannel* ch);

/**
 * Compute burstiness metric: ratio of conditional loss probability
 * P(loss|prev_loss) / P(loss), measuring temporal correlation.
 * Values > 1 indicate bursty behavior.
 */
double pl_gilbert_elliott_burstiness(const GilbertElliottChannel* ch);

void pl_gilbert_elliott_reset(GilbertElliottChannel* ch);

/* ===== K-State Markov Channel ===== */

MarkovChannel* pl_markov_create(int n_states, unsigned long seed);
void pl_markov_free(MarkovChannel* ch);

/** Set one entry in the transition matrix. Row must sum to 1 after all entries set. */
void pl_markov_set_transition(MarkovChannel* ch, int from, int to, double prob);

/** Set the loss rate for a given Markov state. */
void pl_markov_set_loss_rate(MarkovChannel* ch, int state, double rate);

/**
 * Compute stationary distribution π by solving π = π·P.
 * Uses power iteration method (iterative multiplication).
 */
void pl_markov_compute_steady_state(MarkovChannel* ch);

PacketStatus pl_markov_transmit(MarkovChannel* ch);
void pl_markov_reset(MarkovChannel* ch);

/* ===== Burst Channel ===== */

/**
 * Create a burst channel.
 * @param mean_burst: Expected number of consecutive lost packets
 * @param mean_gap: Expected number of consecutive successful deliveries
 */
BurstChannel* pl_burst_create(double mean_burst, double mean_gap, unsigned long seed);
void pl_burst_free(BurstChannel* ch);
PacketStatus pl_burst_transmit(BurstChannel* ch);
void pl_burst_reset(BurstChannel* ch);

/* ===== Packet Lifecycle ===== */

/**
 * Create a network packet with space for payload_size doubles.
 * Sequence number and timestamp should be set by caller.
 */
NetworkPacket* pl_packet_create(int payload_size);
void pl_packet_free(NetworkPacket* pkt);

/** Copy data into packet payload. Requires payload_size ≤ allocated size. */
void pl_packet_set_payload(NetworkPacket* pkt, const double* data, int size);

/** Get read-only pointer to payload data. Returns NULL if no payload. */
const double* pl_packet_get_payload(const NetworkPacket* pkt);

/** Increment retransmission counter (TCP-like behavior). */
void pl_packet_increment_retry(NetworkPacket* pkt);

/* ===== Packet History ===== */

/**
 * Create a transmission history log buffer.
 * @param capacity: Maximum entries to store (ring-buffer behavior)
 */
PacketHistory* pl_history_create(unsigned long capacity);
void pl_history_free(PacketHistory* hist);
void pl_history_record(PacketHistory* hist, unsigned long seq, double time, PacketStatus st);

/** Compute empirical loss rate over the last `window` entries. */
double pl_history_loss_rate(const PacketHistory* hist, double window);

/**
 * Compute burst index: fraction of losses that follow a previous loss.
 * Measures temporal correlation: B = P(loss_i | loss_{i-1}).
 */
double pl_history_burst_index(const PacketHistory* hist);

void pl_history_print(const PacketHistory* hist);

/* ===== Channel Generic Interface ===== */

PacketChannel* pl_channel_create(ChannelModelType type);
void pl_channel_free(PacketChannel* ch);
PacketStatus pl_channel_transmit(PacketChannel* ch);
const char* pl_channel_type_name(ChannelModelType type);
double pl_channel_empirical_loss(const PacketChannel* ch);

/* ===== Statistical Utilities ===== */

/**
 * xorshift32 PRNG (Marsaglia, 2003).
 * Fast, high-quality pseudo-random number generation.
 * Period: 2^32 - 1.
 */
unsigned long pl_xorshift32(unsigned long* state);

/** Uniform random number in [0, 1) using xorshift32. */
double pl_uniform(unsigned long* state);

/** Exponential random variate: -ln(U)/lambda (inverse CDF method). */
double pl_exponential(unsigned long* state, double lambda);

/** Geometric random variate for burst/gap length modeling. */
int pl_geometric(unsigned long* state, double p);

#endif /* PACKET_LOSS_CORE_H */