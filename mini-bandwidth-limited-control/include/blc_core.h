/**
 * blc_core.h — Bandwidth-Limited Control: Core Types and Definitions
 *
 * This module implements the theory of control systems operating under
 * communication bandwidth constraints, as developed in the seminal works
 * of Wong & Brockett (1997, 1999), Nair & Evans (2000, 2004), and
 * Tatikonda & Mitter (2004).
 *
 * Key references:
 *   - Wong & Brockett, "Systems with finite communication bandwidth
 *     constraints" IEEE TAC, 1997, 1999
 *   - Nair & Evans, "Stabilizability of stochastic linear systems with
 *     finite feedback data rates" SICON, 2004
 *   - Tatikonda & Mitter, "Control under communication constraints"
 *     IEEE TAC, 2004
 *   - Baillieul & Antsaklis, "Control and communication challenges in
 *     networked real-time systems" Proc IEEE, 2007
 *
 * Knowledge coverage: L1 (Definitions), L2 (Core Concepts), L3 (Math Structures)
 */

#ifndef BLC_CORE_H
#define BLC_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>

/* ================================================================
 * L1: Core Definitions
 * ================================================================ */

/** Shannon-Hartley channel capacity: C = B * log2(1 + SNR)
 *  @ref Shannon (1948), "A Mathematical Theory of Communication"
 */
#define BLC_MAX_CHANNELS       128
#define BLC_MAX_STATES          32
#define BLC_MAX_PACKET_BYTES   256
#define BLC_MAX_QUANT_LEVELS  4096
#define BLC_MAX_SCHEDULE_SLOTS 64

/** Bandwidth (Hz) — the frequency range available for communication.
 *  In control systems, bandwidth limits constrain the rate at which
 *  sensor measurements and control commands can be exchanged.
 */
typedef struct {
    double bandwidth_hz;         /** Available bandwidth in Hz */
    double snr;                  /** Signal-to-noise ratio (linear, not dB) */
    double capacity_bps;         /** Derived: C = B * log2(1+SNR) bits/sec */
    double latency_ms;           /** One-way communication latency */
    double packet_loss_rate;     /** Packet drop probability [0, 1] */
    double jitter_ms;            /** Delay variation (standard deviation) */
} BLCChannel;

/** Quantizer: maps continuous values to discrete levels.
 *  For a uniform quantizer with q levels over [-U, U]:
 *    quantization step Δ = 2*U / q
 *    quantization error bound |e| ≤ Δ/2 = U/q
 *
 *  Quantization is the fundamental link between continuous control
 *  signals and the discrete, band-limited communication channel.
 */
typedef struct {
    int     levels;              /** Number of quantization levels q */
    double  range_lo;            /** Lower bound of quantizer range */
    double  range_hi;            /** Upper bound of quantizer range */
    double  step;                /** Quantization step Δ = (hi - lo) / levels */
    double  max_error;           /** Maximum quantization error Δ/2 */
    double  overload_prob;       /** Probability of state outside range */
    int     overload_count;      /** Number of overload events observed */
    bool    is_logarithmic;      /** Logarithmic vs. uniform quantization */
    double  log_base;            /** Base for logarithmic quantizer (e.g., 2) */
} BLCQuantizer;

/** Bandwidth-limited plant model.
 *  Continuous-time linear system:  ẋ = Ax + Bu + w
 *  with sensor measurements quantized before transmission.
 *
 *  The fundamental question: what is the minimum bit rate R
 *  required to stabilize this system?
 */
typedef struct {
    int      n_states;                  /** State dimension n */
    int      n_inputs;                  /** Input dimension m */
    int      n_outputs;                 /** Output dimension p */
    double   A[BLC_MAX_STATES][BLC_MAX_STATES]; /** System matrix A ∈ R^{n×n} */
    double   B[BLC_MAX_STATES][BLC_MAX_STATES]; /** Input matrix B ∈ R^{n×m} */
    double   C[BLC_MAX_STATES][BLC_MAX_STATES]; /** Output matrix C ∈ R^{p×n} */
    double   x[BLC_MAX_STATES];               /** Current state vector */
    double   x_hat[BLC_MAX_STATES];           /** Estimated/decoded state at controller */
    double   eigenvalues[BLC_MAX_STATES];      /** Eigenvalues of A (complex stored as re+im) */
    double   eigenvalues_im[BLC_MAX_STATES];   /** Imaginary parts of eigenvalues */
    int      n_unstable;                /** Number of unstable eigenvalues (Re>0) */
    double   min_datarate;             /** Minimum data rate for stabilization (bits/sample) */
} BLCPlant;

/** Packet structure for control communication.
 *  Each packet carries quantized state/control information
 *  across the bandwidth-limited channel.
 */
typedef struct {
    uint8_t  data[BLC_MAX_PACKET_BYTES]; /** Payload */
    int      length_bits;                /** Payload length in bits */
    int      seq_num;                    /** Sequence number */
    double   timestamp;                  /** Transmission time */
    bool     is_measurement;             /** true = sensor→controller, false = controller→actuator */
    bool     is_ack;                     /** ACK/NACK flag */
    double   checksum;                   /** Simple checksum for error detection */
} BLCPacket;

/** Bandwidth-limited controller state.
 *  Encapsulates the full state of a bandwidth-limited control loop.
 */
typedef struct {
    BLCChannel   channel;           /** Communication channel parameters */
    BLCQuantizer sensor_quant;      /** Sensor-side quantizer */
    BLCQuantizer ctrl_quant;        /** Control-side quantizer */
    BLCPlant     plant;             /** Plant model and state */
    double       K[BLC_MAX_STATES][BLC_MAX_STATES]; /** Controller gain matrix */
    double       L[BLC_MAX_STATES][BLC_MAX_STATES]; /** Observer gain matrix (if used) */
    double       sample_period;     /** Sampling period T_s (seconds) */
    int          bit_rate;          /** Current bit rate allocation */
    int          packets_sent;      /** Total packets transmitted */
    int          packets_lost;      /** Total packets lost */
    int          overloads;         /** Total quantization overloads */
    double       control_cost;      /** Accumulated LQR cost */
    bool         is_stable;         /** Stability flag */
    int          sim_steps;         /** Simulation step counter */
    double       sim_time;          /** Current simulation time */
} BLCSystem;

/* ================================================================
 * L2: Core Concepts — Bandwidth-limited stabilization
 * ================================================================
 *
 * The Data Rate Theorem (Nair & Evans, 2004):
 *   For a discrete-time LTI system with state matrix A,
 *   the minimum average data rate required for mean-square
 *   stabilizability is:
 *       R > Σ_{i: |λᵢ(A)| ≥ 1} log₂|λᵢ(A)|  bits/sample
 *
 * For continuous-time with sampling period T_s:
 *   R > Σ_{i: Re(λᵢ(A)) > 0} (2·Re(λᵢ(A))·T_s) / ln(2)  bits/sample
 *
 * This is the control-theoretic analog of Shannon's source coding theorem:
 * the unstable eigenvalues of the plant determine the information rate
 * necessary for stabilization.
 */

/** Convert dB to linear ratio */
#define BLC_DB_TO_LINEAR(db) (pow(10.0, (db) / 10.0))

/** Convert linear ratio to dB */
#define BLC_LINEAR_TO_DB(lin) (10.0 * log10(lin))

/** Shannon capacity: C = B * log2(1 + S/N) */
#define BLC_SHANNON_CAPACITY(B, SNR) ((B) * log2(1.0 + (SNR)))

/** Nyquist rate: minimum sample rate for avoiding aliasing */
#define BLC_NYQUIST_RATE(BW) (2.0 * (BW))

/* ================================================================
 * L3: Mathematical Structures
 * ================================================================ */

/**
 * Discrete-time plant evolution (Euler approximation):
 *   x[k+1] = (I + T_s·A)·x[k] + T_s·B·u[k] + w[k]
 *
 * With quantized state measurement z[k] = Q(C·x[k] + v[k]),
 * the controller receives only z[k] at rate R bits/sample.
 */

/* API: Core System Lifecycle */
BLCSystem* blc_create(int n_states, int n_inputs, int n_outputs);
void       blc_free(BLCSystem* sys);
int        blc_init_channel(BLCSystem* sys, double bw_hz, double snr,
                            double latency_ms, double loss_rate);
int        blc_init_plant(BLCSystem* sys, const double* A, const double* B,
                          const double* C);
int        blc_set_initial_state(BLCSystem* sys, const double* x0);
int        blc_set_controller_gain(BLCSystem* sys, const double* K);
int        blc_set_observer_gain(BLCSystem* sys, const double* L);

/* API: Quantizer Configuration */
int        blc_quantizer_init(BLCQuantizer* q, int levels,
                              double lo, double hi, bool logarithmic);
double     blc_quantize(BLCQuantizer* q, double value);
int        blc_quantize_to_index(BLCQuantizer* q, double value);
double     blc_dequantize(BLCQuantizer* q, int index);
double     blc_quantization_error(BLCQuantizer* q, double value);
int        blc_quantizer_resize(BLCQuantizer* q, double new_lo,
                                double new_hi);
double     blc_quantizer_max_error(BLCQuantizer* q, double range);

/* API: Information-Theoretic Analysis */
double     blc_channel_capacity(const BLCChannel* ch);
double     blc_minimum_datarate(const BLCPlant* plant);
double     blc_minimum_datarate_ct(const BLCPlant* plant, double Ts);
int        blc_count_unstable_eigenvalues(BLCPlant* plant);
double     blc_max_unstable_eigenvalue(const BLCPlant* plant);
double     blc_eigenvalue_sum_log(const BLCPlant* plant);
double     blc_entropy_rate(const BLCPlant* plant, double Ts);

/* API: Packet Operations */
BLCPacket* blc_packet_create(int length_bits, int seq, double ts,
                             bool is_measurement);
void       blc_packet_free(BLCPacket* pkt);
int        blc_packet_encode_state(const BLCSystem* sys, BLCPacket* pkt);
int        blc_packet_decode_state(BLCSystem* sys, const BLCPacket* pkt);
int        blc_packet_encode_control(const BLCSystem* sys, BLCPacket* pkt,
                                     const double* u);
int        blc_packet_decode_control(BLCSystem* sys, const BLCPacket* pkt,
                                     double* u);
bool       blc_packet_is_corrupt(const BLCPacket* pkt);

/* API: Simulation */
int        blc_simulate_step(BLCSystem* sys);
int        blc_simulate(BLCSystem* sys, double duration, double* t_out,
                        double* x_out, int max_steps);
double     blc_compute_lqr_cost(const BLCSystem* sys, const double* Q,
                                const double* R);
int        blc_check_stability(BLCSystem* sys);

/* API: Accessors */
void       blc_get_state(const BLCSystem* sys, double* x);
void       blc_get_estimated_state(const BLCSystem* sys, double* x_hat);
double     blc_get_control_effort(const BLCSystem* sys, double* u);
int        blc_get_packet_stats(const BLCSystem* sys, int* sent,
                                int* lost, int* overloads);
double     blc_get_simulation_time(const BLCSystem* sys);
void       blc_print_state(const BLCSystem* sys, FILE* stream);

#endif /* BLC_CORE_H */