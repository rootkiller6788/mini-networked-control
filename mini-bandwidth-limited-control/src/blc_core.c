/**
 * blc_core.c — Bandwidth-Limited Control Core Implementation
 *
 * Core system lifecycle, quantizer operations, packet I/O,
 * and simulation engine for bandwidth-limited control systems.
 *
 * Implements the fundamental concepts from:
 *   - Wong & Brockett (1997, 1999): finite bandwidth constraints
 *   - Nair & Evans (2004): stabilizability with finite data rates
 *   - Tatikonda & Mitter (2004): control under communication constraints
 *
 * Knowledge coverage: L1 Definitions, L2 Core Concepts, L3 Math Structures
 */

#include "blc_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * System Lifecycle
 * ================================================================ */

BLCSystem* blc_create(int n_states, int n_inputs, int n_outputs) {
    if (n_states <= 0 || n_states > BLC_MAX_STATES) return NULL;
    if (n_inputs <= 0 || n_inputs > BLC_MAX_STATES) return NULL;
    if (n_outputs <= 0 || n_outputs > BLC_MAX_STATES) return NULL;

    BLCSystem* sys = (BLCSystem*)calloc(1, sizeof(BLCSystem));
    if (!sys) return NULL;

    sys->plant.n_states  = n_states;
    sys->plant.n_inputs  = n_inputs;
    sys->plant.n_outputs = n_outputs;
    sys->sample_period   = 0.01;  /** Default 10ms sampling */
    sys->is_stable       = true;
    sys->sim_time        = 0.0;
    sys->sim_steps       = 0;

    /** Default quantization: 256 levels over [-10, 10] */
    blc_quantizer_init(&sys->sensor_quant, 256, -10.0, 10.0, false);
    blc_quantizer_init(&sys->ctrl_quant, 256, -10.0, 10.0, false);

    return sys;
}

void blc_free(BLCSystem* sys) {
    if (!sys) return;
    free(sys);
}

int blc_init_channel(BLCSystem* sys, double bw_hz, double snr,
                     double latency_ms, double loss_rate) {
    if (!sys) return -1;
    if (bw_hz <= 0.0 || snr <= 0.0) return -2;
    if (loss_rate < 0.0 || loss_rate > 1.0) return -3;

    sys->channel.bandwidth_hz     = bw_hz;
    sys->channel.snr              = snr;
    sys->channel.capacity_bps     = bw_hz * log2(1.0 + snr);
    sys->channel.latency_ms       = latency_ms;
    sys->channel.packet_loss_rate  = loss_rate;
    sys->channel.jitter_ms        = 0.0;
    return 0;
}

int blc_init_plant(BLCSystem* sys, const double* A, const double* B,
                   const double* C) {
    if (!sys || !A || !B || !C) return -1;
    int n = sys->plant.n_states;
    int m = sys->plant.n_inputs;
    int p = sys->plant.n_outputs;

    /** Copy A matrix (row-major, n×n) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sys->plant.A[i][j] = A[i * n + j];
        }
    }

    /** Copy B matrix (row-major, n×m) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            sys->plant.B[i][j] = B[i * m + j];
        }
    }

    /** Copy C matrix (row-major, p×n) */
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < n; j++) {
            sys->plant.C[i][j] = C[i * n + j];
        }
    }

    return 0;
}

int blc_set_initial_state(BLCSystem* sys, const double* x0) {
    if (!sys || !x0) return -1;
    int n = sys->plant.n_states;
    for (int i = 0; i < n; i++) {
        sys->plant.x[i] = x0[i];
        sys->plant.x_hat[i] = x0[i];  /** Assume perfect initial knowledge */
    }
    return 0;
}

int blc_set_controller_gain(BLCSystem* sys, const double* K) {
    if (!sys || !K) return -1;
    int n = sys->plant.n_inputs;
    int m = sys->plant.n_states;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            sys->K[i][j] = K[i * m + j];
        }
    }
    return 0;
}

int blc_set_observer_gain(BLCSystem* sys, const double* L) {
    if (!sys || !L) return -1;
    int n = sys->plant.n_states;
    int p = sys->plant.n_outputs;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            sys->L[i][j] = L[i * p + j];
        }
    }
    return 0;
}

/* ================================================================
 * Quantizer Operations
 * ================================================================ */

int blc_quantizer_init(BLCQuantizer* q, int levels, double lo, double hi,
                       bool logarithmic) {
    if (!q || levels < 2 || lo >= hi) return -1;
    if (levels > BLC_MAX_QUANT_LEVELS) return -2;

    q->levels         = levels;
    q->range_lo       = lo;
    q->range_hi       = hi;
    q->step           = (hi - lo) / (double)levels;
    q->max_error      = q->step * 0.5;
    q->overload_prob   = 0.0;
    q->overload_count  = 0;
    q->is_logarithmic  = logarithmic;
    q->log_base        = logarithmic ? 2.0 : 0.0;
    return 0;
}

double blc_quantize(BLCQuantizer* q, double value) {
    if (!q) return 0.0;

    /** Clamp to range */
    if (value < q->range_lo) {
        q->overload_count++;
        value = q->range_lo;
    } else if (value > q->range_hi) {
        q->overload_count++;
        value = q->range_hi;
    }

    if (q->is_logarithmic) {
        /** Logarithmic quantization: q_round(ln(value/lo) / ln(base) * levels) */
        if (value <= 0) value = q->range_lo;
        double normalized = log(value / q->range_lo) / log(q->log_base);
        int idx = (int)round(normalized);
        if (idx < 0) idx = 0;
        if (idx >= q->levels) idx = q->levels - 1;
        return q->range_lo * pow(q->log_base, (double)idx);
    } else {
        /** Uniform quantization: round to nearest step */
        int idx = (int)round((value - q->range_lo) / q->step);
        if (idx < 0) idx = 0;
        if (idx >= q->levels) idx = q->levels - 1;
        return q->range_lo + (double)idx * q->step;
    }
}

int blc_quantize_to_index(BLCQuantizer* q, double value) {
    if (!q) return -1;

    if (value < q->range_lo) {
        q->overload_count++;
        return 0;
    }
    if (value > q->range_hi) {
        q->overload_count++;
        return q->levels - 1;
    }

    if (q->is_logarithmic) {
        if (value <= 0) return 0;
        double normalized = log(value / q->range_lo) / log(q->log_base);
        int idx = (int)round(normalized);
        if (idx < 0) return 0;
        if (idx >= q->levels) return q->levels - 1;
        return idx;
    } else {
        int idx = (int)round((value - q->range_lo) / q->step);
        if (idx < 0) return 0;
        if (idx >= q->levels) return q->levels - 1;
        return idx;
    }
}

double blc_dequantize(BLCQuantizer* q, int index) {
    if (!q || index < 0 || index >= q->levels) return 0.0;

    if (q->is_logarithmic) {
        return q->range_lo * pow(q->log_base, (double)index);
    } else {
        return q->range_lo + (double)index * q->step;
    }
}

double blc_quantization_error(BLCQuantizer* q, double value) {
    if (!q) return 0.0;
    double q_val = blc_quantize(q, value);
    return value - q_val;
}

int blc_quantizer_resize(BLCQuantizer* q, double new_lo, double new_hi) {
    if (!q || new_lo >= new_hi) return -1;
    q->range_lo  = new_lo;
    q->range_hi  = new_hi;
    q->step      = (new_hi - new_lo) / (double)q->levels;
    q->max_error = q->step * 0.5;
    return 0;
}

double blc_quantizer_max_error(BLCQuantizer* q, double range) {
    if (!q) return 0.0;
    return range / (double)q->levels;
}

/* ================================================================
 * Information-Theoretic Analysis
 * ================================================================ */

double blc_channel_capacity(const BLCChannel* ch) {
    if (!ch) return 0.0;
    return ch->bandwidth_hz * log2(1.0 + ch->snr);
}

double blc_minimum_datarate(const BLCPlant* plant) {
    /** Discrete-time Data Rate Theorem:
     *  R_min = Σ log₂|λᵢ| for |λᵢ| > 1
     *  (bits per sample)
     */
    if (!plant) return 0.0;
    double rate = 0.0;
    for (int i = 0; i < plant->n_states; i++) {
        double mag = sqrt(plant->eigenvalues[i] * plant->eigenvalues[i] +
                          plant->eigenvalues_im[i] * plant->eigenvalues_im[i]);
        if (mag > 1.0) {
            rate += log2(mag);
        }
    }
    return rate;
}

double blc_minimum_datarate_ct(const BLCPlant* plant, double Ts) {
    /** Continuous-time Data Rate Theorem:
     *  R_min = Σ (2·Re(λᵢ)·T_s) / ln(2) for Re(λᵢ) > 0
     *  (bits per second, lower bound)
     *
     *  More precisely, for sampling period T_s:
     *  The discrete-time eigenvalues are e^{λᵢ·T_s}.
     *  |e^{λᵢ·T_s}| > 1 ⇔ Re(λᵢ) > 0.
     *  Σ log₂(e^{Re(λᵢ)·T_s}) = Σ Re(λᵢ)·T_s / ln(2)
     *
     *  The factor 2 accounts for complex-conjugate pairs.
     */
    if (!plant || Ts <= 0.0) return 0.0;
    double rate = 0.0;
    for (int i = 0; i < plant->n_states; i++) {
        if (plant->eigenvalues[i] > 0.0) {
            /** For real unstable eigenvalues: λ·T_s / ln(2) */
            rate += plant->eigenvalues[i] / log(2.0);
        }
    }
    /** Convert from bits/sample to bits/sec */
    return rate / Ts;
}

int blc_count_unstable_eigenvalues(BLCPlant* plant) {
    if (!plant) return 0;
    int count = 0;
    for (int i = 0; i < plant->n_states; i++) {
        if (plant->eigenvalues[i] > 1e-10) {
            count++;
        }
    }
    plant->n_unstable = count;
    return count;
}

double blc_max_unstable_eigenvalue(const BLCPlant* plant) {
    if (!plant) return 0.0;
    double max_ev = 0.0;
    for (int i = 0; i < plant->n_states; i++) {
        if (plant->eigenvalues[i] > max_ev) {
            max_ev = plant->eigenvalues[i];
        }
    }
    return max_ev;
}

double blc_eigenvalue_sum_log(const BLCPlant* plant) {
    /** Σ log₂|λᵢ| for |λᵢ| > 1 — the data rate sum */
    if (!plant) return 0.0;
    double sum_log = 0.0;
    for (int i = 0; i < plant->n_states; i++) {
        double mag = sqrt(plant->eigenvalues[i] * plant->eigenvalues[i] +
                          plant->eigenvalues_im[i] * plant->eigenvalues_im[i]);
        if (mag > 1.0) {
            sum_log += log2(mag);
        }
    }
    return sum_log;
}

double blc_entropy_rate(const BLCPlant* plant, double Ts) {
    /** Entropy rate of the plant's state process.
     *  H = Σ_{i: Re(λᵢ)>0} 2·Re(λᵢ)  (nats/sec)
     *  Converted to bits: H_bits = H / ln(2)
     *
     *  This is the minimum information that must be transmitted
     *  per unit time to track the state with bounded error.
     */
    if (!plant) return 0.0;
    double H = 0.0;
    for (int i = 0; i < plant->n_states; i++) {
        if (plant->eigenvalues[i] > 1e-10) {
            H += 2.0 * plant->eigenvalues[i];
        }
    }
    return H / log(2.0);  /** Convert nats to bits */
    (void)Ts;
}

/* ================================================================
 * Packet Operations
 * ================================================================ */

BLCPacket* blc_packet_create(int length_bits, int seq, double ts,
                             bool is_measurement) {
    BLCPacket* pkt = (BLCPacket*)calloc(1, sizeof(BLCPacket));
    if (!pkt) return NULL;

    pkt->length_bits    = length_bits;
    pkt->seq_num        = seq;
    pkt->timestamp      = ts;
    pkt->is_measurement  = is_measurement;
    pkt->is_ack         = false;
    pkt->checksum       = 0.0;
    return pkt;
}

void blc_packet_free(BLCPacket* pkt) {
    free(pkt);
}

int blc_packet_encode_state(const BLCSystem* sys, BLCPacket* pkt) {
    /** Encode plant state vector into packet payload.
     *  Each state component is quantized and packed.
     *  Returns number of bits used. */
    if (!sys || !pkt) return -1;
    int n = sys->plant.n_states;
    int bits_per_state = (int)ceil(log2((double)sys->sensor_quant.levels));

    memset(pkt->data, 0, BLC_MAX_PACKET_BYTES);
    int bit_pos = 0;
    for (int i = 0; i < n; i++) {
        int idx = blc_quantize_to_index(
            (BLCQuantizer*)&sys->sensor_quant, sys->plant.x[i]);
        /** Pack idx into bit stream (simple bit packing) */
        for (int b = 0; b < bits_per_state && bit_pos < BLC_MAX_PACKET_BYTES * 8; b++) {
            if (idx & (1 << b)) {
                int byte_idx = bit_pos / 8;
                int bit_idx  = bit_pos % 8;
                pkt->data[byte_idx] |= (uint8_t)(1 << bit_idx);
            }
            bit_pos++;
        }
    }
    pkt->length_bits = bit_pos;
    return bit_pos;
}

int blc_packet_decode_state(BLCSystem* sys, const BLCPacket* pkt) {
    /** Decode packet payload back to estimated state.
     *  Returns number of bits consumed. */
    if (!sys || !pkt) return -1;
    int n = sys->plant.n_states;
    int bits_per_state = (int)ceil(log2((double)sys->sensor_quant.levels));
    int total_bits = bits_per_state * n;

    if (pkt->length_bits < total_bits) return -2;

    int bit_pos = 0;
    for (int i = 0; i < n; i++) {
        int idx = 0;
        for (int b = 0; b < bits_per_state; b++) {
            if (bit_pos >= pkt->length_bits) break;
            int byte_idx = bit_pos / 8;
            int bit_idx  = bit_pos % 8;
            if (pkt->data[byte_idx] & (1 << bit_idx)) {
                idx |= (1 << b);
            }
            bit_pos++;
        }
        sys->plant.x_hat[i] = blc_dequantize(
            (BLCQuantizer*)&sys->sensor_quant, idx);
    }
    return bit_pos;
}

int blc_packet_encode_control(const BLCSystem* sys, BLCPacket* pkt,
                              const double* u) {
    if (!sys || !pkt || !u) return -1;
    int m = sys->plant.n_inputs;
    int bits_per_ctrl = (int)ceil(log2((double)sys->ctrl_quant.levels));

    memset(pkt->data, 0, BLC_MAX_PACKET_BYTES);
    int bit_pos = 0;
    for (int i = 0; i < m; i++) {
        int idx = blc_quantize_to_index(
            (BLCQuantizer*)&sys->ctrl_quant, u[i]);
        for (int b = 0; b < bits_per_ctrl && bit_pos < BLC_MAX_PACKET_BYTES * 8; b++) {
            if (idx & (1 << b)) {
                int byte_idx = bit_pos / 8;
                int bit_idx  = bit_pos % 8;
                pkt->data[byte_idx] |= (uint8_t)(1 << bit_idx);
            }
            bit_pos++;
        }
    }
    pkt->length_bits = bit_pos;
    return bit_pos;
}

int blc_packet_decode_control(BLCSystem* sys, const BLCPacket* pkt,
                              double* u) {
    if (!sys || !pkt || !u) return -1;
    int m = sys->plant.n_inputs;
    int bits_per_ctrl = (int)ceil(log2((double)sys->ctrl_quant.levels));

    int bit_pos = 0;
    for (int i = 0; i < m; i++) {
        int idx = 0;
        for (int b = 0; b < bits_per_ctrl; b++) {
            if (bit_pos >= pkt->length_bits) break;
            int byte_idx = bit_pos / 8;
            int bit_idx  = bit_pos % 8;
            if (pkt->data[byte_idx] & (1 << bit_idx)) {
                idx |= (1 << b);
            }
            bit_pos++;
        }
        u[i] = blc_dequantize((BLCQuantizer*)&sys->ctrl_quant, idx);
    }
    return bit_pos;
}

bool blc_packet_is_corrupt(const BLCPacket* pkt) {
    /** Simple checksum verification: sum of all bytes.
     *  In a real system, this would be a CRC. */
    if (!pkt) return true;
    double sum = 0.0;
    int n_bytes = (pkt->length_bits + 7) / 8;
    for (int i = 0; i < n_bytes && i < BLC_MAX_PACKET_BYTES; i++) {
        sum += (double)pkt->data[i];
    }
    return fabs(sum - pkt->checksum) > 1e-6;
}

/* ================================================================
 * Simulation Engine
 * ================================================================ */

int blc_simulate_step(BLCSystem* sys) {
    if (!sys) return -1;
    int n = sys->plant.n_states;
    int m = sys->plant.n_inputs;
    double dt = sys->sample_period;

    /** Compute control: u = -K * x_hat (quantized state feedback) */
    double u[BLC_MAX_STATES] = {0};
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            u[i] -= sys->K[i][j] * sys->plant.x_hat[j];
        }
    }

    /** Euler integration: x[k+1] = x[k] + dt * (A*x[k] + B*u[k]) */
    double dx[BLC_MAX_STATES] = {0};
    for (int i = 0; i < n; i++) {
        double Ax = 0.0, Bu = 0.0;
        for (int j = 0; j < n; j++) {
            Ax += sys->plant.A[i][j] * sys->plant.x[j];
        }
        for (int j = 0; j < m; j++) {
            Bu += sys->plant.B[i][j] * u[j];
        }
        dx[i] = Ax + Bu;
    }
    for (int i = 0; i < n; i++) {
        sys->plant.x[i] += dt * dx[i];
    }

    /** Update estimator (open-loop prediction for now) */
    for (int i = 0; i < n; i++) {
        double Ax_hat = 0.0, Bu = 0.0;
        for (int j = 0; j < n; j++) {
            Ax_hat += sys->plant.A[i][j] * sys->plant.x_hat[j];
        }
        for (int j = 0; j < m; j++) {
            Bu += sys->plant.B[i][j] * u[j];
        }
        sys->plant.x_hat[i] += dt * (Ax_hat + Bu);
    }

    sys->sim_time += dt;
    sys->sim_steps++;
    sys->packets_sent++;

    return 0;
}

int blc_simulate(BLCSystem* sys, double duration, double* t_out,
                 double* x_out, int max_steps) {
    if (!sys || duration <= 0.0) return -1;
    int n = sys->plant.n_states;
    int steps = (int)(duration / sys->sample_period);
    if (max_steps > 0 && steps > max_steps) steps = max_steps;

    for (int k = 0; k < steps; k++) {
        blc_simulate_step(sys);
        if (t_out) t_out[k] = sys->sim_time;
        if (x_out) {
            for (int i = 0; i < n; i++) {
                x_out[k * n + i] = sys->plant.x[i];
            }
        }
    }
    return steps;
}

double blc_compute_lqr_cost(const BLCSystem* sys, const double* Q,
                            const double* R) {
    if (!sys) return 0.0;
    int n = sys->plant.n_states;
    int m = sys->plant.n_inputs;
    double cost = 0.0;

    /** State cost: x'Qx */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cost += sys->plant.x[i] * Q[i * n + j] * sys->plant.x[j];
        }
    }

    /** Input cost: u'Ru — compute u from x_hat and K */
    double u[BLC_MAX_STATES] = {0};
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            u[i] -= sys->K[i][j] * sys->plant.x_hat[j];
        }
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            cost += u[i] * R[i * m + j] * u[j];
        }
    }

    return cost;
}

int blc_check_stability(BLCSystem* sys) {
    if (!sys) return -1;
    int n = sys->plant.n_states;
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += sys->plant.x[i] * sys->plant.x[i];
    }
    norm = sqrt(norm);

    /** Simple check: state norm bounded */
    if (norm > 1e6) {
        sys->is_stable = false;
        return 0;
    }
    sys->is_stable = true;
    return 1;
}

/* ================================================================
 * Accessors
 * ================================================================ */

void blc_get_state(const BLCSystem* sys, double* x) {
    if (!sys || !x) return;
    int n = sys->plant.n_states;
    for (int i = 0; i < n; i++) {
        x[i] = sys->plant.x[i];
    }
}

void blc_get_estimated_state(const BLCSystem* sys, double* x_hat) {
    if (!sys || !x_hat) return;
    int n = sys->plant.n_states;
    for (int i = 0; i < n; i++) {
        x_hat[i] = sys->plant.x_hat[i];
    }
}

double blc_get_control_effort(const BLCSystem* sys, double* u) {
    if (!sys) return 0.0;
    int n = sys->plant.n_states;
    int m = sys->plant.n_inputs;
    double effort = 0.0;
    for (int i = 0; i < m; i++) {
        u[i] = 0.0;
        for (int j = 0; j < n; j++) {
            u[i] -= sys->K[i][j] * sys->plant.x_hat[j];
        }
        effort += u[i] * u[i];
    }
    return sqrt(effort);
}

int blc_get_packet_stats(const BLCSystem* sys, int* sent,
                         int* lost, int* overloads) {
    if (!sys) return -1;
    if (sent)      *sent      = sys->packets_sent;
    if (lost)      *lost      = sys->packets_lost;
    if (overloads) *overloads = sys->overloads;
    return 0;
}

double blc_get_simulation_time(const BLCSystem* sys) {
    if (!sys) return -1.0;
    return sys->sim_time;
}

void blc_print_state(const BLCSystem* sys, FILE* stream) {
    if (!sys) return;
    FILE* out = stream ? stream : stdout;
    int n = sys->plant.n_states;
    fprintf(out, "BLCSystem t=%.4f s, stable=%s, steps=%d\n",
            sys->sim_time, sys->is_stable ? "yes" : "NO", sys->sim_steps);
    fprintf(out, "  State x:  ");
    for (int i = 0; i < n; i++) fprintf(out, "%+.4f ", sys->plant.x[i]);
    fprintf(out, "\n  Estimate: ");
    for (int i = 0; i < n; i++) fprintf(out, "%+.4f ", sys->plant.x_hat[i]);
    fprintf(out, "\n  Packets: sent=%d lost=%d overloads=%d\n",
            sys->packets_sent, sys->packets_lost, sys->overloads);
    fprintf(out, "  Channel: BW=%.2f Hz SNR=%.1f dB Cap=%.2f bps\n",
            sys->channel.bandwidth_hz,
            10.0 * log10(sys->channel.snr),
            sys->channel.capacity_bps);
}