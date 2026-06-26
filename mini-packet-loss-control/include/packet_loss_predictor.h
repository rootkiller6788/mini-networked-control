#ifndef PACKET_LOSS_PREDICTOR_H
#define PACKET_LOSS_PREDICTOR_H

#include "packet_loss_core.h"
#include <stdbool.h>

/* ============================================================================
 * Packet Loss Predictor — Predictive & Proactive Compensation
 *
 * Implements predictive strategies to compensate for packet loss:
 *   1. Packetized predictive control (PPC) — send a sequence of future
 *      controls in each packet; the actuator uses the most recent
 *      prediction if current packet is lost.
 *   2. Model-based prediction — predict missing sensor data using
 *      the process model.
 *   3. Hold strategies with optimal prediction.
 *
 * Key references:
 *   - Bemporad (1998): "Predictive Control of Teleoperated Constrained
 *     Systems with Bounded Time-Varying Delays"
 *   - Quevedo & Nesic (2012): "Robust Stability of Packetized Predictive
 *     Control of Nonlinear Systems with Disturbances and Markovian
 *     Packet Losses"
 *   - Tang & de Silva (2006): "Compensation for Transmission Delays
 *     in an Ethernet-based Control Network Using Variable-Horizon
 *     Predictive Control"
 *   - Quevedo, Jurado, Silva (2012): "Packetized Predictive Control
 *     over Erasure Channels"
 * ============================================================================ */

/* --- Packetized Predictive Control (PPC) --- */

/**
 * Packetized Predictive Control: Instead of sending a single control
 * signal u_k, the controller sends a buffer of predicted future controls:
 *   [u_{k|k}, u_{k+1|k}, ..., u_{k+H-1|k}]
 *
 * If a packet is lost, the actuator uses the next prediction
 * from the most recently received packet (i.e., it shifts the buffer).
 *
 * Horizon H: number of future steps included in each packet.
 * Larger H → more resilience against consecutive losses.
 *
 * Buffer consumption:
 *   At time k, if packet arrives → use u_{k|k} from new buffer
 *   If packet is lost → use u_{k|k-1} from previous buffer (shift)
 *
 * This guarantees stability as long as the number of consecutive
 * losses is ≤ H (assuming the system is open-loop stable over H steps).
 *
 * For unstable systems, PPC provides graceful degradation.
 */
typedef struct {
    int horizon;                    /* Number of future controls per packet */
    double** buffer;                /* 2D buffer: [horizon][m] */
    double* current_input;          /* Currently applied input, size m */
    int buffer_index;               /* How far into the buffer we've consumed */
    int consecutive_losses;         /* Current consecutive loss count */
    int max_packet_loss_tolerated;  /* = horizon (buffer exhaustion protection) */
    bool buffer_exhausted;          /* All predictions have been consumed */
    double* emergency_input;        /* Safe input when buffer exhausted */

    /* Plant model for prediction */
    const double* A;                /* n×n system matrix */
    const double* B;                /* n×m input matrix */
    const double* L;                /* LQR gain (m×n) for generating predictions */
    int n;                          /* State dimension */
    int m;                          /* Input dimension */

    /* Current state estimate for prediction */
    double* x_current;              /* n-vector */

    /* Statistics */
    unsigned long total_packets_sent;
    unsigned long total_packets_used_from_buffer;
    int max_consecutive_losses_seen;
    double avg_buffer_utilization;
} PacketizedPredictiveControl;

/* --- Model-Based Sensor Prediction --- */

/**
 * Model-based missing sensor predictor.
 *
 * When a sensor packet is lost, the receiver (controller) cannot directly
 * update its estimate. However, it can predict the missing measurement
 * using the system model:
 *
 *   ŷ_k = C · x̂_{k|k-1}  (open-loop prediction of measurement)
 *
 * More sophisticated approaches use the innovation sequence
 * statistics to detect when the model is drifting.
 *
 * The prediction quality degrades over consecutive losses since
 * the open-loop prediction error grows: ||x̂_pred - x_true|| ≈ ρ(A)^N
 * where N is the number of consecutive losses.
 */
typedef struct {
    /* System model */
    const double* A;                /* n×n */
    const double* C;                /* p×n */
    int n;                          /* State dimension */
    int p;                          /* Measurement dimension */

    /* Predictor state */
    double* x_pred;                 /* Predicted state, size n */
    double* y_pred;                 /* Predicted measurement, size p */
    double* P_pred;                 /* Predicted covariance, n×n */

    /* Prediction quality tracking */
    int consecutive_losses;
    double prediction_error;        /* ||y_true - y_pred|| when measurement arrives */
    double accumulated_drift;       /* Cumulative prediction error */

    /* Adaptive correction */
    double correction_gain;         /* How much to trust prediction correction */
    double drift_threshold;         /* When to flag model mismatch */
    bool model_valid;               /* Is the model still adequate? */
} ModelSensorPredictor;

/* --- Optimal Hold Strategy Selector --- */

/**
 * Adaptive hold strategy selector.
 *
 * Chooses the best hold strategy based on:
 *   - Current loss pattern (isolated vs. burst)
 *   - System dynamics (stable vs. unstable open-loop)
 *   - Available computational resources
 *
 * Strategy selection logic:
 *   - Open-loop stable + isolated loss → ZERO_ORDER (simple and effective)
 *   - Open-loop stable + burst loss → PREDICTIVE (compensate for burst)
 *   - Open-loop unstable + isolated loss → PREDICTIVE (need model)
 *   - Open-loop unstable + burst loss → LQG_OPTIMAL (best possible)
 *   - Very high loss → ZERO_INPUT (safety)
 */
typedef struct {
    HoldStrategy current_strategy;
    double zero_order_cost;         /* Estimated cost of ZOH */
    double predictive_cost;         /* Estimated cost of predictive */
    double lqg_optimal_cost;        /* Estimated cost of LQG-optimal */
    double safety_cost;             /* Estimated cost of zero-input */

    /* Loss pattern characterization */
    double recent_loss_rate;        /* Loss rate in recent window */
    double burst_index;             /* Burstiness measure */
    int consecutive_losses;

    /* System characterization */
    bool open_loop_stable;
    double spectral_radius;

    /* Performance history per strategy */
    double cost_zoh_history[100];
    double cost_pred_history[100];
    double cost_lqg_history[100];
    int history_idx;
} HoldStrategySelector;

/* --- Network-Aware Reference Governor --- */

/**
 * Reference governor for NCS with packet loss.
 *
 * Modifies the reference signal to ensure constraint satisfaction
 * despite packet loss. The governor predicts future system evolution
 * under worst-case loss scenarios and adjusts the reference accordingly.
 *
 * Approach (Bemporad, 1998):
 *   At time k, compute the maximum reachable set over the next
 *   H_max_loss steps assuming all packets are lost.
 *   If this set violates constraints, scale back the reference.
 *
 * This prevents constraint violations caused by packet loss induced
 * overshoot or instability.
 */
typedef struct {
    double* reference;              /* Original (desired) reference */
    double* modified_reference;     /* Governor-adjusted reference */
    int ref_size;

    /* Constraint set */
    double* y_min;                  /* Lower bound on output */
    double* y_max;                  /* Upper bound on output */
    int n_constraints;

    /* Predictor for constraint enforcement */
    const double* A;                /* System matrix */
    const double* B;                /* Input matrix */
    const double* C;                /* Output matrix */
    int n;
    int m;
    int p;

    int max_consecutive_losses;     /* Worst-case consecutive losses to protect */
    double safety_factor;           /* Scaling factor ∈ (0,1] */
    double* x_pred;                 /* Predicted state buffer */
} ReferenceGovernor;

/* --- Loss-Constrained Optimal Controller --- */

/**
 * Controller that explicitly optimizes over loss scenarios.
 *
 * Formulation (Gupta et al., 2007 for UDP-like):
 *   min_u E[ Σ (x_k'Qx_k + u_k'Ru_k) | info at k ]
 *
 * The expectation is over the random loss process. The optimal control
 * is a function of the belief state (distribution over plant state).
 *
 * For finite-horizon with Bernoulli loss:
 *   The cost-to-go is piecewise quadratic in the state estimate.
 *   The optimal policy is linear in the state estimate but the
 *   gain depends on the time-to-go and loss probability.
 */
typedef struct {
    const double* A;
    const double* B;
    const double* Q;                /* State cost, n×n */
    const double* R;                /* Control cost, m×m */
    int n;
    int m;

    /* Time-varying LQR gains for finite horizon */
    double** gain_sequence;         /* K[t]: m×n gain at time t */
    double** cost_to_go;            /* P[t]: n×n cost-to-go matrix at time t */
    int horizon;

    /* Loss model */
    double loss_probability;        /* Bernoulli loss rate, p ∈ [0,1] */
    TransportProtocol protocol;     /* TCP vs UDP-like */

    /* Current control */
    double* u_optimal;              /* Optimal control for current step */
} LossConstrainedController;

/* ============================================================================
 * Packetized Predictive Control API
 * ============================================================================ */

/**
 * Create PPC controller with given prediction horizon.
 * @param horizon: Number of future controls to pack per packet
 * @param A, B: System matrices (references)
 * @param L: LQR gain matrix (m×n, row-major)
 * @param x0: Initial state estimate
 */
PacketizedPredictiveControl* pl_ppc_create(int horizon,
                                             const double* A, const double* B,
                                             const double* L,
                                             int n, int m,
                                             const double* x0);
void pl_ppc_free(PacketizedPredictiveControl* ppc);

/**
 * Generate a new packet buffer of predicted future controls.
 * Fills buffer[t] = -L · (A+B·L)^t · x_current for t = 0, ..., H-1.
 *
 * The generated sequence is open-loop optimal if no further measurements
 * arrive (certainty-equivalent prediction).
 */
void pl_ppc_generate_buffer(PacketizedPredictiveControl* ppc);

/**
 * Consume the next control input from the buffer.
 *
 * If new_packet_arrived: reset buffer index, use u_{0|k}
 * If lost: advance buffer index, use u_{buffer_index|k-H_last}
 * If buffer exhausted: apply emergency input
 *
 * Returns pointer to the applied control signal (size m).
 */
const double* pl_ppc_consume(PacketizedPredictiveControl* ppc,
                              bool new_packet_arrived);

/** Update the state estimate used for generating predictions. */
void pl_ppc_update_state(PacketizedPredictiveControl* ppc, const double* x_new);

void pl_ppc_reset(PacketizedPredictiveControl* ppc, const double* x0);
void pl_ppc_print(const PacketizedPredictiveControl* ppc);

/* ============================================================================
 * Model-Based Sensor Predictor API
 * ============================================================================ */

ModelSensorPredictor* pl_msp_create(const double* A, const double* C,
                                      int n, int p,
                                      const double* x0, const double* P0);
void pl_msp_free(ModelSensorPredictor* msp);

/**
 * Predict missing measurement: ŷ = C · Â·x_pred (open-loop prediction).
 * Updates internal state and returns pointer to predicted measurement.
 */
const double* pl_msp_predict_missing(ModelSensorPredictor* msp);

/**
 * Correct the predictor when a measurement actually arrives.
 * Uses the innovation y_true - y_pred to correct the model.
 */
void pl_msp_correct(ModelSensorPredictor* msp, const double* y_true);

/**
 * Combined step: if arrived → correct; else → predict_missing.
 * Returns the measurement (true or predicted).
 */
const double* pl_msp_step(ModelSensorPredictor* msp,
                           const double* y, bool arrived);

double pl_msp_get_prediction_error(const ModelSensorPredictor* msp);
bool pl_msp_is_model_valid(const ModelSensorPredictor* msp);
void pl_msp_print(const ModelSensorPredictor* msp);

/* ============================================================================
 * Hold Strategy Selector API
 * ============================================================================ */

HoldStrategySelector* pl_hss_create(bool open_loop_stable, double spectral_radius);
void pl_hss_free(HoldStrategySelector* hss);

/**
 * Select the best hold strategy based on current conditions.
 * Updates internal cost estimates and returns the recommended strategy.
 */
HoldStrategy pl_hss_select(HoldStrategySelector* hss,
                            double loss_rate, double burst_index,
                            int consecutive_losses,
                            double* estimated_cost_out);

void pl_hss_record_performance(HoldStrategySelector* hss,
                                HoldStrategy strategy, double cost);
void pl_hss_print(const HoldStrategySelector* hss);

/* ============================================================================
 * Reference Governor API
 * ============================================================================ */

ReferenceGovernor* pl_rg_create(const double* A, const double* B, const double* C,
                                  int n, int m, int p,
                                  const double* reference, int ref_size,
                                  const double* y_min, const double* y_max,
                                  int max_loss);
void pl_rg_free(ReferenceGovernor* rg);

/**
 * Compute the governor-adjusted reference for the next step.
 * Ensures that even under worst-case (all packets lost for max_loss steps),
 * the output stays within [y_min, y_max].
 *
 * Returns pointer to modified reference.
 */
const double* pl_rg_adjust(ReferenceGovernor* rg, const double* x_current);

void pl_rg_print(const ReferenceGovernor* rg);

/* ============================================================================
 * Loss-Constrained Controller API
 * ============================================================================ */

LossConstrainedController* pl_lcc_create(const double* A, const double* B,
                                           const double* Q, const double* R,
                                           int n, int m, int horizon,
                                           double loss_prob,
                                           TransportProtocol protocol);
void pl_lcc_free(LossConstrainedController* lcc);

/**
 * Compute the loss-constrained optimal control for the current state.
 *
 * For TCP-like: the optimal control is certainty-equivalent LQR.
 * For UDP-like: the gains are modified to account for delivery uncertainty.
 *
 * Reference: Gupta et al. (2007), Theorem 1 (UDP-like structure).
 */
const double* pl_lcc_compute(LossConstrainedController* lcc,
                              const double* x_estimate);

void pl_lcc_print(const LossConstrainedController* lcc);

#endif /* PACKET_LOSS_PREDICTOR_H */