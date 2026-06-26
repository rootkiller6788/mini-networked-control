/**
 * cloud_control_delay.h - Network Delay Models & Compensation
 *
 * Delay is the defining challenge of cloud control. This module provides
 * delay models (constant, time-varying, stochastic), stability analysis
 * under delay, and compensation strategies including the Smith predictor,
 * model predictive control with delay, and Kalman filtering with
 * out-of-sequence measurements.
 *
 * Domain: Networked Control / Time-Delay Systems
 * References:
 *   - Zhang et al., "Stability of Networked Control Systems", IEEE CSM, 2001
 *   - Smith, "A controller to overcome dead time", ISA Journal, 1959
 *   - Schenato, "Optimal estimation in networked control with packet losses",
 *     Automatica, 2008
 *   - Richard, "Time-delay systems: an overview", Automatica, 2003
 *
 * Knowledge Coverage:
 *   L1: NetworkDelay, DelayModel, SmithPredictor, DelayCompensator
 *   L2: Delay stability, MATI (Maximum Allowable Transfer Interval)
 *   L3: Lyapunov-Krasovskii functional, LMI conditions
 *   L4: Smith predictor principle, delay-dependent stability
 *   L5: Delay compensation algorithms, out-of-order measurement handling
 */

#ifndef CLOUD_CONTROL_DELAY_H
#define CLOUD_CONTROL_DELAY_H

#include "cloud_control_core.h"

/* ============================================================================
 * L1: Delay Model Types & Constants
 * ============================================================================ */

#define CCS_MAX_DELAY_HISTORY   1024
#define CCS_MAX_SMITH_STATES    64
#define CCS_DEFAULT_DELAY_US    5000.0
#define CCS_MAX_DELAY_US        100000.0

/** Statistical model for network delay distribution */
typedef enum {
    DELAY_MODEL_CONSTANT  = 0,  /* Fixed deterministic delay */
    DELAY_MODEL_UNIFORM   = 1,  /* Uniform distribution [min, max] */
    DELAY_MODEL_NORMAL    = 2,  /* Truncated normal distribution */
    DELAY_MODEL_EXPONENTIAL = 3, /* Exponential (memoryless) */
    DELAY_MODEL_GAMMA     = 4,  /* Gamma distribution (sum of exponentials) */
    DELAY_MODEL_MEASURED  = 5,  /* Empirically measured trace */
    DELAY_MODEL_PERIODIC  = 6,  /* Sinusoidal or periodic pattern */
    DELAY_MODEL_MARKOV    = 7   /* Markov-modulated (bursty) */
} DelayModelType;

/**
 * DelayStatistics — statistical characterization of network delay
 *
 * Collected from real network measurements or synthesized from models.
 * Used to parameterize stability analysis and compensation.
 */
typedef struct {
    double   mean_us;
    double   variance_us;
    double   min_us;
    double   max_us;
    double   median_us;
    double   p95_us;          /* 95th percentile */
    double   p99_us;          /* 99th percentile */
    double   jitter_us;       /* Mean absolute deviation */
    double   autocorrelation; /* Lag-1 autocorrelation */
    int      sample_count;
    int      outliers_count;  /* Samples beyond 3 sigma */
} DelayStatistics;

/**
 * NetworkDelay — single delay measurement or model configuration
 */
typedef struct {
    DelayModelType model_type;
    double         value_us;       /* Current delay value */
    double         param1;         /* Model parameter 1 (e.g., min for uniform) */
    double         param2;         /* Model parameter 2 (e.g., max for uniform) */
    double         timestamp;      /* Measurement timestamp */
    int            source_id;      /* Which network path */
    int            is_measured;    /* 1 = real measurement, 0 = synthetic */
} NetworkDelay;

/* ============================================================================
 * L1 / L3: Smith Predictor Structure
 * ============================================================================ */

/**
 * SmithPredictor — delay-compensating controller
 *
 * The Smith predictor (Smith, 1959) is the fundamental delay compensation
 * technique for cloud control. It uses an internal plant model to predict
 * the delay-free output, enabling the controller to act on predicted rather
 * than delayed measurements.
 *
 * Structure:
 *   G(s) = plant, G0(s) = delay-free plant model
 *   C(s) = primary controller
 *   tau  = estimated delay
 *
 * The predictor subtracts the delayed model output from the actual delayed
 * measurement and adds the delay-free model output:
 *   y_compensated = y_model(t) + [y_measured(t) - y_model(t - tau)]
 *
 * Mathematical structure (L3):
 *   Plant:           x_dot = A x + B u,   y = C x
 *   Delay-free model: x_m_dot = A_m x_m + B_m u,  y_m = C_m x_m
 *   Delayed model:    x_d_dot = A_m x_d + B_m u(t-tau),  y_d = C_m x_d
 */
typedef struct {
    CloudControlSystem *ccs;   /* Reference to parent system */

    /* Delay-free plant model (may differ from actual plant) */
    double A_m[CCS_MAX_STATES][CCS_MAX_STATES];
    double B_m[CCS_MAX_STATES][CCS_MAX_INPUTS];
    double C_m[CCS_MAX_OUTPUTS][CCS_MAX_STATES];
    int    n_m, m_m, p_m;

    /* Model states */
    double x_model[CCS_MAX_STATES];        /* Delay-free model state */
    double x_delayed[CCS_MAX_STATES];      /* Delayed model state */
    double u_history[CCS_MAX_DELAY_HISTORY][CCS_MAX_INPUTS];
    double u_time[CCS_MAX_DELAY_HISTORY];
    int    u_history_count;

    /* Delay estimate */
    double estimated_delay_us;
    double delay_adaptation_rate;  /* For adaptive delay estimation */

    /* Performance */
    double prediction_error;       /* y_actual - y_predicted */
    double compensation_ratio;     /* Improvement ratio over no compensation */

    /* Buffer for out-of-order measurements */
    double meas_buffer[CCS_MAX_HISTORY][CCS_MAX_OUTPUTS];
    double meas_ts[CCS_MAX_HISTORY];
    int    meas_count;
} SmithPredictor;

/* ============================================================================
 * L1: Delay Compensation Configuration
 * ============================================================================ */

/**
 * DelayCompensatorConfig — configuration for delay compensation strategy
 */
typedef enum {
    COMPENSATE_NONE          = 0,  /* No compensation */
    COMPENSATE_SMITH         = 1,  /* Smith predictor */
    COMPENSATE_MPC_DELAY     = 2,  /* MPC with explicit delay model */
    COMPENSATE_KALMAN_DELAY  = 3,  /* Kalman filter with delayed measurements */
    COMPENSATE_ADAPTIVE      = 4   /* Adaptive compensation */
} CompensationStrategy;

/* ============================================================================
 * L2: Delay Analysis API
 * ============================================================================ */

/**
 * delay_stats_compute — compute statistics from a delay trace
 *
 * @param delays     Array of delay measurements (microseconds)
 * @param count      Number of measurements
 * @param stats_out  Output statistics structure
 * @return           0 on success, -1 if count == 0
 */
int delay_stats_compute(const double *delays, int count,
                         DelayStatistics *stats_out);

/**
 * delay_stats_print — print delay statistics in human-readable format
 */
void delay_stats_print(const DelayStatistics *stats);

/**
 * delay_is_stable — determine if system is stable under given delay
 *
 * Uses eigenvalue analysis of the delayed closed-loop system.
 * For a system x_dot = (A - BK)x(t - tau), stability requires
 * all roots of det(sI - (A - BK)e^{-s*tau}) = 0 have negative real parts.
 *
 * Implements the Walton-Marshall (1987) method for delay-independent
 * stability check and the Zhang (2001) eigenvalue method for delay-dependent
 * stability.
 *
 * Complexity: O(n^3) for eigenvalue computation
 *
 * @param ccs       Cloud control system
 * @param delay_us  Network delay in microseconds
 * @return          1 if stable, 0 if unstable, -1 on error
 */
int delay_is_stable(const CloudControlSystem *ccs, double delay_us);

/**
 * delay_mati_compute — compute Maximum Allowable Transfer Interval
 *
 * The MATI is the maximum time between successive control updates
 * for which the system remains stable. Based on Walsh et al. (2002)
 * and Nesic & Teel (2004) for nonlinear NCS.
 *
 * For linear systems with constant delay tau:
 *   MATI satisfies: there exists P > 0 such that
 *   (A - BK)^T P + P (A - BK) + tau * ... < 0 (LMI condition)
 *
 * @param ccs  Cloud control system
 * @return     MATI in microseconds, or -1 if unbounded
 */
double delay_mati_compute(const CloudControlSystem *ccs);

/**
 * delay_lyapunov_krasovskii — verify stability via LK functional
 *
 * Lyapunov-Krasovskii functional for time-delay systems:
 *   V(x_t) = x^T(t) P x(t) + integral_{t-tau}^{t} x^T(s) Q x(s) ds
 *          + integral_{-tau}^{0} integral_{t+theta}^{t} x_dot^T(s) R x_dot(s) ds dtheta
 *
 * Stability condition: V_dot < 0 for all non-zero trajectories.
 *
 * This function checks the LMI conditions derived from V_dot < 0.
 *
 * @param ccs       Cloud control system
 * @param delay_us  Network delay
 * @return          1 if LK-stable, 0 if not, -1 on error
 */
int delay_lyapunov_krasovskii(const CloudControlSystem *ccs, double delay_us);

/* ============================================================================
 * L5: Delay Compensation API
 * ============================================================================ */

/**
 * smith_create — allocate and initialize a Smith predictor
 *
 * Creates a Smith predictor for the given cloud control system.
 * The internal model is initialized from the plant model in ccs.
 *
 * @param ccs  Parent cloud control system
 * @return     Initialized SmithPredictor, or NULL on failure
 */
SmithPredictor* smith_create(CloudControlSystem *ccs);

/**
 * smith_free — release resources associated with a Smith predictor
 */
void smith_free(SmithPredictor *sp);

/**
 * smith_set_delay — set the estimated delay for the Smith predictor
 *
 * @param sp        Smith predictor
 * @param delay_us  Estimated round-trip delay in microseconds
 */
void smith_set_delay(SmithPredictor *sp, double delay_us);

/**
 * smith_predict — compute delay-compensated output
 *
 * The core Smith predictor algorithm:
 *   1. Advance delay-free model: x_m_dot = A_m x_m + B_m u
 *   2. Advance delayed model:    x_d_dot = A_m x_d + B_m u_hist(t-tau)
 *   3. Compute compensation:     y_comp = y_meas - C_m x_d + C_m x_m
 *
 * This effectively "removes" the delay from the measurement, allowing
 * the controller to act on predicted current plant state.
 *
 * Complexity: O(n^2 + nm + pn) per call
 *
 * @param sp       Smith predictor
 * @param y_meas   Delayed measurement vector (dimension p)
 * @param u        Current control input (dimension m)
 * @param dt       Time step
 * @param y_comp   Output: delay-compensated measurement (pre-allocated)
 * @return         0 on success
 */
int smith_predict(SmithPredictor *sp, const double *y_meas,
                   const double *u, double dt, double *y_comp);

/**
 * smith_adapt_delay — adapt delay estimate based on prediction error
 *
 * Uses gradient descent on prediction error to refine the delay estimate.
 *   tau_new = tau_old + adaptation_rate * prediction_error * sensitivity
 *
 * @param sp  Smith predictor
 * @return    Updated delay estimate in microseconds
 */
double smith_adapt_delay(SmithPredictor *sp);

/**
 * smith_get_prediction_error — return current prediction error magnitude
 */
double smith_get_prediction_error(const SmithPredictor *sp);

/**
 * smith_handle_oos_measurement — handle out-of-sequence measurement
 *
 * When a delayed measurement arrives out of order, the Smith predictor
 * must re-process historical data to maintain consistency. This function
 * inserts the OOS measurement at the correct position in the history buffer
 * and re-computes state estimates forward from that point.
 *
 * Based on: Schenato (2008), "Optimal estimation in networked control
 * systems subject to random delay and packet drop"
 *
 * @param sp         Smith predictor
 * @param y_meas     Out-of-sequence measurement
 * @param timestamp  Actual measurement timestamp
 * @return           0 on success, -1 if measurement too old
 */
int smith_handle_oos_measurement(SmithPredictor *sp,
                                  const double *y_meas, double timestamp);

/* ============================================================================
 * L5: Delay Generation & Simulation
 * ============================================================================ */

/**
 * delay_generate — generate a synthetic network delay sample
 *
 * Produces a single delay value from the specified statistical model.
 * This enables Monte Carlo simulation of cloud control performance
 * under various network conditions.
 *
 * @param model  Delay model type
 * @param param1 Model parameter (e.g., mean for exponential, min for uniform)
 * @param param2 Model parameter (e.g., rate for exponential, max for uniform)
 * @return       Generated delay in microseconds (>= 0)
 */
double delay_generate(DelayModelType model, double param1, double param2);

/**
 * delay_trace_generate — generate a sequence of delay samples
 *
 * Fills a pre-allocated array with delay samples from the specified
 * model. Useful for simulating network conditions over a control scenario.
 *
 * @param model   Delay model type
 * @param param1  Model parameter 1
 * @param param2  Model parameter 2
 * @param count   Number of samples to generate
 * @param out     Output array (pre-allocated, size >= count)
 */
void delay_trace_generate(DelayModelType model, double param1, double param2,
                           int count, double *out);

/**
 * delay_estimate_online — online delay estimation from measurements
 *
 * Maintains a running exponential moving average (EMA) of network delay
 * for adaptive compensation.
 *
 * @param current_ema    Current EMA estimate (updated in place)
 * @param new_measurement New delay measurement
 * @param alpha          Smoothing factor [0,1] (higher = more weight on new)
 * @return               Updated EMA
 */
double delay_estimate_online(double *current_ema, double new_measurement,
                              double alpha);

/**
 * delay_jitter_compute — compute delay jitter from a trace
 *
 * Jitter is defined as the mean absolute difference between consecutive
 * delay measurements: J = mean(|d_i - d_{i-1}|).
 *
 * @param delays  Array of delay measurements
 * @param count   Number of measurements
 * @return        Mean jitter in microseconds
 */
double delay_jitter_compute(const double *delays, int count);

/**
 * delay_packet_loss_correlate — correlate packet loss with delay
 *
 * In many networks, packet loss increases with delay (buffer overflow).
 * This function computes the correlation between delay and loss events.
 *
 * @param delays    Array of delay measurements
 * @param lost      Array of loss indicators (1 = lost, 0 = received)
 * @param count     Number of samples
 * @return          Pearson correlation coefficient [-1, 1]
 */
double delay_packet_loss_correlate(const double *delays, const int *lost,
                                     int count);

/* ============================================================================
 * L2: Combined Delay-Compensated Control
 * ============================================================================ */

/**
 * delay_compensated_step — execute one control step with delay compensation
 *
 * Combines delay modeling, compensation, and control computation into
 * a single step for cloud control simulation.
 *
 * Flow:
 *   1. Generate/measure network delay
 *   2. Apply delay compensation (Smith predictor)
 *   3. Compute control action using compensated state
 *   4. Simulate delay in control application
 *   5. Advance plant model
 *   6. Update performance metrics
 *
 * @param ccs            Cloud control system
 * @param sp             Smith predictor (NULL to skip compensation)
 * @param measurement    Latest plant measurement
 * @param ts             Measurement timestamp
 * @param dt             Control time step
 * @param strategy       Compensation strategy to use
 * @return               0 on success
 */
int delay_compensated_step(CloudControlSystem *ccs, SmithPredictor *sp,
                            const double *measurement, double ts,
                            double dt, CompensationStrategy strategy);


/* Extended analysis functions */
int delay_fit_exponential(const double *delays, int count, double *rate_out);
int delay_fit_normal(const double *delays, int count, double *mean_out, double *std_out);
double delay_predict_bound(const DelayStatistics *stats, double confidence);
double delay_adaptive_threshold(const DelayStatistics *stats, double safety_factor);
int delay_is_anomalous(const DelayStatistics *baseline, double new_delay, double sigma_threshold);
double delay_network_quality_score(const DelayStatistics *stats, double target_latency_us);
double delay_control_cost(const CloudControlSystem *ccs, const double *delays, int count);
int delay_process_batch(const double *delays, int count, DelayStatistics *stats, double *fitted_rate);

#endif /* CLOUD_CONTROL_DELAY_H */
