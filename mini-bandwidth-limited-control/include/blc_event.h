/**
 * blc_event.h — Event-Triggered and Aperiodic Control under Bandwidth Limits
 *
 * Event-triggered control reduces bandwidth consumption by transmitting
 * sensor data only when "interesting" events occur, rather than at a
 * fixed period. This is the dual of data compression: instead of
 * compressing every sample, skip uninformative samples entirely.
 *
 * Fundamental Concepts:
 *
 * 1. Send-on-Delta (SOD): transmit when |x(t) - x_last| > δ
 *    Simple, effective, requires only local comparison.
 *
 * 2. Lebesgue Sampling (Astrom & Bernhardsson, 2002):
 *    Transmit when the state crosses a predefined level set.
 *    For first-order stable systems, Lebesgue sampling with
 *    logarithmic level spacing is optimal.
 *
 * 3. Event-Triggered State Feedback (Tabuada, 2007):
 *    Transmit when ||e(t)|| > σ||x(t)|| where e(t) = x(t) - x(t_k)
 *    and σ is related to the ISS gain of the closed-loop system.
 *
 * 4. Self-Triggered Control (Velasco et al., 2003; Wang & Lemmon, 2009):
 *    The next transmission time is pre-computed at the current one,
 *    eliminating the need for continuous monitoring.
 *
 * 5. Periodic Event-Triggered Control (Heemels et al., 2012):
 *    Check event condition only at sampling instants, combining
 *    the simplicity of periodic sampling with the efficiency of
 *    event-triggering.
 *
 * Knowledge coverage: L5 (Algorithms), L2 (Core Concepts)
 */

#ifndef BLC_EVENT_H
#define BLC_EVENT_H

#include "blc_core.h"

/* ================================================================
 * Event-Triggered Transmission Types
 * ================================================================ */

/** Send-on-Delta event detector
 *
 *  Trigger condition: ||x - x_last||_∞ > δ
 *  or:              ||x - x_last||_2 > δ  (configurable norm)
 *
 *  Delta is the quantization-aware threshold:
 *    δ = α · Δ_max   where Δ_max is the maximum quantization error
 *    and α ∈ [0.5, 2.0] is a tuning parameter.
 *
 *  Smaller δ → more transmissions, better fidelity
 *  Larger δ  → fewer transmissions, coarser control
 */
typedef struct {
    double  delta;              /** Transmission threshold δ */
    int     norm_type;          /** 1=L1, 2=L2, INF=Chebyshev */
    double  x_last[BLC_MAX_STATES]; /** Last transmitted state */
    double  error_norm;         /** Current error norm */
    int     transmissions;      /** Total transmission count */
    int     suppressed;         /** Suppressed samples (saved bandwidth) */
    double  bandwidth_saved;    /** Percentage of bandwidth saved */
    double  min_interval;       /** Minimum inter-transmission interval */
    double  last_tx_time;       /** Time of last transmission */
    bool    forced_tx;          /** Whether last TX was forced */
} BLCSendOnDelta;

/** Lebesgue sampler
 *
 *  Level set: L_i = {x : μ^{i-1} < |x| ≤ μ^i}
 *  where μ > 1 is the level spacing factor.
 *
 *  Transmission occurs when x crosses from one level set to another.
 *  For stable scalar systems with pole λ < 0, the optimal spacing is:
 *    μ = e^{|λ|·T_desired}  where T_desired is the desired
 *    average inter-sample time.
 *
 *  @ref Astrom & Bernhardsson (2002), "Comparison of Riemann and
 *       Lebesgue sampling for first order stochastic systems", IEEE CDC.
 */
typedef struct {
    double  mu;                 /** Level spacing factor μ > 1 */
    double  levels[32];         /** Pre-computed level thresholds */
    int     n_levels;           /** Number of levels */
    int     current_level;      /** Current level index */
    int     last_level;         /** Last transmitted level */
    int     transmissions;      /** Count of transmissions */
    double  avg_interval;       /** Average inter-sample time */
    double  last_tx_time;       /** Nominal sample period for comparison */
} BLCLebesgueSampler;

/** Event-triggered state feedback controller (Tabuada, 2007)
 *
 *  Control law: u(t) = K x(t_k)  for t ∈ [t_k, t_{k+1})
 *
 *  Trigger: t_{k+1} = inf{ t > t_k : ||e(t)|| > σ ||x(t)|| }
 *  where e(t) = x(t_k) - x(t) is the measurement error.
 *
 *  Theorem (Tabuada, 2007): If the system is ISS with respect to
 *  measurement errors, then there exists σ > 0 such that the
 *  event-triggered implementation is asymptotically stable
 *  AND has a positive minimum inter-event time.
 */
typedef struct {
    double  sigma;              /** Trigger parameter σ (ISS-related) */
    double  beta;               /** ISS gain */
    double  lambda_cl;          /** Closed-loop convergence rate */
    double  min_inter_event;    /** Theoretical minimum inter-event time */
    double  last_event_time;    /** Time of last event */
    int     event_count;        /** Number of events triggered */
    double  x_last[BLC_MAX_STATES]; /** State at last event */
    double  e[BLC_MAX_STATES];  /** Current error e(t) */
    double  V_last;             /** Lyapunov function at last event */
    bool    is_iss;             /** ISS property verified */
} BLCEventTriggeredFeedback;

/** Self-triggered controller
 *
 *  At transmission time t_k, compute the next transmission time:
 *    t_{k+1} = t_k + τ(x(t_k))
 *  where τ(·) is a state-dependent holding time function.
 *
 *  For linear systems with quadratic Lyapunov function V(x) = x'Px:
 *    τ(x) = min{ t : V(e^{A_c t} x) ≤ ρ V(x) }
 *  where ρ ∈ (0, 1) is a design parameter.
 *
 *  Advantage: no continuous monitoring needed between transmissions.
 *  Disadvantage: more conservative (larger inter-event times would be possible
 *  with continuous event-triggering).
 */
typedef struct {
    double  rho;                /** Decay parameter ρ ∈ (0, 1) */
    double  tau_min;            /** Minimum holding time (hardware limit) */
    double  tau_max;            /** Maximum holding time (safety) */
    double  next_event_time;    /** Pre-computed next event time */
    double  P[BLC_MAX_STATES][BLC_MAX_STATES]; /** Lyapunov matrix */
    double  A_c[BLC_MAX_STATES][BLC_MAX_STATES]; /** Closed-loop matrix */
    int     event_count;        /** Number of events */
    double  avg_interval;       /** Average inter-event time */
} BLCSelfTriggered;

/* ================================================================
 * Send-on-Delta API
 * ================================================================ */

/** Initialize send-on-delta detector */
int     blc_sod_init(BLCSendOnDelta* sod, double delta, int norm_type);

/** Check if transmission should occur.
 *  @param sod Detector state
 *  @param x Current state vector
 *  @param n Number of states
 *  @param current_time Current simulation time
 *  @param force If true, force transmission regardless
 *  @return true if transmission should occur */
bool    blc_sod_should_transmit(BLCSendOnDelta* sod, const double* x,
                                 int n, double current_time, bool force);

/** Update after transmission (reset last known state) */
void    blc_sod_transmitted(BLCSendOnDelta* sod, const double* x,
                             double current_time);

/** Get bandwidth saving statistics */
void    blc_sod_stats(const BLCSendOnDelta* sod, double* pct_saved,
                       int* tx_count, int* suppressed);

/** Compute optimal delta from quantization step.
 *  δ_opt = k · Δ  where k ∈ [1, 3] trades off fidelity vs bandwidth */
double  blc_sod_optimal_delta(double quant_step, double tradeoff_factor);

/* ================================================================
 * Lebesgue Sampler API
 * ================================================================ */

/** Initialize Lebesgue sampler with logarithmic level spacing.
 *  @param ls Sampler to initialize
 *  @param mu Level spacing factor (> 1)
 *  @param n_levels Number of levels
 *  @param max_val Maximum absolute value to cover */
int     blc_lebesgue_init(BLCLebesgueSampler* ls, double mu,
                           int n_levels, double max_val);

/** Check if Lebesgue crossing occurred.
 *  @param ls Sampler state
 *  @param x Current signal value (scalar)
 *  @param time Current time
 *  @return true if level crossing detected */
bool    blc_lebesgue_check(BLCLebesgueSampler* ls, double x, double time);

/** Get current level index for a value */
int     blc_lebesgue_get_level(const BLCLebesgueSampler* ls, double x);

/** Get Lebesgue sampling statistics */
void    blc_lebesgue_stats(const BLCLebesgueSampler* ls,
                            double* avg_interval, int* tx_count);

/* ================================================================
 * Event-Triggered Feedback API
 * ================================================================ */

/** Initialize event-triggered feedback controller.
 *  @param etf Controller structure
 *  @param sigma Trigger parameter (typically 0.01 to 0.5)
 *  @param Ac Closed-loop system matrix A-BK
 *  @param n State dimension
 *  @return 0 on success, -1 if not ISS */
int     blc_etf_init(BLCEventTriggeredFeedback* etf, double sigma,
                      const double* Ac, int n);

/** Check event trigger condition.
 *  @param etf Controller
 *  @param x Current state
 *  @param n State dimension
 *  @param time Current time
 *  @return true if event (transmission) should be triggered */
bool    blc_etf_check(BLCEventTriggeredFeedback* etf, const double* x,
                       int n, double time);

/** Update after event (reset error to zero) */
void    blc_etf_triggered(BLCEventTriggeredFeedback* etf,
                           const double* x, double time);

/** Compute inter-event time lower bound.
 *  τ_min = (1/λ) · ln(1 + σ·β)
 *  where λ is the closed-loop convergence rate and β is the ISS gain.
 *  @param etf Controller (must be initialized)
 *  @return Minimum possible inter-event time */
double  blc_etf_min_inter_event(const BLCEventTriggeredFeedback* etf);

/** Get event-triggering statistics */
void    blc_etf_stats(const BLCEventTriggeredFeedback* etf,
                       int* event_count, double* avg_interval);

/* ================================================================
 * Self-Triggered Control API
 * ================================================================ */

/** Initialize self-triggered controller.
 *  @param st Controller structure
 *  @param rho Decay parameter (0 < rho < 1)
 *  @param tau_min Minimum holding time in seconds
 *  @param tau_max Maximum holding time in seconds
 *  @param Ac Closed-loop matrix A-BK (n×n, row-major)
 *  @param P Lyapunov matrix (n×n, row-major)
 *  @param n State dimension */
int     blc_self_trig_init(BLCSelfTriggered* st, double rho,
                            double tau_min, double tau_max,
                            const double* Ac, const double* P, int n);

/** Compute next event time based on current state.
 *  @param st Self-triggered controller
 *  @param x Current state
 *  @param n State dimension
 *  @param current_time Current time
 *  @return Next scheduled event time */
double  blc_self_trig_next_event(BLCSelfTriggered* st, const double* x,
                                  int n, double current_time);

/** Get self-triggering statistics */
void    blc_self_trig_stats(const BLCSelfTriggered* st,
                             double* avg_interval, int* event_count);

#endif /* BLC_EVENT_H */