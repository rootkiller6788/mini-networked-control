/**
 * blc_scheduling.h — Bandwidth Scheduling and Multi-Loop Resource Allocation
 *
 * When multiple control loops share a common communication network,
 * the available bandwidth must be scheduled among them. This is the
 * networked control version of the real-time scheduling problem.
 *
 * Key results:
 *
 * 1. Rate-Monotonic Scheduling for Control (Zhang et al., 2001):
 *    Assign sampling periods proportional to control loop bandwidth needs.
 *    Higher unstable eigenvalue → higher priority → shorter period.
 *
 * 2. Maximum Allowable Transfer Interval (MATI, Walsh et al., 2002):
 *    The maximum time between successive transmissions that still
 *    guarantees stability. The scheduler must ensure each loop's MATI
 *    is not exceeded.
 *
 * 3. Try-Once-Discard (TOD) Protocol (Walsh & Ye, 2001):
 *    In case of contention, the node with the largest error transmits.
 *    This is the optimal policy for maximizing overall stability.
 *
 * 4. Optimal Bandwidth Allocation (Branicky et al., 2002):
 *    For N control loops sharing bandwidth B_total:
 *      maximize Σ α_i log(performance_i)
 *      s.t. Σ bandwidth_i ≤ B_total
 *    The solution allocates bandwidth proportional to the "control value"
 *    (unstable eigenvalue magnitude) of each loop.
 *
 * 5. Adaptive Bandwidth Scheduling (Velasco et al., 2004):
 *    Dynamically reallocate bandwidth based on current control performance.
 *    Loops that are near instability get more bandwidth temporarily.
 *
 * Knowledge coverage: L5 (Scheduling Algorithms), L7 (Applications)
 */

#ifndef BLC_SCHEDULING_H
#define BLC_SCHEDULING_H

#include "blc_core.h"

/* ================================================================
 * Bandwidth Scheduling Structures
 * ================================================================ */

/** Control loop descriptor for scheduling */
typedef struct {
    int      loop_id;           /** Loop identifier */
    double   required_rate;     /** Minimum data rate for stabilization (bps) */
    double   allocated_rate;    /** Currently allocated data rate (bps) */
    double   priority;          /** Scheduling priority (higher = more urgent) */
    double   mati;              /** Maximum Allowable Transfer Interval */
    double   last_tx_time;      /** Time of last successful transmission */
    double   current_error;     /** Current control error norm */
    double   instability_risk;  /** Risk of instability [0, 1] */
    bool     is_critical;       /** Critical loop flag (cannot be skipped) */
    int      n_unstable;        /** Number of unstable eigenvalues */
    double   max_unstable_ev;   /** Maximum unstable eigenvalue magnitude */
} BLCLoopDescriptor;

/** Time-Division Multiple Access (TDMA) schedule.
 *  Divides the communication cycle into fixed time slots,
 *  each assigned to a specific control loop.
 *
 *  The cycle time T_cycle must satisfy:
 *    T_cycle ≤ min_i MATI_i / guard_factor
 *  where guard_factor ≥ 2 accounts for jitter and retransmissions.
 */
typedef struct {
    int      n_slots;            /** Number of slots in one cycle */
    int      slot_owner[BLC_MAX_SCHEDULE_SLOTS]; /** Loop ID assigned to each slot */
    double   slot_duration;      /** Duration of each slot (seconds) */
    double   cycle_time;         /** Total cycle time = n_slots * slot_duration */
    double   guard_time;         /** Guard time between slots */
    double   utilization;        /** Schedule utilization (0 to 1) */
    int      current_slot;       /** Current active slot index */
    double   cycle_start_time;   /** Start time of current cycle */
} BLCTDMASchedule;

/** Priority-based bandwidth allocator.
 *  Uses weighted fair queuing with control-theoretic weights
 *  based on the Data Rate Theorem.
 *
 *  Weight of loop i: w_i = Σ_{j: |λ_ij| > 1} log₂|λ_ij|
 *  (sum of logs of unstable eigenvalue magnitudes)
 */
typedef struct {
    int      n_loops;            /** Number of control loops */
    BLCLoopDescriptor loops[BLC_MAX_CHANNELS]; /** Loop descriptors */
    double   total_bandwidth;    /** Total available bandwidth (bps) */
    double   allocated_bandwidth;/** Sum of all allocations */
    double   fairness_index;    /** Jain's fairness index */
    double   stability_margin;  /** Worst-case stability margin */
    int      starving_loops;    /** Number of loops below required rate */
} BLCBandwidthAllocator;

/** Adaptive bandwidth scheduler with performance feedback.
 *
 *  At each scheduling epoch:
 *   1. Assess each loop's current stability margin
 *   2. Compute bandwidth demand increase/decrease
 *   3. Reallocate using water-filling algorithm
 *
 *  Water-filling: allocate bandwidth to the loop with the
 *  highest marginal control value:
 *    dV_i/dR_i ∝ λ_i · e^{-λ_i · T_s(R_i)}
 *  where T_s(R_i) is the achievable sampling period given rate R_i.
 */
typedef struct {
    BLCLoopDescriptor loops[BLC_MAX_CHANNELS];
    int      n_loops;
    double   total_bandwidth;
    double   min_rate[BLC_MAX_CHANNELS];     /** Minimum acceptable rate per loop */
    double   max_rate[BLC_MAX_CHANNELS];     /** Maximum useful rate per loop */
    double   current_cost;                   /** Current aggregate control cost */
    double   previous_cost;                  /** Previous aggregate control cost */
    double   cost_gradient[BLC_MAX_CHANNELS];/** Marginal benefit of extra bandwidth */
    double   learning_rate;                  /** Gradient step size for reallocation */
    double   epoch_duration;                 /** Time between scheduling decisions */
    int      epoch_count;                    /** Number of scheduling epochs */
    bool     converged;                      /** Allocation convergence flag */
} BLCAdaptiveScheduler;

/** Try-Once-Discard (TOD) medium access control state.
 *
 *  In TOD, when multiple nodes want to transmit, the node with
 *  the largest weighted error gets the channel. Others discard
 *  their packets and try next time.
 *
 *  This is provably optimal for maximizing the stability region
 *  under bandwidth constraints (Walsh & Ye, 2001).
 */
typedef struct {
    int      n_nodes;            /** Number of competing nodes */
    double   errors[BLC_MAX_CHANNELS]; /** Current error at each node */
    double   weights[BLC_MAX_CHANNELS]; /** Priority weights (eigenvalue-based) */
    int      last_winner;        /** Node that won last contention */
    int      contention_count;   /** Number of contention events */
    int      wins[BLC_MAX_CHANNELS]; /** Win count per node */
    int      discards[BLC_MAX_CHANNELS]; /** Discard count per node */
    double   max_error;          /** Maximum observed error */
} BLCTODProtocol;

/* ================================================================
 * Loop Descriptor API
 * ================================================================ */

/** Initialize loop descriptor from plant model.
 *  Computes required_rate from Data Rate Theorem.
 *  @param loop Loop descriptor
 *  @param loop_id Unique loop identifier
 *  @param plant Plant model with eigenvalues computed */
int     blc_loop_init(BLCLoopDescriptor* loop, int loop_id,
                       const BLCPlant* plant);

/** Compute Maximum Allowable Transfer Interval (MATI).
 *  MATI = τ_MATI where stability is guaranteed for any
 *  transmission intervals τ ≤ τ_MATI.
 *
 *  For linear systems, using the result of Walsh et al. (2002):
 *    τ_MATI = 1 / (16·||A_c||·sqrt(n)·((1+α)/(1-α)))
 *  where α relates to the TOD protocol parameter.
 *
 *  @param loop Loop descriptor with plant data
 *  @param Ac Closed-loop matrix A - BK
 *  @param n State dimension
 *  @param alpha Protocol parameter (0 < α < 1)
 *  @return MATI in seconds */
double  blc_loop_compute_mati(const BLCLoopDescriptor* loop,
                               const double* Ac, int n, double alpha);

/** Update loop error and instability risk estimate.
 *  @param loop Loop descriptor
 *  @param error_norm Current control error norm
 *  @param time Current time
 *  @param last_tx_time Last successful transmission time */
void    blc_loop_update(BLCLoopDescriptor* loop, double error_norm,
                         double time, double last_tx_time);

/* ================================================================
 * TDMA Schedule API
 * ================================================================ */

/** Build a static TDMA schedule.
 *  Assigns slots proportionally to required data rates.
 *  @param tdma Schedule structure
 *  @param loops Loop descriptors
 *  @param n_loops Number of loops
 *  @param cycle_time Total cycle time (seconds)
 *  @param guard_time Guard time between slots (seconds)
 *  @return 0 on success */
int     blc_tdma_build(BLCTDMASchedule* tdma, const BLCLoopDescriptor* loops,
                        int n_loops, double cycle_time, double guard_time);

/** Advance TDMA schedule to next slot.
 *  @param tdma Schedule state
 *  @param current_time Current simulation time
 *  @return Loop ID that owns the current slot, -1 if idle */
int     blc_tdma_next_slot(BLCTDMASchedule* tdma, double current_time);

/** Get the current slot owner */
int     blc_tdma_current_owner(const BLCTDMASchedule* tdma);

/** Compute TDMA schedule utilization */
double  blc_tdma_utilization(const BLCTDMASchedule* tdma);

/** Check if all loops are feasible within one cycle.
 *  Returns true if all MATI constraints are satisfied. */
bool    blc_tdma_is_feasible(const BLCTDMASchedule* tdma,
                              const BLCLoopDescriptor* loops, int n_loops);

/* ================================================================
 * Priority-Based Bandwidth Allocation API
 * ================================================================ */

/** Initialize bandwidth allocator.
 *  @param alloc Allocator structure
 *  @param total_bandwidth Total available bandwidth (bps)
 *  @param loops Loop descriptors
 *  @param n_loops Number of control loops
 *  @return 0 on success */
int     blc_alloc_init(BLCBandwidthAllocator* alloc, double total_bandwidth,
                        const BLCLoopDescriptor* loops, int n_loops);

/** Allocate bandwidth proportionally to control-theoretic weights.
 *  Weight_i = Σ log₂(|λ_ij|) for |λ_ij| > 1
 *
 *  This is the optimal static allocation that maximizes
 *  the minimum stability margin across all loops.
 *  @param alloc Allocator state
 *  @return 0 on success, -1 if infeasible */
int     blc_alloc_proportional(BLCBandwidthAllocator* alloc);

/** Allocate bandwidth using max-min fairness.
 *  Each loop gets at least its required rate; excess is distributed
 *  equally to loops that can use more bandwidth.
 *  @param alloc Allocator state
 *  @return 0 on success */
int     blc_alloc_maxmin_fair(BLCBandwidthAllocator* alloc);

/** Compute Jain's fairness index:
 *  J = (Σ r_i)² / (n · Σ r_i²)
 *  where r_i = allocated_i / required_i
 *  J = 1 ⇒ perfectly fair, J = 1/n ⇒ completely unfair
 *  @param alloc Allocator state
 *  @return Fairness index [0, 1] */
double  blc_alloc_fairness(const BLCBandwidthAllocator* alloc);

/** Check feasibility: does every loop get ≥ its required rate? */
bool    blc_alloc_is_feasible(const BLCBandwidthAllocator* alloc);

/** Get bandwidth allocation for a specific loop */
double  blc_alloc_get_rate(const BLCBandwidthAllocator* alloc, int loop_id);

/** Print bandwidth allocation table */
void    blc_alloc_print(const BLCBandwidthAllocator* alloc, FILE* stream);

/* ================================================================
 * Adaptive Bandwidth Scheduling API
 * ================================================================ */

/** Initialize adaptive scheduler.
 *  @param sched Scheduler structure
 *  @param n_loops Number of control loops
 *  @param total_bandwidth Total shared bandwidth (bps)
 *  @param epoch_duration Reallocation interval (seconds)
 *  @param learning_rate Gradient step size
 *  @return 0 on success */
int     blc_adaptive_init(BLCAdaptiveScheduler* sched, int n_loops,
                           double total_bandwidth, double epoch_duration,
                           double learning_rate);

/** Register a control loop with the adaptive scheduler.
 *  @param sched Scheduler
 *  @param loop_id Loop identifier
 *  @param plant Plant model
 *  @param min_rate Minimum acceptable rate
 *  @param max_rate Maximum useful rate
 *  @return 0 on success */
int     blc_adaptive_register_loop(BLCAdaptiveScheduler* sched,
                                    int loop_id, const BLCPlant* plant,
                                    double min_rate, double max_rate);

/** Execute one scheduling epoch (reallocate bandwidth).
 *  Uses water-filling based on marginal control value.
 *  @param sched Scheduler
 *  @param current_errors Array of current control errors per loop
 *  @return Total control cost after reallocation */
double  blc_adaptive_epoch(BLCAdaptiveScheduler* sched,
                            const double* current_errors);

/** Get current allocation for a loop */
double  blc_adaptive_get_rate(const BLCAdaptiveScheduler* sched,
                               int loop_id);

/** Check if allocation has converged */
bool    blc_adaptive_has_converged(const BLCAdaptiveScheduler* sched);

/** Print scheduler state */
void    blc_adaptive_print(const BLCAdaptiveScheduler* sched, FILE* stream);

/* ================================================================
 * Try-Once-Discard Protocol API
 * ================================================================ */

/** Initialize TOD protocol state.
 *  @param tod TOD state
 *  @param n_nodes Number of competing nodes */
int     blc_tod_init(BLCTODProtocol* tod, int n_nodes);

/** Set node weights (typically proportional to unstable eigenvalue magnitude) */
void    blc_tod_set_weights(BLCTODProtocol* tod, const double* weights);

/** Resolve contention: winner has max weight*error.
 *  @param tod TOD state
 *  @param errors Current errors at each node
 *  @param n_active Number of active contenders
 *  @return Winning node ID, -1 if no contention */
int     blc_tod_resolve(BLCTODProtocol* tod, const double* errors,
                          int n_active);

/** Update after successful transmission (reset winner's error) */
void    blc_tod_transmitted(BLCTODProtocol* tod, int node_id);

/** Get TOD statistics */
void    blc_tod_stats(const BLCTODProtocol* tod, int* wins,
                       int* discards, int* contentions);

#endif /* BLC_SCHEDULING_H */