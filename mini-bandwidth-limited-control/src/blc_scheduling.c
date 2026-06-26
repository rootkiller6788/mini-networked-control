/**
 * blc_scheduling.c — Bandwidth Scheduling and Resource Allocation
 *
 * Implementation of bandwidth scheduling algorithms for multi-loop
 * networked control systems:
 *  - TDMA schedule construction and slot management
 *  - Priority-based proportional bandwidth allocation
 *  - Max-min fair allocation
 *  - Adaptive bandwidth scheduling with water-filling
 *  - Try-Once-Discard (TOD) medium access protocol
 *
 * Key results:
 *  - Walsh, Ye, Bushnell (2002): MATI bounds for NCS stability
 *  - Branicky, Phillips, Zhang (2002): optimal bandwidth allocation
 *    maximizes the minimum stability margin
 *  - Walsh & Ye (2001): TOD is optimal for maximizing stability region
 *
 * Knowledge coverage: L5 (Scheduling Algorithms), L7 (Applications)
 */

#include "blc_core.h"
#include "blc_datarate.h"
#include "blc_scheduling.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * Loop Descriptor
 * ================================================================ */

int blc_loop_init(BLCLoopDescriptor* loop, int loop_id,
                   const BLCPlant* plant) {
    if (!loop || !plant) return -1;

    loop->loop_id       = loop_id;
    loop->allocated_rate = 0.0;
    loop->priority       = 1.0;
    loop->last_tx_time   = -1.0;
    loop->current_error  = 0.0;
    loop->instability_risk = 0.0;
    loop->is_critical    = false;

    /** Compute required rate from Data Rate Theorem */
    loop->required_rate = blc_datarate_min_ct(plant);
    if (loop->required_rate < 1.0) loop->required_rate = 1.0;  /** Minimum */

    /** Priority = sum of unstable eigenvalue magnitudes */
    loop->n_unstable = 0;
    loop->max_unstable_ev = 0.0;
    for (int i = 0; i < plant->n_states; i++) {
        if (plant->eigenvalues[i] > 1e-10) {
            loop->n_unstable++;
            if (plant->eigenvalues[i] > loop->max_unstable_ev) {
                loop->max_unstable_ev = plant->eigenvalues[i];
            }
        }
    }
    loop->priority = loop->max_unstable_ev + 0.1;
    if (loop->priority < 0.1) loop->priority = 0.1;

    /** Default MATI */
    loop->mati = 0.1;

    return 0;
}

double blc_loop_compute_mati(const BLCLoopDescriptor* loop,
                              const double* Ac, int n, double alpha) {
    /** MATI bound from Walsh et al. (2002):
     *  τ_MATI ≤ 1 / (16 · ||A_c|| · √n · (1+α)/(1-α))
     *
     *  where α ∈ (0, 1) is the TOD protocol parameter.
     *  For practical purposes, use α = 1/2 → (1+α)/(1-α) = 3. */
    if (!loop || !Ac || n < 1) return 0.001;
    if (alpha <= 0.0 || alpha >= 1.0) alpha = 0.5;

    double Anorm = 0.0;
    for (int i = 0; i < n * n; i++) {
        Anorm += Ac[i] * Ac[i];
    }
    Anorm = sqrt(Anorm);

    double factor = (1.0 + alpha) / (1.0 - alpha);
    double mati = 1.0 / (16.0 * Anorm * sqrt((double)n) * factor);

    /** Practical lower bound */
    if (mati < 0.001) mati = 0.001;
    if (mati > 1.0) mati = 1.0;

    return mati;
}

void blc_loop_update(BLCLoopDescriptor* loop, double error_norm,
                      double time, double last_tx_time) {
    if (!loop) return;

    loop->current_error = error_norm;

    /** Instability risk: time since last TX / MATI */
    if (loop->mati > 0 && last_tx_time >= 0) {
        double elapsed = time - last_tx_time;
        loop->instability_risk = elapsed / loop->mati;
        if (loop->instability_risk > 1.0) loop->instability_risk = 1.0;
    }
}

/* ================================================================
 * TDMA Schedule
 *
 * Builds a static TDMA schedule where each control loop gets
 * slots proportional to its required data rate.
 *
 * The cycle time T_cycle must satisfy:
 *   T_cycle ≤ min_i MATI_i / guard_factor
 *
 * Each slot duration: T_slot = T_cycle / n_slots
 * Number of slots for loop i: n_i = ceil(required_rate_i / R_slot)
 * where R_slot = bits_per_slot / T_cycle
 * ================================================================ */

int blc_tdma_build(BLCTDMASchedule* tdma, const BLCLoopDescriptor* loops,
                    int n_loops, double cycle_time, double guard_time) {
    if (!tdma || !loops || n_loops < 1 || n_loops > BLC_MAX_SCHEDULE_SLOTS
        || cycle_time <= 0.0) return -1;

    tdma->n_slots        = n_loops;
    tdma->slot_duration   = cycle_time / (double)n_loops;
    tdma->cycle_time      = cycle_time;
    tdma->guard_time      = guard_time;
    tdma->current_slot    = 0;
    tdma->cycle_start_time = 0.0;

    /** Compute total required rate */
    double total_rate = 0.0;
    for (int i = 0; i < n_loops; i++) {
        total_rate += loops[i].required_rate;
    }

    if (total_rate <= 0.0) total_rate = 1.0;

    /** Assign slots proportionally to required rates */
    double cum_fraction[BLC_MAX_SCHEDULE_SLOTS + 1];
    cum_fraction[0] = 0.0;
    for (int i = 0; i < n_loops; i++) {
        cum_fraction[i + 1] = cum_fraction[i]
                              + loops[i].required_rate / total_rate;
    }

    tdma->n_slots = BLC_MAX_SCHEDULE_SLOTS;
    if (tdma->n_slots > n_loops * 4) tdma->n_slots = n_loops * 4;

    int current_owner = 0;
    for (int s = 0; s < tdma->n_slots; s++) {
        double slot_frac = (double)s / (double)tdma->n_slots;
        /** Find which loop owns this fraction */
        while (current_owner < n_loops - 1 &&
               cum_fraction[current_owner + 1] < slot_frac) {
            current_owner++;
        }
        tdma->slot_owner[s] = loops[current_owner].loop_id;
    }

    tdma->utilization = total_rate > 0 ? 1.0 : 0.0;

    return 0;
}

int blc_tdma_next_slot(BLCTDMASchedule* tdma, double current_time) {
    if (!tdma) return -1;

    tdma->current_slot = (tdma->current_slot + 1) % tdma->n_slots;
    if (tdma->current_slot == 0) {
        tdma->cycle_start_time = current_time;
    }

    return tdma->slot_owner[tdma->current_slot];
}

int blc_tdma_current_owner(const BLCTDMASchedule* tdma) {
    if (!tdma || tdma->n_slots == 0) return -1;
    return tdma->slot_owner[tdma->current_slot];
}

double blc_tdma_utilization(const BLCTDMASchedule* tdma) {
    return tdma ? tdma->utilization : 0.0;
}

bool blc_tdma_is_feasible(const BLCTDMASchedule* tdma,
                           const BLCLoopDescriptor* loops, int n_loops) {
    if (!tdma || !loops) return false;

    /** Each loop must get at least one slot within its MATI */
    for (int i = 0; i < n_loops; i++) {
        if (loops[i].mati > 0 && tdma->cycle_time > loops[i].mati) {
            return false;
        }
    }
    return true;
}

/* ================================================================
 * Priority-Based Bandwidth Allocator
 *
 * Proportional allocation: bandwidth_i ∝ w_i
 * where w_i = Σ log₂|λ_ij| for |λ_ij| > 1
 *
 * This is the optimal static allocation (Branicky et al., 2002).
 * ================================================================ */

int blc_alloc_init(BLCBandwidthAllocator* alloc, double total_bandwidth,
                    const BLCLoopDescriptor* loops, int n_loops) {
    if (!alloc || total_bandwidth <= 0.0 || !loops || n_loops < 1
        || n_loops > BLC_MAX_CHANNELS) return -1;

    alloc->n_loops = n_loops;
    alloc->total_bandwidth = total_bandwidth;
    alloc->allocated_bandwidth = 0.0;
    alloc->fairness_index = 0.0;
    alloc->stability_margin = 0.0;
    alloc->starving_loops = 0;

    for (int i = 0; i < n_loops; i++) {
        alloc->loops[i] = loops[i];
    }

    return 0;
}

int blc_alloc_proportional(BLCBandwidthAllocator* alloc) {
    if (!alloc) return -1;
    int N = alloc->n_loops;

    /** Compute weights w_i */
    double weights[BLC_MAX_CHANNELS];
    double total_weight = 0.0;
    for (int i = 0; i < N; i++) {
        weights[i] = alloc->loops[i].priority;
        total_weight += weights[i];
    }

    if (total_weight <= 0.0) {
        /** Equal allocation */
        double each = alloc->total_bandwidth / (double)N;
        for (int i = 0; i < N; i++) {
            alloc->loops[i].allocated_rate = each;
        }
    } else {
        /** Proportional */
        for (int i = 0; i < N; i++) {
            alloc->loops[i].allocated_rate =
                alloc->total_bandwidth * weights[i] / total_weight;
        }
    }

    /** Update summary */
    alloc->allocated_bandwidth = 0.0;
    alloc->starving_loops = 0;
    for (int i = 0; i < N; i++) {
        alloc->allocated_bandwidth += alloc->loops[i].allocated_rate;
        if (alloc->loops[i].allocated_rate < alloc->loops[i].required_rate) {
            alloc->starving_loops++;
        }
    }

    /** Compute Jain's fairness index */
    alloc->fairness_index = blc_alloc_fairness(alloc);

    /** Stability margin = min(allocated/required) */
    double min_margin = DBL_MAX;
    for (int i = 0; i < N; i++) {
        double margin = alloc->loops[i].allocated_rate
                        / (alloc->loops[i].required_rate + 1e-10);
        if (margin < min_margin) min_margin = margin;
    }
    alloc->stability_margin = min_margin;

    return 0;
}

int blc_alloc_maxmin_fair(BLCBandwidthAllocator* alloc) {
    /** Max-min fair allocation:
     *  1. Each loop gets its required rate (if possible)
     *  2. Excess bandwidth is distributed equally
     */
    if (!alloc) return -1;
    int N = alloc->n_loops;
    double B = alloc->total_bandwidth;

    /** Step 1: give each loop its required rate */
    double remaining = B;
    int unsatisfied = 0;
    bool satisfied[BLC_MAX_CHANNELS] = {false};

    for (int i = 0; i < N; i++) {
        double req = alloc->loops[i].required_rate;
        if (req <= remaining / (double)(N - i)) {
            alloc->loops[i].allocated_rate = req;
            remaining -= req;
            satisfied[i] = true;
        } else {
            unsatisfied++;
        }
    }

    /** Step 2: distribute remaining equally among unsatisfied */
    if (unsatisfied > 0 && remaining > 0) {
        double each = remaining / (double)unsatisfied;
        for (int i = 0; i < N; i++) {
            if (!satisfied[i]) {
                alloc->loops[i].allocated_rate = each;
            }
        }
    }

    alloc->starving_loops = 0;
    alloc->allocated_bandwidth = 0.0;
    for (int i = 0; i < N; i++) {
        alloc->allocated_bandwidth += alloc->loops[i].allocated_rate;
        if (alloc->loops[i].allocated_rate < alloc->loops[i].required_rate) {
            alloc->starving_loops++;
        }
    }

    alloc->fairness_index = blc_alloc_fairness(alloc);

    return 0;
}

double blc_alloc_fairness(const BLCBandwidthAllocator* alloc) {
    if (!alloc || alloc->n_loops == 0) return 0.0;
    int N = alloc->n_loops;

    /** Jain's fairness index: J = (Σr_i)² / (N · Σr_i²)
     *  where r_i = allocated_i / required_i */
    double sum_r = 0.0, sum_r2 = 0.0;
    for (int i = 0; i < N; i++) {
        double ri = alloc->loops[i].allocated_rate
                     / (alloc->loops[i].required_rate + 1e-10);
        sum_r  += ri;
        sum_r2 += ri * ri;
    }
    if (sum_r2 < 1e-15) return 1.0;

    double J = (sum_r * sum_r) / ((double)N * sum_r2);
    return J;
}

bool blc_alloc_is_feasible(const BLCBandwidthAllocator* alloc) {
    return alloc ? (alloc->starving_loops == 0) : false;
}

double blc_alloc_get_rate(const BLCBandwidthAllocator* alloc, int loop_id) {
    if (!alloc) return 0.0;
    for (int i = 0; i < alloc->n_loops; i++) {
        if (alloc->loops[i].loop_id == loop_id) {
            return alloc->loops[i].allocated_rate;
        }
    }
    return 0.0;
}

void blc_alloc_print(const BLCBandwidthAllocator* alloc, FILE* stream) {
    if (!alloc) return;
    FILE* out = stream ? stream : stdout;

    fprintf(out, "Bandwidth Allocation Report\n");
    fprintf(out, "Total bandwidth: %.2f bps | Allocated: %.2f bps\n",
            alloc->total_bandwidth, alloc->allocated_bandwidth);
    fprintf(out, "Fairness (Jain): %.4f | Stability margin: %.4f\n",
            alloc->fairness_index, alloc->stability_margin);
    fprintf(out, "Starving loops: %d\n", alloc->starving_loops);
    fprintf(out, "%-6s %12s %12s %12s %10s\n",
            "Loop", "Required", "Allocated", "Priority", "Risk");
    fprintf(out, "------ ------------ ------------ ------------ ----------\n");
    for (int i = 0; i < alloc->n_loops; i++) {
        fprintf(out, "%-6d %12.2f %12.2f %12.2f %10.3f\n",
                alloc->loops[i].loop_id,
                alloc->loops[i].required_rate,
                alloc->loops[i].allocated_rate,
                alloc->loops[i].priority,
                alloc->loops[i].instability_risk);
    }
}

/* ================================================================
 * Adaptive Bandwidth Scheduler
 *
 * Water-filling algorithm for dynamic reallocation:
 *   At each epoch, estimate dCost/dRate for each loop.
 *   Allocate bandwidth to loops with highest marginal benefit.
 *
 * Marginal benefit of bandwidth for loop i:
 *   dV/dR ∝ λ_i · exp(-λ_i / R_i)
 * where λ_i is the unstable eigenvalue and R_i is the bit rate.
 * ================================================================ */

int blc_adaptive_init(BLCAdaptiveScheduler* sched, int n_loops,
                       double total_bandwidth, double epoch_duration,
                       double learning_rate) {
    if (!sched || n_loops < 1 || n_loops > BLC_MAX_CHANNELS
        || total_bandwidth <= 0.0) return -1;

    sched->n_loops         = n_loops;
    sched->total_bandwidth  = total_bandwidth;
    sched->current_cost     = 0.0;
    sched->previous_cost    = 0.0;
    sched->learning_rate    = learning_rate;
    sched->epoch_duration   = epoch_duration;
    sched->epoch_count      = 0;
    sched->converged        = false;

    memset(sched->loops, 0, sizeof(sched->loops));
    memset(sched->min_rate, 0, sizeof(sched->min_rate));
    memset(sched->max_rate, 0, sizeof(sched->max_rate));
    memset(sched->cost_gradient, 0, sizeof(sched->cost_gradient));

    return 0;
}

int blc_adaptive_register_loop(BLCAdaptiveScheduler* sched,
                                int loop_id, const BLCPlant* plant,
                                double min_rate, double max_rate) {
    if (!sched || !plant || loop_id < 0 || loop_id >= BLC_MAX_CHANNELS)
        return -1;

    BLCLoopDescriptor* loop = &sched->loops[loop_id];
    blc_loop_init(loop, loop_id, plant);

    sched->min_rate[loop_id] = min_rate;
    sched->max_rate[loop_id] = max_rate;

    /** Initial allocation: equal share */
    loop->allocated_rate = sched->total_bandwidth / (double)sched->n_loops;

    return 0;
}

double blc_adaptive_epoch(BLCAdaptiveScheduler* sched,
                           const double* current_errors) {
    if (!sched) return 0.0;

    sched->previous_cost = sched->current_cost;
    sched->current_cost  = 0.0;

    /** Compute marginal benefit for each loop */
    for (int i = 0; i < sched->n_loops; i++) {
        double lambda = sched->loops[i].priority;
        double R_i = sched->loops[i].allocated_rate;

        /** Cost = error², marginal benefit ≈ λ·e^{-R/λ} */
        if (R_i > 1e-10 && lambda > 1e-10) {
            sched->cost_gradient[i] = lambda * exp(-R_i / lambda);
        } else {
            sched->cost_gradient[i] = 1.0;
        }

        if (current_errors) {
            sched->current_cost += current_errors[i] * current_errors[i];
        }
    }

    /** Water-filling: allocate to highest gradient */
    /** Sort by gradient descending */
    int order[BLC_MAX_CHANNELS];
    for (int i = 0; i < sched->n_loops; i++) order[i] = i;

    /** Simple bubble sort */
    for (int i = 0; i < sched->n_loops - 1; i++) {
        for (int j = 0; j < sched->n_loops - i - 1; j++) {
            if (sched->cost_gradient[order[j]] < sched->cost_gradient[order[j+1]]) {
                int tmp = order[j];
                order[j] = order[j+1];
                order[j+1] = tmp;
            }
        }
    }

    /** Allocate: each loop gets at least min_rate, gradient decides extra */
    double remaining = sched->total_bandwidth;
    for (int i = 0; i < sched->n_loops; i++) {
        remaining -= sched->min_rate[i];
        sched->loops[i].allocated_rate = sched->min_rate[i];
    }

    /** Distribute remaining to highest gradients */
    for (int idx = 0; idx < sched->n_loops && remaining > 0; idx++) {
        int i = order[idx];
        double extra = remaining * sched->cost_gradient[i];
        double total_grad = 0.0;
        for (int j = idx; j < sched->n_loops; j++) {
            total_grad += sched->cost_gradient[order[j]];
        }
        if (total_grad > 0) extra = remaining * sched->cost_gradient[i] / total_grad;

        double max_extra = sched->max_rate[i] - sched->loops[i].allocated_rate;
        if (extra > max_extra) extra = max_extra;
        if (extra < 0) extra = 0;

        sched->loops[i].allocated_rate += extra;
        remaining -= extra;
    }

    sched->epoch_count++;

    /** Convergence check */
    if (sched->epoch_count > 1) {
        double delta = fabs(sched->current_cost - sched->previous_cost)
                       / (sched->previous_cost + 1e-10);
        if (delta < 1e-4) sched->converged = true;
    }

    return sched->current_cost;
}

double blc_adaptive_get_rate(const BLCAdaptiveScheduler* sched,
                              int loop_id) {
    if (!sched || loop_id < 0 || loop_id >= sched->n_loops) return 0.0;
    return sched->loops[loop_id].allocated_rate;
}

bool blc_adaptive_has_converged(const BLCAdaptiveScheduler* sched) {
    return sched ? sched->converged : false;
}

void blc_adaptive_print(const BLCAdaptiveScheduler* sched, FILE* stream) {
    if (!sched) return;
    FILE* out = stream ? stream : stdout;

    fprintf(out, "Adaptive Scheduler State (epoch %d, %s)\n",
            sched->epoch_count,
            sched->converged ? "converged" : "active");
    fprintf(out, "Total BW: %.2f bps | Cost: %e\n",
            sched->total_bandwidth, sched->current_cost);
    fprintf(out, "%-6s %10s %10s %10s\n",
            "Loop", "Rate", "Gradient", "Min/Max");
    fprintf(out, "------ ---------- ---------- ----------\n");
    for (int i = 0; i < sched->n_loops; i++) {
        fprintf(out, "%-6d %10.2f %10.4f %8.1f/%4.1f\n",
                sched->loops[i].loop_id,
                sched->loops[i].allocated_rate,
                sched->cost_gradient[i],
                sched->min_rate[i], sched->max_rate[i]);
    }
}

/* ================================================================
 * Try-Once-Discard Protocol
 *
 * TOD resolves contention by granting the channel to the node
 * with the largest weighted error. This is provably optimal
 * for maximizing the stability region (Walsh & Ye, 2001).
 *
 * Error weight: w_i = ||C_i|| · e^{λ_i_max · T_s}
 * where λ_i_max is the maximum unstable eigenvalue of loop i.
 * ================================================================ */

int blc_tod_init(BLCTODProtocol* tod, int n_nodes) {
    if (!tod || n_nodes < 1 || n_nodes > BLC_MAX_CHANNELS) return -1;

    tod->n_nodes          = n_nodes;
    tod->last_winner      = -1;
    tod->contention_count  = 0;
    tod->max_error         = 0.0;

    memset(tod->errors, 0, sizeof(tod->errors));
    memset(tod->weights, 0, sizeof(tod->weights));
    memset(tod->wins, 0, sizeof(tod->wins));
    memset(tod->discards, 0, sizeof(tod->discards));

    /** Default: equal weights */
    for (int i = 0; i < n_nodes; i++) {
        tod->weights[i] = 1.0;
    }

    return 0;
}

void blc_tod_set_weights(BLCTODProtocol* tod, const double* weights) {
    if (!tod || !weights) return;
    for (int i = 0; i < tod->n_nodes; i++) {
        tod->weights[i] = weights[i];
    }
}

int blc_tod_resolve(BLCTODProtocol* tod, const double* errors,
                     int n_active) {
    if (!tod || !errors || n_active < 1) return -1;
    if (n_active > tod->n_nodes) n_active = tod->n_nodes;

    tod->contention_count++;

    /** Update errors */
    for (int i = 0; i < n_active; i++) {
        tod->errors[i] = errors[i];
        if (errors[i] > tod->max_error) {
            tod->max_error = errors[i];
        }
    }

    /** Find winner: argmax(w_i * e_i) */
    double max_score = -1.0;
    int winner = -1;
    for (int i = 0; i < n_active; i++) {
        double score = tod->weights[i] * fabs(tod->errors[i]);
        if (score > max_score) {
            max_score = score;
            winner = i;
        }
    }

    /** Update statistics */
    if (winner >= 0) {
        tod->wins[winner]++;
        tod->last_winner = winner;
    }

    /** All others discard */
    for (int i = 0; i < n_active; i++) {
        if (i != winner) {
            tod->discards[i]++;
        }
    }

    return winner;
}

void blc_tod_transmitted(BLCTODProtocol* tod, int node_id) {
    /** Reset error for the winning node after successful transmission */
    if (!tod || node_id < 0 || node_id >= tod->n_nodes) return;
    tod->errors[node_id] = 0.0;
}

void blc_tod_stats(const BLCTODProtocol* tod, int* wins,
                    int* discards, int* contentions) {
    if (!tod) return;
    if (wins) {
        for (int i = 0; i < tod->n_nodes; i++) wins[i] = tod->wins[i];
    }
    if (discards) {
        for (int i = 0; i < tod->n_nodes; i++) discards[i] = tod->discards[i];
    }
    if (contentions) *contentions = tod->contention_count;
}