/**
 * example3_multi_loop_scheduling.c
 *
 * Demonstrates bandwidth scheduling for multiple control loops
 * sharing a common communication network.
 *
 * Scenario: A smart grid with 5 distributed control loops:
 *   - Loop 0: Generator frequency regulation (critical)
 *   - Loop 1: Voltage regulator at substation
 *   - Loop 2: Battery storage charge control
 *   - Loop 3: Solar inverter power tracking
 *   - Loop 4: Load demand response
 *
 * All loops share a 2 Mbps SCADA network with TDMA scheduling.
 * The scheduler allocates bandwidth based on each loop's
 * stability requirements (Data Rate Theorem).
 *
 * Real-world relevance:
 *   - Smart grid SCADA (IEC 61850) with bandwidth-limited WAN
 *   - NASA Deep Space Network: multiple spacecraft share limited bandwidth
 *   - Toyota production line: multiple robots share EtherCAT
 *   - Fukushima Daiichi: post-accident monitoring with degraded comms
 *
 * Key findings:
 *   1. Proportional allocation based on unstable eigenvalues
 *      is optimal for maximizing minimum stability margin
 *   2. Adaptive scheduling can respond to changing conditions
 *   3. TOD protocol handles contention optimally
 */

#include "blc_core.h"
#include "blc_scheduling.h"
#include "blc_datarate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/** Loop definitions for the smart grid scenario */
typedef struct {
    const char* name;
    double      unstable_ev;  /** Unstable eigenvalue magnitude */
    const char* application;
    const char* risk;          /** Consequence of instability */
} LoopSpec;

static LoopSpec loops_spec[5] = {
    {"Gen Freq",    3.0, "Generator frequency control",
     "Grid blackout (NERC CIP-014)"},
    {"Voltage Reg", 1.5, "Substation voltage regulation",
     "Equipment damage, brownout"},
    {"Battery BMS", 0.8, "Battery storage charge control",
     "Thermal runaway risk"},
    {"Solar MPPT",  0.5, "Solar inverter power tracking",
     "Efficiency loss, revenue"},
    {"Demand Resp", 0.3, "Load demand response",
     "Customer discomfort"}
};

int main(void) {
    printf("=== Example 3: Multi-Loop Bandwidth Scheduling ===\n\n");
    printf("Smart Grid SCADA — 5 control loops, 2 Mbps shared bus\n\n");

    /** Create loop descriptors */
    BLCPlant plants[5];
    BLCLoopDescriptor loops[5];

    for (int i = 0; i < 5; i++) {
        memset(&plants[i], 0, sizeof(BLCPlant));
        plants[i].n_states = 2;
        plants[i].eigenvalues[0] = loops_spec[i].unstable_ev;
        plants[i].eigenvalues[1] = -2.0;  /** Stable eigenvalue */

        blc_loop_init(&loops[i], i, &plants[i]);

        printf("Loop %d: %s\n", i, loops_spec[i].name);
        printf("  Application: %s\n", loops_spec[i].application);
        printf("  Unstable λ:  %.1f\n", loops_spec[i].unstable_ev);
        printf("  Required rate: %.2f bps\n", loops[i].required_rate);
        printf("  Priority:    %.2f\n", loops[i].priority);
        printf("  Risk:        %s\n\n", loops_spec[i].risk);
    }

    /** --- Static Proportional Allocation --- */
    printf("=== Static Proportional Allocation ===\n\n");

    BLCBandwidthAllocator alloc;
    blc_alloc_init(&alloc, 2e6, loops, 5);  /** 2 Mbps total */
    blc_alloc_proportional(&alloc);

    printf("Total bandwidth: %.2f Mbps\n", alloc.total_bandwidth / 1e6);
    printf("Allocated:       %.2f Mbps\n", alloc.allocated_bandwidth / 1e6);
    printf("Fairness (Jain): %.4f\n", alloc.fairness_index);
    printf("Stability margin: %.4f\n", alloc.stability_margin);
    printf("Starving loops:  %d\n\n", alloc.starving_loops);

    blc_alloc_print(&alloc, stdout);

    /** Check sufficiency for each loop */
    printf("\nRate Sufficiency:\n");
    printf("%-15s %12s %12s %s\n", "Loop", "Required", "Allocated", "Status");
    printf("--------------- ------------ ------------ ------\n");
    for (int i = 0; i < 5; i++) {
        bool sufficient = loops[i].required_rate <= alloc.loops[i].allocated_rate;
        printf("%-15s %12.2f %12.2f %s\n",
               loops_spec[i].name,
               loops[i].required_rate,
               alloc.loops[i].allocated_rate,
               sufficient ? "OK ✓" : "LOW ⚠");
    }

    /** --- Max-Min Fair Allocation --- */
    printf("\n=== Max-Min Fair Allocation ===\n\n");

    BLCBandwidthAllocator alloc_mm;
    blc_alloc_init(&alloc_mm, 2e6, loops, 5);
    blc_alloc_maxmin_fair(&alloc_mm);

    printf("Fairness (Jain): %.4f (vs %.4f proportional)\n",
           alloc_mm.fairness_index, alloc.fairness_index);
    printf("Stability margin: %.4f\n", alloc_mm.stability_margin);
    blc_alloc_print(&alloc_mm, stdout);

    /** --- TDMA Schedule Construction --- */
    printf("\n=== TDMA Schedule ===\n\n");

    BLCTDMASchedule tdma;
    blc_tdma_build(&tdma, loops, 5, 0.1, 0.001);  /** 100ms cycle */

    printf("Cycle time: %.3f sec\n", tdma.cycle_time);
    printf("Number of slots: %d\n", tdma.n_slots);
    printf("Slot duration: %.4f sec\n", tdma.slot_duration);
    printf("Guard time: %.4f sec\n", tdma.guard_time);
    printf("Utilization: %.1f%%\n\n", tdma.utilization * 100.0);

    /** Check MATI feasibility */
    printf("MATI Feasibility:\n");
    double Ac[4] = {-1.0, 0.0, 0.0, -2.0};  /** Stable A_c */
    for (int i = 0; i < 5; i++) {
        double mati = blc_loop_compute_mati(&loops[i], Ac, 2, 0.5);
        printf("  Loop %d (%s): MATI=%.4f sec, cycle=%.4f sec → %s\n",
               i, loops_spec[i].name, mati, tdma.cycle_time,
               tdma.cycle_time <= mati ? "FEASIBLE ✓" : "TOO SLOW ⚠");
        loops[i].mati = mati;
    }

    /** --- TOD Protocol Simulation --- */
    printf("\n=== Try-Once-Discard Contention Resolution ===\n\n");

    BLCTODProtocol tod;
    blc_tod_init(&tod, 5);

    /** Set weights based on priority (unstable eigenvalue magnitude) */
    double tod_weights[5];
    for (int i = 0; i < 5; i++) {
        tod_weights[i] = loops[i].priority;
    }
    blc_tod_set_weights(&tod, tod_weights);

    /** Simulate several contention epochs */
    printf("Epoch  Contenders  Winner  Reason\n");
    printf("-----  ----------  ------  ------\n");

    double test_errors[5][5] = {
        {0.5, 0.2, 0.1, 0.05, 0.02},   /** Normal */
        {0.1, 2.0, 0.3, 0.1, 0.05},    /** Voltage regulator error spike */
        {3.0, 0.2, 0.1, 0.1, 0.1},     /** Generator frequency spike */
        {0.1, 0.1, 0.1, 0.1, 5.0},     /** Demand response anomaly */
        {0.5, 1.5, 1.0, 0.5, 0.8},     /** Multiple elevated errors */
    };

    for (int e = 0; e < 5; e++) {
        int winner = blc_tod_resolve(&tod, test_errors[e], 5);
        printf("%5d  %10d  %6d  %s (w*e=%.2f)\n",
               e, 5, winner,
               loops_spec[winner].name,
               tod_weights[winner] * test_errors[e][winner]);
        blc_tod_transmitted(&tod, winner);
    }

    /** TOD statistics */
    int wins[5], discards[5], contentions;
    blc_tod_stats(&tod, wins, discards, &contentions);
    printf("\nTOD Statistics:\n");
    printf("%-15s %8s %8s %8s\n", "Loop", "Wins", "Loses", "Win%%");
    printf("--------------- -------- -------- --------\n");
    for (int i = 0; i < 5; i++) {
        printf("%-15s %8d %8d %7.1f%%\n",
               loops_spec[i].name, wins[i], discards[i],
               (wins[i]+discards[i]) > 0
               ? 100.0*wins[i]/(double)(wins[i]+discards[i]) : 0.0);
    }
    printf("Total contentions: %d\n", contentions);

    /** --- Adaptive Scheduling --- */
    printf("\n=== Adaptive Bandwidth Scheduling ===\n\n");

    BLCAdaptiveScheduler asched;
    blc_adaptive_init(&asched, 5, 2e6, 1.0, 0.05);

    for (int i = 0; i < 5; i++) {
        blc_adaptive_register_loop(&asched, i, &plants[i],
                                    loops[i].required_rate * 0.5,
                                    loops[i].required_rate * 3.0);
    }

    /** Run adaptation over varying error conditions */
    double error_scenarios[3][5] = {
        {0.1, 0.1, 0.1, 0.1, 0.1},   /** Normal */
        {2.0, 0.2, 0.1, 0.1, 0.1},   /** Generator stress */
        {0.1, 3.0, 0.1, 2.0, 0.1},   /** Voltage + Solar issues */
    };

    for (int s = 0; s < 3; s++) {
        printf("Scenario %d errors: [%.1f, %.1f, %.1f, %.1f, %.1f]\n",
               s, error_scenarios[s][0], error_scenarios[s][1],
               error_scenarios[s][2], error_scenarios[s][3],
               error_scenarios[s][4]);

        double cost = blc_adaptive_epoch(&asched, error_scenarios[s]);
        printf("  Cost=%.4f, Converged=%s\n",
               cost, blc_adaptive_has_converged(&asched) ? "yes" : "no");
        blc_adaptive_print(&asched, stdout);
        printf("\n");
    }

    printf("=== Summary ===\n\n");
    printf("This example demonstrated:\n");
    printf("  1. Proportional bandwidth allocation based on Data Rate Theorem\n");
    printf("  2. Max-min fair allocation ensuring baseline for all loops\n");
    printf("  3. TDMA schedule construction and MATI feasibility checking\n");
    printf("  4. Try-Once-Discard protocol for optimal contention resolution\n");
    printf("  5. Adaptive scheduling responding to changing control errors\n");
    printf("\nKey insight: In bandwidth-limited NCS, allocating bandwidth\n");
    printf("proportional to unstable eigenvalue magnitude maximizes the\n");
    printf("minimum stability margin across all control loops.\n");

    return 0;
}