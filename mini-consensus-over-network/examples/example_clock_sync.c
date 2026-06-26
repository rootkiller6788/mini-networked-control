#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_algorithms.h"

/* Example 3: Distributed Clock Synchronization via Consensus
 *
 * 6 sensor nodes each have a local hardware clock with random skew and offset.
 * Using consensus on the virtual clock parameters, all nodes converge
 * to a common logical time without a central time server.
 *
 * This implements the Average TimeSynch protocol (Schenato & Fiorentin 2011).
 */
typedef struct {
    double skew;
    double offset;
    double virtual_skew;
    double virtual_offset;
    double logical_time;
} SyncNode;

int main(void) {
    consensus_seed(77);
    printf("=== Example 3: Clock Synchronization via Consensus ===\n\n");

    int N = 6;
    SyncNode nodes[6];
    double real_time = 0.0;

    /* Initialize clocks with manufacturing imperfections */
    printf("Initial clock parameters:\n");
    for (int i = 0; i < N; i++) {
        nodes[i].skew = 1.0 + consensus_rand_uniform(-0.005, 0.005);
        nodes[i].offset = consensus_rand_uniform(-0.5, 0.5);
        nodes[i].virtual_skew = 1.0;
        nodes[i].virtual_offset = 0.0;
        nodes[i].logical_time = 0.0;
        printf("  Node %d: skew=%.4f, offset=%.4f\n",
               i, nodes[i].skew, nodes[i].offset);
    }

    /* Build communication graph: ring topology */
    ConsensusGraph* graph = consensus_graph_create_cycle(N);
    consensus_graph_build_adjacency(graph);
    consensus_graph_build_degree(graph);
    consensus_graph_build_perron_matrix(graph, WEIGHT_METROPOLIS);

    printf("\nNetwork: ring topology with %d nodes\n", N);

    /* Synchronization rounds */
    printf("\nSynchronization progress:\n");
    for (int round = 0; round < 50; round++) {
        real_time += 0.1;

        /* Update logical times */
        for (int i = 0; i < N; i++) {
            nodes[i].logical_time =
                nodes[i].virtual_skew * (nodes[i].skew * real_time +
                                          nodes[i].offset) +
                nodes[i].virtual_offset;
        }

        /* Consensus on virtual skew (rate compensation) */
        double new_skew[6], new_offset[6];
        for (int i = 0; i < N; i++) {
            new_skew[i] = 0.0;
            new_offset[i] = 0.0;
            for (int j = 0; j < N; j++) {
                double w = graph->perron_matrix[i * N + j];
                new_skew[i] += w * nodes[j].virtual_skew;
                new_offset[i] += w * nodes[j].logical_time;
            }
        }
        for (int i = 0; i < N; i++) {
            nodes[i].virtual_skew = new_skew[i];
            nodes[i].virtual_offset += 0.1 * (new_offset[i] -
                                               nodes[i].logical_time);
        }

        if (round % 10 == 0) {
            double t_min = nodes[0].logical_time;
            double t_max = nodes[0].logical_time;
            for (int i = 1; i < N; i++) {
                if (nodes[i].logical_time < t_min) t_min = nodes[i].logical_time;
                if (nodes[i].logical_time > t_max) t_max = nodes[i].logical_time;
            }
            printf("  Round %2d: real_t=%.1f, logical range=[%.4f, %.4f], "
                   "error=%.4f\n",
                   round, real_time, t_min, t_max, t_max - t_min);
        }
    }

    /* Final state */
    printf("\nFinal logical times:\n");
    for (int i = 0; i < N; i++)
        printf("  Node %d: %.6f\n", i, nodes[i].logical_time);

    double t_min = nodes[0].logical_time;
    double t_max = nodes[0].logical_time;
    for (int i = 1; i < N; i++) {
        if (nodes[i].logical_time < t_min) t_min = nodes[i].logical_time;
        if (nodes[i].logical_time > t_max) t_max = nodes[i].logical_time;
    }
    printf("Sync error (max-min): %.6f seconds\n", t_max - t_min);

    consensus_graph_free(graph);
    printf("\nExample 3 complete.\n");
    return 0;
}