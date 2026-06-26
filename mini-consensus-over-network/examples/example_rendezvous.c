#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"
#include "consensus_algorithms.h"

/* Example 2: Multi-Robot Rendezvous in 2D
 *
 * 8 robots scattered in a 2D plane must meet at a common point.
 * Each robot only communicates with robots within a radius (proximity graph).
 * Using continuous consensus on x and y coordinates,
 * all robots converge to their initial positions' average.
 */
int main(void) {
    consensus_seed(99);
    printf("=== Example 2: Multi-Robot Rendezvous (2D) ===\n\n");

    int N = 8;
    ConsensusState* cs = consensus_state_create(N, 2);
    cs->graph = consensus_graph_create(N, TOPO_UNDIRECTED);

    /* Initial positions: randomly scattered in [-10, 10] x [-10, 10] */
    double init_x[8], init_y[8];
    double sum_x = 0, sum_y = 0;
    printf("Initial positions:\n");
    for (int i = 0; i < N; i++) {
        init_x[i] = consensus_rand_uniform(-10, 10);
        init_y[i] = consensus_rand_uniform(-10, 10);
        sum_x += init_x[i]; sum_y += init_y[i];
        cs->state_matrix[i * 2] = init_x[i];
        cs->state_matrix[i * 2 + 1] = init_y[i];
        cs->agents[i].state[0] = init_x[i];
        cs->agents[i].state[1] = init_y[i];
        printf("  Robot %d: (%.2f, %.2f)\n", i, init_x[i], init_y[i]);
    }

    double avg_x = sum_x / N, avg_y = sum_y / N;
    printf("\nRendezvous target (average): (%.4f, %.4f)\n", avg_x, avg_y);

    /* Build proximity graph: connect robots within radius 8.0 */
    double* positions = (double*)malloc((size_t)N * 2 * sizeof(double));
    for (int i = 0; i < N; i++) {
        positions[i * 2] = init_x[i];
        positions[i * 2 + 1] = init_y[i];
    }
    consensus_graph_set_proximity(cs->graph, positions, 2, 8.0);
    free(positions);
    printf("Communication links: %d edges\n", cs->graph->n_edges / 2);
    printf("Graph connected: %s\n",
           consensus_graph_is_connected(cs->graph) ? "YES" : "NO");

    /* Run consensus */
    consensus_state_set_tolerance(cs, 0.01);
    cs->step_size = 0.05;
    consensus_continuous_simulate(cs, 200.0, 0.05);

    printf("\nFinal positions (t=%.2f):\n", cs->time_elapsed);
    for (int i = 0; i < N; i++) {
        printf("  Robot %d: (%.4f, %.4f) dist_to_target=%.4f\n",
               i, cs->state_matrix[i * 2], cs->state_matrix[i * 2 + 1],
               sqrt(pow(cs->state_matrix[i * 2] - avg_x, 2) +
                    pow(cs->state_matrix[i * 2 + 1] - avg_y, 2)));
    }

    printf("\nRendezvous achieved: %s\n",
           cs->has_converged ? "YES" : "NO");
    printf("Final max disagreement: %.4f\n", cs->max_disagreement);

    consensus_state_free(cs);
    printf("\nExample 2 complete.\n");
    return 0;
}