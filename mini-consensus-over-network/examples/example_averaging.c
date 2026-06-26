#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"

/* Example 1: Distributed Averaging via Consensus
 *
 * Scenario: 10 sensor nodes each measure local temperature.
 * They must compute the global average without a central server.
 * Using consensus, each node communicates only with neighbors
 * and all converge to the true average.
 */
int main(void) {
    consensus_seed(42);
    printf("=== Example 1: Distributed Averaging via Consensus ===\n\n");

    int N = 10;
    ConsensusState* cs = consensus_state_create(N, 1);

    /* Ring topology - sparse but connected */
    cs->graph = consensus_graph_create_cycle(N);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_graph_build_laplacian(cs->graph, LAPLACIAN_COMBINATORIAL);
    consensus_graph_compute_eigenvalues(cs->graph);

    printf("Graph: cycle on %d nodes\n", N);
    printf("Algebraic connectivity lambda_2 = %.4f\n",
           cs->graph->algebraic_connectivity);

    /* Initialize with "measurements" e.g., temperatures */
    double measurements[10] = {23.1, 19.5, 27.3, 21.0, 25.7,
                                18.9, 30.2, 22.4, 26.1, 20.8};
    double true_avg = 0.0;
    for (int i = 0; i < N; i++) true_avg += measurements[i];
    true_avg /= (double)N;

    printf("True average: %.4f\n", true_avg);
    printf("Initial measurements:\n");
    for (int i = 0; i < N; i++) {
        cs->state_matrix[i] = measurements[i];
        cs->agents[i].state[0] = measurements[i];
        cs->agents[i].initial[0] = measurements[i];
        printf("  Node %2d: %.2f\n", i, measurements[i]);
    }

    /* Run continuous-time consensus */
    consensus_state_set_tolerance(cs, 1e-6);
    cs->step_size = 0.01;
    int steps = consensus_continuous_simulate(cs, 100.0, 0.01);

    printf("\nAfter %d steps (t=%.2f):\n", steps, cs->time_elapsed);
    for (int i = 0; i < N; i++)
        printf("  Node %2d: %.6f\n", i, cs->state_matrix[i]);

    printf("\nConverged: %s\n", cs->has_converged ? "YES" : "NO");
    printf("Final disagreement: %.2e\n", cs->max_disagreement);
    printf("Average preservation error: %.2e\n",
           consensus_average_preservation_error(cs));

    /* All nodes should agree on the true average */
    double max_error = 0.0;
    for (int i = 0; i < N; i++) {
        double err = fabs(cs->state_matrix[i] - true_avg);
        if (err > max_error) max_error = err;
    }
    printf("Max error from true average: %.6f\n", max_error);

    consensus_state_free(cs);
    printf("\nExample 1 complete.\n");
    return 0;
}