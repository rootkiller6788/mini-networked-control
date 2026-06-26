#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"

/* Example 4: UAV Formation Control via Consensus
 *
 * 5 drones must form a V-formation while maintaining consensus on velocity.
 * Each drone adjusts its position relative to the formation center using
 * consensus-based formation control with predetermined offsets.
 *
 * Reference: Ren & Beard (2008), Ch. 5, "Formation Control"
 */
int main(void) {
    consensus_seed(123);
    printf("=== Example 4: UAV Formation Control ===\n\n");

    int N = 5;
    int dim = 2;

    /* Desired V-formation offsets (relative to formation center) */
    double offsets[5][2] = {
        { 0.0,  0.0},   /* Leader */
        {-3.0, -2.0},   /* Left wing 1 */
        {-6.0, -4.0},   /* Left wing 2 */
        { 3.0, -2.0},   /* Right wing 1 */
        { 6.0, -4.0}    /* Right wing 2 */
    };

    printf("Desired formation (V-shape):\n");
    for (int i = 0; i < N; i++)
        printf("  Drone %d: offset (%.1f, %.1f)\n", i, offsets[i][0], offsets[i][1]);

    /* Initialize positions with random perturbations */
    ConsensusState* cs = consensus_state_create(N, dim);
    cs->graph = consensus_graph_create_complete(N);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_graph_build_laplacian(cs->graph, LAPLACIAN_COMBINATORIAL);

    printf("\nInitial positions (perturbed from desired):\n");
    for (int i = 0; i < N; i++) {
        double px = offsets[i][0] + consensus_rand_uniform(-3, 3);
        double py = offsets[i][1] + consensus_rand_uniform(-3, 3);
        cs->state_matrix[i * dim] = px;
        cs->state_matrix[i * dim + 1] = py;
        cs->agents[i].state[0] = px;
        cs->agents[i].state[1] = py;
        cs->agents[i].initial[0] = px;
        cs->agents[i].initial[1] = py;
        printf("  Drone %d: (%.2f, %.2f)\n", i, px, py);
    }

    /* Consensus on offset-compensated positions:
     * x_i converges to (1/N) sum (x_j(0) - offset_j) + offset_i
     * This achieves the desired formation shape at the rendezvous point. */
    consensus_state_set_tolerance(cs, 0.05);
    cs->step_size = 0.05;

    /* Use compensated consensus: ẋ_i = -sum a_ij ((x_i - o_i) - (x_j - o_j)) */
    for (int step = 0; step < 5000; step++) {
        for (int i = 0; i < N; i++) {
            cs->agents[i].input[0] = 0.0;
            cs->agents[i].input[1] = 0.0;
            for (int j = 0; j < N; j++) {
                double a_ij = cs->graph->adjacency[i * N + j];
                if (a_ij > 0) {
                    cs->agents[i].input[0] -= a_ij *
                        ((cs->state_matrix[i * dim] - offsets[i][0]) -
                         (cs->state_matrix[j * dim] - offsets[j][0]));
                    cs->agents[i].input[1] -= a_ij *
                        ((cs->state_matrix[i * dim + 1] - offsets[i][1]) -
                         (cs->state_matrix[j * dim + 1] - offsets[j][1]));
                }
            }
        }
        for (int i = 0; i < N; i++) {
            cs->state_matrix[i * dim] += cs->agents[i].input[0] * 0.05;
            cs->state_matrix[i * dim + 1] += cs->agents[i].input[1] * 0.05;
        }
    }

    printf("\nFinal positions:\n");
    for (int i = 0; i < N; i++)
        printf("  Drone %d: (%.4f, %.4f)\n",
               i, cs->state_matrix[i * dim], cs->state_matrix[i * dim + 1]);

    /* Check formation: relative positions should match desired offsets */
    printf("\nFormation errors:\n");
    double max_err = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx_actual = cs->state_matrix[i * dim] - cs->state_matrix[j * dim];
            double dy_actual = cs->state_matrix[i * dim + 1] - cs->state_matrix[j * dim + 1];
            double dx_desired = offsets[i][0] - offsets[j][0];
            double dy_desired = offsets[i][1] - offsets[j][1];
            double err = sqrt(pow(dx_actual - dx_desired, 2) +
                             pow(dy_actual - dy_desired, 2));
            if (err > max_err) max_err = err;
        }
    }
    printf("  Max formation error: %.4f\n", max_err);

    consensus_state_free(cs);
    printf("\nExample 4 complete.\n");
    return 0;
}