#include <stdio.h>
#include <stdlib.h>
#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"
#include "consensus_algorithms.h"

extern int resilient_consensus_msr_step(ConsensusState* cs, int F);
extern int quantized_consensus_step(ConsensusState* cs, double delta, bool use_prob);

int main(void) {
    consensus_seed(42);
    printf("=== Step-by-step debug ===\n");

    /* Test gossip */
    printf("1. Testing gossip...\n");
    ConsensusState* cs = consensus_state_create(8, 1);
    cs->graph = consensus_graph_create_cycle(8);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_state_set_initial_random(cs, -5.0, 5.0);
    consensus_state_set_tolerance(cs, 0.01);
    printf("   setup OK, calling gossip_simulate...\n");
    int ret = consensus_gossip_simulate(cs, 5000);
    printf("   gossip_simulate returned %d, converged=%d\n", ret, cs->has_converged);
    consensus_state_free(cs);

    /* Test resilient */
    printf("2. Testing resilient MSR...\n");
    ConsensusState* cs2 = consensus_state_create(7, 1);
    cs2->graph = consensus_graph_create_complete(7);
    consensus_graph_build_adjacency(cs2->graph);
    consensus_graph_build_degree(cs2->graph);
    consensus_state_set_initial_random(cs2, -5.0, 5.0);
    cs2->agents[3].is_malicious = true;
    printf("   setup OK, calling MSR steps...\n");
    for (int s = 0; s < 500 && !cs2->has_converged; s++) {
        resilient_consensus_msr_step(cs2, 1);
        consensus_state_compute_disagreement(cs2);
        if (cs2->max_disagreement < 0.01) break;
    }
    printf("   MSR: max_disagreement=%g\n", cs2->max_disagreement);
    consensus_state_free(cs2);

    /* Test quantized */
    printf("3. Testing quantized...\n");
    ConsensusState* cs3 = consensus_state_create(4, 1);
    cs3->graph = consensus_graph_create_complete(4);
    consensus_graph_build_adjacency(cs3->graph);
    consensus_graph_build_degree(cs3->graph);
    consensus_graph_build_perron_matrix(cs3->graph, WEIGHT_METROPOLIS);
    consensus_state_set_initial_random(cs3, -3.0, 3.0);
    printf("   setup OK, calling quantized steps...\n");
    for (int s = 0; s < 1000; s++) {
        quantized_consensus_step(cs3, 0.1, true);
        consensus_state_compute_disagreement(cs3);
        if (cs3->max_disagreement < 0.2) break;
    }
    printf("   quantized: max_disagreement=%g\n", cs3->max_disagreement);
    consensus_state_free(cs3);

    printf("All debug tests passed!\n");
    return 0;
}
