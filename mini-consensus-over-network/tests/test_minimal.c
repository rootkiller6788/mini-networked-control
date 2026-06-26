#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"
#include "consensus_algorithms.h"

int main(void) {
    printf("Test 1: agent_create\n");
    ConsensusAgent* a = consensus_agent_create(0, "test", 2);
    assert(a != NULL);
    consensus_agent_free(a);
    printf("  OK\n");

    printf("Test 2: vector\n");
    ConsensusVector* v = consensus_vector_create(3);
    assert(v != NULL);
    consensus_vector_free(v);
    printf("  OK\n");

    printf("Test 3: graph path\n");
    ConsensusGraph* g = consensus_graph_create_path(3);
    assert(g != NULL);
    printf("  Nodes: %d, Edges: %d\n", g->n_nodes, g->n_edges);
    consensus_graph_free(g);
    printf("  OK\n");

    printf("Test 4: consensus state\n");
    ConsensusState* cs = consensus_state_create(3, 1);
    assert(cs != NULL);
    cs->graph = consensus_graph_create_path(3);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_graph_build_laplacian(cs->graph, LAPLACIAN_COMBINATORIAL);
    consensus_state_set_initial_random(cs, -1.0, 1.0);
    int ret = consensus_continuous_step(cs, 0.01);
    printf("  step ret=%d\n", ret);
    consensus_continuous_simulate(cs, 1.0, 0.01);
    printf("  converged=%d, max_disagreement=%g\n", cs->has_converged, cs->max_disagreement);
    consensus_state_free(cs);
    printf("  OK\n");

    printf("Test 5: discrete consensus\n");
    ConsensusState* cs2 = consensus_state_create(3, 1);
    cs2->graph = consensus_graph_create_cycle(3);
    consensus_graph_build_adjacency(cs2->graph);
    consensus_graph_build_degree(cs2->graph);
    consensus_graph_build_perron_matrix(cs2->graph, WEIGHT_METROPOLIS);
    consensus_state_set_initial_random(cs2, -2.0, 2.0);
    consensus_state_set_protocol(cs2, PROTO_DISCRETE_TIME);
    consensus_discrete_simulate(cs2, 100);
    printf("  converged=%d\n", cs2->has_converged);
    consensus_state_free(cs2);
    printf("  OK\n");

    printf("All minimal tests passed!\n");
    return 0;
}
