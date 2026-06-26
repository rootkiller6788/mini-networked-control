#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"
#include "consensus_algorithms.h"

/* extern declarations for resilience functions (L8) */
extern int resilient_consensus_msr_step(ConsensusState* cs, int F);
extern int quantized_consensus_step(ConsensusState* cs, double delta, bool use_probabilistic);

static int tests_passed = 0;
static int tests_failed = 0;
#define TEST(name) printf("  TEST: %s ... ", name)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)

/* ============================================================================
 * L1 Tests: Type Definitions
 * ============================================================================ */
static void test_agent_create(void) {
    TEST("agent_create");
    ConsensusAgent* a = consensus_agent_create(0, "test_agent", 3);
    CHECK(a != NULL, "agent allocation");
    CHECK(a->id == 0, "agent id");
    CHECK(a->state_dim == 3, "state dimension");
    CHECK(a->n_neighbors == 0, "initial neighbors");
    CHECK(!a->is_leader, "not leader");
    CHECK(!a->is_malicious, "not malicious");
    consensus_agent_free(a);
}

static void test_agent_state(void) {
    TEST("agent_set_state");
    ConsensusAgent* a = consensus_agent_create(1, "state_test", 2);
    double s[2] = {1.5, -2.5};
    consensus_agent_set_state(a, s, 2);
    CHECK(fabs(a->state[0] - 1.5) < 1e-12, "state[0]");
    CHECK(fabs(a->state[1] + 2.5) < 1e-12, "state[1]");
    consensus_agent_free(a);
}

static void test_vector_operations(void) {
    TEST("vector_operations");
    ConsensusVector* v = consensus_vector_create(4);
    consensus_vector_set(v, 2.0);
    double nrm = consensus_vector_norm(v);
    CHECK(fabs(nrm - 4.0) < 1e-10, "norm of constant-2 vector");
    ConsensusVector* w = consensus_vector_create(4);
    consensus_vector_set(w, 3.0);
    double dot = consensus_vector_dot(v, w);
    CHECK(fabs(dot - 24.0) < 1e-10, "dot product");
    consensus_vector_axpy(1.0, v, w);
    for (int i = 0; i < 4; i++)
        assert(fabs(w->data[i] - 5.0) < 1e-10);
    PASS();
    consensus_vector_free(v);
    consensus_vector_free(w);
}

static void test_matrix_operations(void) {
    TEST("matrix_operations");
    ConsensusMatrix* A = consensus_matrix_create(3, 3);
    consensus_matrix_eye(A);
    CHECK(fabs(consensus_matrix_get(A, 1, 1) - 1.0) < 1e-12, "identity diagonal");
    CHECK(fabs(consensus_matrix_get(A, 0, 1)) < 1e-12, "identity off-diag");
    double tr = consensus_matrix_trace(A);
    CHECK(fabs(tr - 3.0) < 1e-10, "trace of I_3");
    consensus_matrix_free(A);
}

/* ============================================================================
 * L3 Tests: Graph Construction and Laplacian
 * ============================================================================ */
static void test_graph_create(void) {
    TEST("graph_create");
    ConsensusGraph* g = consensus_graph_create_path(5);
    CHECK(g != NULL, "path allocation");
    CHECK(g->n_nodes == 5, "node count");
    CHECK(g->n_edges == 8, "edge count (4 undirected = 8 directed)");
    double lambda2 = consensus_graph_algebraic_connectivity(g);
    CHECK(lambda2 > 0, "algebraic connectivity > 0");
    consensus_graph_free(g);
}

static void test_graph_laplacian(void) {
    TEST("graph_laplacian");
    ConsensusGraph* g = consensus_graph_create_path(3);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    /* L should be:
     * [1  -1   0]
     * [-1  2  -1]
     * [0  -1   1] */
    CHECK(fabs(g->laplacian[0*3+0] - 1.0) < 1e-10, "L[0,0]=1");
    CHECK(fabs(g->laplacian[0*3+1] + 1.0) < 1e-10, "L[0,1]=-1");
    CHECK(fabs(g->laplacian[1*3+1] - 2.0) < 1e-10, "L[1,1]=2");
    /* Row sums should be 0 */
    for (int i = 0; i < 3; i++) {
        double rs = 0;
        for (int j = 0; j < 3; j++) rs += g->laplacian[i*3+j];
        assert(fabs(rs) < 1e-10);
    }
    PASS();
    consensus_graph_free(g);
}

static void test_graph_connectivity(void) {
    TEST("graph_connectivity");
    ConsensusGraph* g_path = consensus_graph_create_path(6);
    CHECK(consensus_graph_is_connected(g_path), "path is connected");
    consensus_graph_free(g_path);

    ConsensusGraph* g_complete = consensus_graph_create_complete(5);
    CHECK(consensus_graph_is_connected(g_complete), "complete is connected");
    double lambda2 = consensus_graph_algebraic_connectivity(g_complete);
    CHECK(lambda2 > 0, "complete lambda_2 > 0");
    CHECK(fabs(lambda2 - 5.0) < 1.0, "complete lambda_2 ~ N");
    consensus_graph_free(g_complete);
}

static void test_graph_eigenvalues(void) {
    TEST("graph_eigenvalues");
    ConsensusGraph* g = consensus_graph_create_cycle(4);
    /* Cycle C4: Laplacian eigenvalues are {0, 2, 2, 4} */
    consensus_graph_compute_eigenvalues(g);
    CHECK(g->n_eigen == 4, "eigenvalue count after compute");
    CHECK(fabs(g->eigenvalues[0]) < 1e-10, "lambda_1 ~ 0");
    CHECK(fabs(g->eigenvalues[1] - 2.0) < 0.5, "lambda_2 ~ 2 for C4");
    CHECK(g->algebraic_connectivity > 0, "C4 is connected");
    consensus_graph_free(g);
}

static void test_graph_generators(void) {
    TEST("graph_generators");
    ConsensusGraph* star = consensus_graph_create_star(5);
    CHECK(consensus_graph_is_connected(star), "star connected");
    CHECK(fabs(star->degree[0] - 4.0) < 0.1, "star center degree = 4");
    consensus_graph_free(star);

    ConsensusGraph* er = consensus_graph_create_erdos_renyi(20, 0.3);
    CHECK(er->n_nodes == 20, "ER node count");
    CHECK(er->n_edges > 0, "ER has edges");
    consensus_graph_free(er);

    ConsensusGraph* sw = consensus_graph_create_small_world(20, 4, 0.1);
    CHECK(sw->n_nodes == 20, "SW node count");
    consensus_graph_free(sw);
}

/* ============================================================================
 * L4 Tests: Consensus Dynamics and Theorems
 * ============================================================================ */
static void test_continuous_consensus(void) {
    TEST("continuous_consensus");
    ConsensusState* cs = consensus_state_create(5, 1);
    cs->graph = consensus_graph_create_path(5);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_graph_build_laplacian(cs->graph, LAPLACIAN_COMBINATORIAL);
    consensus_state_set_initial_random(cs, -5.0, 5.0);
    consensus_state_set_tolerance(cs, 1e-4);
    cs->step_size = 0.01;

    consensus_continuous_simulate(cs, 50.0, 0.01);
    CHECK(cs->has_converged, "continuous consensus converges");
    CHECK(cs->max_disagreement < 1e-3, "disagreement below threshold");
    CHECK(cs->time_elapsed > 0, "elapsed time positive");

    consensus_state_free(cs);
}

static void test_discrete_consensus(void) {
    TEST("discrete_consensus");
    ConsensusState* cs2 = consensus_state_create(5, 1);
    cs2->graph = consensus_graph_create_cycle(5);
    consensus_graph_build_adjacency(cs2->graph);
    consensus_graph_build_degree(cs2->graph);
    consensus_graph_build_perron_matrix(cs2->graph, WEIGHT_METROPOLIS);
    consensus_state_set_initial_random(cs2, -10.0, 10.0);
    consensus_state_set_tolerance(cs2, 1e-4);
    consensus_state_set_protocol(cs2, PROTO_DISCRETE_TIME);

    consensus_discrete_simulate(cs2, 5000);
    CHECK(cs2->has_converged, "discrete consensus converges");

    /* Check average is preserved */
    double avg_err = consensus_average_preservation_error(cs2);
    CHECK(avg_err < 1e-4, "average preservation");

    consensus_state_free(cs2);
}

static void test_disagreement_energy(void) {
    TEST("disagreement_energy");
    ConsensusState* cs = consensus_state_create(4, 1);
    cs->graph = consensus_graph_create_complete(4);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_graph_build_laplacian(cs->graph, LAPLACIAN_COMBINATORIAL);
    consensus_state_set_initial_random(cs, -3.0, 3.0);

    double e0 = consensus_disagreement_energy(cs);
    CHECK(e0 >= 0, "energy non-negative");

    consensus_continuous_step(cs, 0.01);
    double e1 = consensus_disagreement_energy(cs);
    CHECK(e1 <= e0 + 1e-10, "energy non-increasing");

    consensus_state_free(cs);
}

/* ============================================================================
 * L5 Tests: Algorithms
 * ============================================================================ */
static void test_gossip_consensus(void) {
    TEST("gossip_consensus");
    ConsensusState* cs = consensus_state_create(8, 1);
    cs->graph = consensus_graph_create_cycle(8);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_state_set_initial_random(cs, -5.0, 5.0);
    consensus_state_set_tolerance(cs, 0.01);

    int gossips = consensus_gossip_simulate(cs, 5000);
    CHECK(cs->has_converged, "gossip converges");
    CHECK(gossips > 0, "gossip steps positive");
    CHECK(cs->max_disagreement < 0.02, "gossip disagreement small");

    consensus_state_free(cs);
}

static void test_max_consensus(void) {
    TEST("max_consensus");
    double values[5] = {3.0, 7.0, 1.0, 9.0, 4.0};
    ConsensusGraph* g = consensus_graph_create_path(5);
    double result = 0.0;
    int iters = 0;
    consensus_compute_max_consensus_value(values, 5, g, &result, &iters);
    CHECK(fabs(result - 9.0) < 1e-10, "max consensus finds max");
    CHECK(iters <= 5, "converges in <= diam steps");

    consensus_compute_min_consensus_value(values, 5, g, &result, &iters);
    CHECK(fabs(result - 1.0) < 1e-10, "min consensus finds min");

    consensus_graph_free(g);
}

static void test_weight_design(void) {
    TEST("weight_design");
    ConsensusGraph* g = consensus_graph_create_path(4);
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);

    WeightMatrix* wm = weight_matrix_create(4, WEIGHT_METROPOLIS, g);
    CHECK(wm != NULL, "weight matrix created");
    CHECK(wm->is_row_stochastic, "Metropolis is row-stochastic");

    weight_matrix_free(wm);
    consensus_graph_free(g);
}

/* ============================================================================
 * L8 Tests: Advanced Topics
 * ============================================================================ */
static void test_resilient_consensus(void) {
    TEST("resilient_consensus");
    ConsensusState* cs = consensus_state_create(7, 1);
    cs->graph = consensus_graph_create_complete(7);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_state_set_initial_random(cs, -5.0, 5.0);

    /* Make one agent malicious (index 3 keeps its original value) */
    cs->agents[3].is_malicious = true;

    /* Run MSR with F=1 for many iterations */
    for (int s = 0; s < 2000; s++) {
        resilient_consensus_msr_step(cs, 1);
        consensus_state_compute_disagreement(cs);
    }
    /* Normal agents should have converged among themselves.
     * With one malicious agent, max_disagreement (which includes malicious)
     * may remain large, but normal agents should be close. Check that
     * non-malicious agents have converged. */
    double max_normal_diff = 0.0;
    for (int i = 0; i < 7; i++) {
        if (cs->agents[i].is_malicious) continue;
        for (int j = i + 1; j < 7; j++) {
            if (cs->agents[j].is_malicious) continue;
            double diff = fabs(cs->state_matrix[i] - cs->state_matrix[j]);
            if (diff > max_normal_diff) max_normal_diff = diff;
        }
    }
    CHECK(max_normal_diff < 0.1, "normal agents converge under MSR");

    consensus_state_free(cs);
}

static void test_quantized_consensus(void) {
    TEST("quantized_consensus");
    ConsensusState* cs = consensus_state_create(4, 1);
    cs->graph = consensus_graph_create_complete(4);
    consensus_graph_build_adjacency(cs->graph);
    consensus_graph_build_degree(cs->graph);
    consensus_graph_build_perron_matrix(cs->graph, WEIGHT_METROPOLIS);
    consensus_state_set_initial_random(cs, -3.0, 3.0);

    for (int s = 0; s < 1000; s++) {
        quantized_consensus_step(cs, 0.1, true);
        consensus_state_compute_disagreement(cs);
        if (cs->max_disagreement < 0.2) break;
    }
    CHECK(cs->max_disagreement < 1.0, "quantized consensus runs");

    consensus_state_free(cs);
}

/* ============================================================================
 * Main Test Runner
 * ============================================================================ */
int main(void) {
    printf("=== mini-consensus-over-network Test Suite ===\n\n");
    consensus_seed(12345);

    printf("[L1] Definitions:\n");
    test_agent_create();
    test_agent_state();
    test_vector_operations();
    test_matrix_operations();

    printf("\n[L3] Graph / Laplacian:\n");
    test_graph_create();
    test_graph_laplacian();
    test_graph_connectivity();
    test_graph_eigenvalues();
    test_graph_generators();

    printf("\n[L4] Consensus Dynamics:\n");
    test_continuous_consensus();
    test_discrete_consensus();
    test_disagreement_energy();

    printf("\n[L5] Algorithms:\n");
    test_gossip_consensus();
    test_max_consensus();
    test_weight_design();

    printf("\n[L8] Advanced Topics:\n");
    test_resilient_consensus();
    test_quantized_consensus();

    printf("\n=== Results: %d passed, %d failed ===\n",
           tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}