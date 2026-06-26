#ifndef CONSENSUS_ALGORITHMS_H
#define CONSENSUS_ALGORITHMS_H
#include "consensus_types.h"
#include "consensus_graph.h"
/* ============================================================================
 * Consensus Algorithms (L5 - Algorithms/Methods)
 *
 * Method catalog for distributed consensus:
 *   A1: Continuous-time Laplacian flow (Olfati-Saber 2004)
 *   A2: Discrete-time weighted averaging (DeGroot 1974)
 *   A3: Gossip algorithms (Boyd-Ghosh-Prabhakar-Shah 2006)
 *   A4: Push-sum for directed graphs (Kempe-Dobra-Gehrke 2003)
 *   A5: Ratio consensus for digraphs
 *   A6: Max/min consensus
 *   A7: Finite-time consensus (Sundaram-Hadjicostis 2007)
 *   A8: Event-triggered consensus (Dimarogonas-Frazzoli-Johansson 2012)
 *
 * Weight design strategies (L5):
 *   W1: Metropolis-Hastings: w_ij = 1/(max(d_i,d_j)+1)
 *   W2: Max-degree: w_ij = 1/(d_max+1)
 *   W3: Laplacian-based: w_ij = epsilon * a_ij
 *   W4: Optimal: minimize rho(W - 11^T/N)
 *   W5: Lazy Metropolis: w_ii = 1 - sum_{j!=i} w_ij
 * ============================================================================ */

/* Gossip algorithm (L5) - pairwise randomized consensus */
typedef struct {
    int active_i;
    int active_j;
    double probability;
    double mixing_rate;
    int total_gossips;
} GossipState;

GossipState* gossip_state_create(void);
void gossip_state_free(GossipState* gs);
int consensus_gossip_step(ConsensusState* cs, GossipState* gs);
int consensus_gossip_simulate(ConsensusState* cs, int num_gossips);
double consensus_gossip_expected_convergence_time(const ConsensusGraph* graph);

/* Push-sum algorithm for directed graphs (L5) */
typedef struct {
    double* push_weights;    /* N*N column-stochastic */
    double* sum_weights;     /* auxiliary sum variable per agent */
    double* sum_estimates;   /* s_i per agent */
} PushSumState;

PushSumState* pushsum_state_create(int N, const ConsensusGraph* graph);
void pushsum_state_free(PushSumState* ps);
int consensus_pushsum_step(ConsensusState* cs, PushSumState* ps);
int consensus_pushsum_simulate(ConsensusState* cs, int max_steps);

/* Max/Min consensus (L5) */
int consensus_max_step(ConsensusState* cs);
int consensus_min_step(ConsensusState* cs);
int consensus_max_simulate(ConsensusState* cs, int max_steps);
int consensus_min_simulate(ConsensusState* cs, int max_steps);
void consensus_compute_max_consensus_value(const double* values, int N,
                                            const ConsensusGraph* graph,
                                            double* result, int* iterations);
void consensus_compute_min_consensus_value(const double* values, int N,
                                            const ConsensusGraph* graph,
                                            double* result, int* iterations);

/* Ratio consensus for directed graphs (L5) */
int consensus_ratio_step(ConsensusState* cs);
int consensus_ratio_simulate(ConsensusState* cs, int max_steps);

/* Finite-time consensus (L5) - exact in finite steps using observability */
int consensus_finite_time_minimal_polynomial(const ConsensusGraph* graph,
                                              double* coefficients, int* degree);
int consensus_finite_time_step(ConsensusState* cs,
                                const double* poly_coeffs, int poly_degree);
bool consensus_finite_time_is_feasible(const ConsensusGraph* graph, int* K);

/* Event-triggered consensus (L5) */
typedef struct {
    double* last_broadcast_state;  /* N x d */
    double* measurement_error;     /* N x d */
    double threshold;              /* sigma */
    double trigger_rate;           /* fraction of agents triggering */
    int total_triggers;
} EventTriggerState;

EventTriggerState* event_trigger_state_create(int N, int dim, double threshold);
void event_trigger_state_free(EventTriggerState* et);
int consensus_event_triggered_step(ConsensusState* cs, EventTriggerState* et, double dt);
bool consensus_event_triggered_should_trigger(const ConsensusAgent* agent,
                                               const EventTriggerState* et, int agent_id);

/* Weight optimization (L5) */
double consensus_optimal_constant_edge_weight(const ConsensusGraph* graph);
int consensus_find_fastest_mixing_chain(const ConsensusGraph* graph,
                                         double* optimal_weights);
double consensus_fastest_mixing_rate(const ConsensusGraph* graph);

/* Distributed estimation / sensor fusion via consensus (L7 application) */
int consensus_distributed_estimation_init(ConsensusState* cs,
                                           const double* measurements,
                                           const double* noise_variance);
int consensus_distributed_estimation_step(ConsensusState* cs);

#endif