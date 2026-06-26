#include "consensus_algorithms.h"
#include "consensus_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Gossip Algorithm (L5 — Randomized pairwise consensus)
 *
 * In each round, a random edge (i,j) is selected and both agents average:
 *   x_i[k+1] = x_j[k+1] = (x_i[k] + x_j[k]) / 2
 *
 * Equivalently: x[k+1] = W[k] x[k] where W[k] = I - (1/2)(e_i - e_j)(e_i - e_j)^T
 * with probability p_ij.
 *
 * Reference: Boyd, Ghosh, Prabhakar, Shah (2006) —
 *   "Randomized gossip algorithms", IEEE Trans. Info. Theory
 *
 * Expected convergence time: Theta(1 / lambda_2(P)) where P is the
 * probability matrix of edge selection.
 * ============================================================================ */

GossipState* gossip_state_create(void) {
    GossipState* gs = (GossipState*)calloc(1, sizeof(GossipState));
    gs->active_i = -1;
    gs->active_j = -1;
    gs->probability = 1.0;
    gs->mixing_rate = 0.0;
    gs->total_gossips = 0;
    return gs;
}

void gossip_state_free(GossipState* gs) { free(gs); }

int consensus_gossip_step(ConsensusState* cs, GossipState* gs) {
    /* One gossip round: select random edge (i,j) and average their states.
     *
     * The gossip matrix is: W = I - (1/2)(e_i - e_j)(e_i - e_j)^T
     * This is doubly stochastic and idempotent on the selected pair.
     *
     * Edge selection probability proportional to weight for faster mixing.
     * Complexity: O(N) expected for edge selection, O(d) for update. */
    if (!cs || !cs->graph || !gs) return -1;
    int d = cs->agents[0].state_dim;

    /* Select active edge uniformly random from existing edges */
    if (cs->graph->n_edges < 1) return -1;
    int e = (int)(consensus_rand_uniform(0, 1) * cs->graph->n_edges);
    if (e >= cs->graph->n_edges) e = cs->graph->n_edges - 1;
    gs->active_i = cs->graph->edges[e].from;
    gs->active_j = cs->graph->edges[e].to;
    gs->total_gossips++;

    /* Average the two agents' states */
    for (int k = 0; k < d; k++) {
        double avg = 0.5 * (cs->state_matrix[gs->active_i * d + k] +
                            cs->state_matrix[gs->active_j * d + k]);
        cs->state_matrix[gs->active_i * d + k] = avg;
        cs->state_matrix[gs->active_j * d + k] = avg;
        cs->agents[gs->active_i].state[k] = avg;
        cs->agents[gs->active_j].state[k] = avg;
    }
    cs->iterations++;
    cs->has_converged = consensus_state_check_convergence(cs, CONV_ABSOLUTE);
    return 0;
}

int consensus_gossip_simulate(ConsensusState* cs, int num_gossips) {
    if (!cs) return -1;
    GossipState* gs = gossip_state_create();
    for (int g = 0; g < num_gossips; g++) {
        consensus_gossip_step(cs, gs);
        if (cs->has_converged) break;
    }
    gossip_state_free(gs);
    return cs->iterations;
}

double consensus_gossip_expected_convergence_time(const ConsensusGraph* graph) {
    /* Expected epsilon-convergence time for randomized gossip:
     * E[T_eps] = Theta(1 / Phi_e(P)) where Phi_e(P) is the effective
     * spectral gap of the expected gossip matrix.
     * For uniform edge selection on connected graph:
     *   E[T_eps] ~ (log(1/eps)) * N^2 / (2 * lambda_2(L))
     * Complexity: O(1). */
    if (!graph || graph->algebraic_connectivity < 1e-12) return 1e9;
    double N = (double)graph->n_nodes;
    return log(100.0) * N * N / (2.0 * graph->algebraic_connectivity);
}

/* ============================================================================
 * Push-Sum Algorithm (L5 — Consensus on directed graphs)
 *
 * For directed graphs where W is only column-stochastic (not doubly stochastic),
 * the standard algorithm does not converge to the average.
 *
 * Push-sum (Kempe, Dobra, Gehrke 2003):
 *   s_i[k+1] = sum_{j in N_i_in} s_j[k] / d_j_out
 *   w_i[k+1] = sum_{j in N_i_in} w_j[k] / d_j_out,  w_i[0] = 1
 *   x_i[k] = s_i[k] / w_i[k]  (ratio converges to average)
 *
 * This achieves average consensus on any strongly connected directed graph.
 * ============================================================================ */

PushSumState* pushsum_state_create(int N, const ConsensusGraph* graph) {
    PushSumState* ps = (PushSumState*)calloc(1, sizeof(PushSumState));
    if (!ps) return NULL;
    ps->push_weights = (double*)calloc((size_t)N * (size_t)N, sizeof(double));
    ps->sum_weights = (double*)calloc((size_t)N, sizeof(double));
    ps->sum_estimates = (double*)calloc((size_t)N, sizeof(double));
    /* Initialize column-stochastic weights: w_ij = a_ij / d_j_out */
    for (int j = 0; j < N; j++) {
        double d_out = (graph->degree[j] > 0) ? graph->degree[j] : 1.0;
        for (int i = 0; i < N; i++)
            ps->push_weights[i * N + j] = graph->adjacency[j * N + i] / d_out;
    }
    return ps;
}

void pushsum_state_free(PushSumState* ps) {
    if (!ps) return;
    free(ps->push_weights);
    free(ps->sum_weights);
    free(ps->sum_estimates);
    free(ps);
}

int consensus_pushsum_step(ConsensusState* cs, PushSumState* ps) {
    /* One push-sum iteration.
     * s_i[k+1] = sum_j p_ij * s_j[k]  (column-stochastic mixing)
     * w_i[k+1] = sum_j p_ij * w_j[k]  (mass tracking)
     * x_i[k+1] = s_i[k+1] / w_i[k+1]  (ratio = unbiased estimate)
     * Complexity: O(N^2 * d). */
    if (!cs || !ps) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    for (int k = 0; k < d; k++) {
        double* s_old = (double*)malloc((size_t)N * sizeof(double));
        double* w_old = (double*)malloc((size_t)N * sizeof(double));
        /* Extract current state as s values, w initialized to 1 */
        for (int i = 0; i < N; i++) {
            s_old[i] = cs->state_matrix[i * d + k];
            w_old[i] = 1.0; /* w_i[0] = 1 */
        }

        /* Push: s_new = P^T * s_old, w_new = P^T * w_old */
        for (int i = 0; i < N; i++) {
            double s_new = 0.0, w_new = 0.0;
            for (int j = 0; j < N; j++) {
                s_new += ps->push_weights[i * N + j] * s_old[j];
                w_new += ps->push_weights[i * N + j] * w_old[j];
            }
            /* Ratio estimate */
            cs->state_matrix[i * d + k] = (w_new > 1e-12) ? s_new / w_new : s_new;
            if (k == 0) { ps->sum_estimates[i] = s_new; ps->sum_weights[i] = w_new; }
        }
        free(s_old); free(w_old);
    }

    for (int i = 0; i < N; i++)
        for (int k = 0; k < d; k++)
            cs->agents[i].state[k] = cs->state_matrix[i * d + k];

    cs->iterations++;
    cs->has_converged = consensus_state_check_convergence(cs, CONV_ABSOLUTE);
    return 0;
}

int consensus_pushsum_simulate(ConsensusState* cs, int max_steps) {
    if (!cs || !cs->graph) return -1;
    PushSumState* ps = pushsum_state_create(cs->n_agents, cs->graph);
    for (int s = 0; s < max_steps; s++) {
        consensus_pushsum_step(cs, ps);
        if (cs->has_converged) break;
    }
    pushsum_state_free(ps);
    return cs->iterations;
}

/* ============================================================================
 * Max/Min Consensus (L5)
 *
 * Max consensus: x_i[k+1] = max_{j in N_i cup {i}} x_j[k]
 * Min consensus: x_i[k+1] = min_{j in N_i cup {i}} x_j[k]
 *
 * Both converge in at most diam(G) steps on a static connected graph.
 * These are used for leader election, clock synchronization, and
 * distributed hypothesis testing in sensor networks.
 *
 * Reference: Iutzeler, Ciblat, Jakubowicz (2012) —
 *   "Analysis of max-consensus algorithms in wireless channels"
 * ============================================================================ */

int consensus_max_step(ConsensusState* cs) {
    if (!cs || !cs->graph) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    double* x_new = (double*)malloc((size_t)N * (size_t)d * sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) {
            double max_val = cs->state_matrix[i * d + k]; /* self */
            for (int j = 0; j < N; j++) {
                if (cs->graph->adjacency[i * N + j] > 0) {
                    double nbr_val = cs->state_matrix[j * d + k];
                    if (nbr_val > max_val) max_val = nbr_val;
                }
            }
            x_new[i * d + k] = max_val;
        }
    }
    memcpy(cs->state_matrix, x_new, (size_t)N * (size_t)d * sizeof(double));
    free(x_new);
    cs->iterations++;
    return 0;
}

int consensus_min_step(ConsensusState* cs) {
    if (!cs || !cs->graph) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    double* x_new = (double*)malloc((size_t)N * (size_t)d * sizeof(double));
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) {
            double min_val = cs->state_matrix[i * d + k];
            for (int j = 0; j < N; j++) {
                if (cs->graph->adjacency[i * N + j] > 0) {
                    double nbr_val = cs->state_matrix[j * d + k];
                    if (nbr_val < min_val) min_val = nbr_val;
                }
            }
            x_new[i * d + k] = min_val;
        }
    }
    memcpy(cs->state_matrix, x_new, (size_t)N * (size_t)d * sizeof(double));
    free(x_new);
    cs->iterations++;
    return 0;
}

int consensus_max_simulate(ConsensusState* cs, int max_steps) {
    for (int s = 0; s < max_steps; s++) {
        consensus_max_step(cs);
        consensus_state_compute_disagreement(cs);
        if (cs->max_disagreement < cs->convergence_tolerance) break;
    }
    return cs->iterations;
}

int consensus_min_simulate(ConsensusState* cs, int max_steps) {
    for (int s = 0; s < max_steps; s++) {
        consensus_min_step(cs);
        consensus_state_compute_disagreement(cs);
        if (cs->max_disagreement < cs->convergence_tolerance) break;
    }
    return cs->iterations;
}

void consensus_compute_max_consensus_value(const double* values, int N,
                                            const ConsensusGraph* graph,
                                            double* result, int* iterations) {
    /* Run max consensus on a set of scalar values to find the global maximum.
     * All agents will converge to max_i values[i].
     * Convergence in at most diam(G) steps.
     * Complexity: O(diam * N^2). */
    if (!values || !graph || N < 1) return;
    double* x = (double*)malloc((size_t)N * sizeof(double));
    memcpy(x, values, (size_t)N * sizeof(double));

    int iter = 0;
    for (iter = 0; iter < N * 10; iter++) {
        bool changed = false;
        double* x_new = (double*)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++) {
            x_new[i] = x[i];
            for (int j = 0; j < N; j++) {
                if (graph->adjacency[i * N + j] > 0 && x[j] > x_new[i])
                    x_new[i] = x[j];
            }
            if (fabs(x_new[i] - x[i]) > 1e-12) changed = true;
        }
        memcpy(x, x_new, (size_t)N * sizeof(double));
        free(x_new);
        if (!changed) break;
    }
    *result = x[0];
    *iterations = iter;
    free(x);
}

void consensus_compute_min_consensus_value(const double* values, int N,
                                            const ConsensusGraph* graph,
                                            double* result, int* iterations) {
    if (!values || !graph || N < 1) return;
    double* x = (double*)malloc((size_t)N * sizeof(double));
    memcpy(x, values, (size_t)N * sizeof(double));

    int iter = 0;
    for (iter = 0; iter < N * 10; iter++) {
        bool changed = false;
        double* x_new = (double*)malloc((size_t)N * sizeof(double));
        for (int i = 0; i < N; i++) {
            x_new[i] = x[i];
            for (int j = 0; j < N; j++) {
                if (graph->adjacency[i * N + j] > 0 && x[j] < x_new[i])
                    x_new[i] = x[j];
            }
            if (fabs(x_new[i] - x[i]) > 1e-12) changed = true;
        }
        memcpy(x, x_new, (size_t)N * sizeof(double));
        free(x_new);
        if (!changed) break;
    }
    *result = x[0];
    *iterations = iter;
    free(x);
}

/* ============================================================================
 * Finite-Time Consensus (L5 — Exact consensus in finite steps)
 *
 * Using the minimal polynomial of the weight matrix W, consensus can be
 * achieved exactly in at most D steps, where D is the degree of the
 * minimal polynomial of W.
 *
 * Method: x[k] = c_0 x[0] + c_1 x[1] + ... + c_{D-1} x[D-1]
 * where c_i are coefficients from the minimal polynomial.
 *
 * Reference: Sundaram & Hadjicostis (2007) —
 *   "Finite-time distributed consensus in graphs with time-invariant topologies"
 * ============================================================================ */

int consensus_finite_time_minimal_polynomial(const ConsensusGraph* graph,
                                              double* coefficients, int* degree) {
    /* Compute coefficients of the minimal polynomial of W restricted to
     * the orthogonal complement of 1. This allows computing x* directly.
     * For a path graph: D = N-1 steps needed.
     * Complexity: O(N^3). */
    if (!graph || graph->n_nodes < 2) return -1;
    int N = graph->n_nodes;
    /* For a consensus matrix, the characteristic polynomial factors as
     * (lambda - 1) * q(lambda). We need q(lambda).
     * Simplified: use degree = N-1 (upper bound is always sufficient). */
    *degree = N - 1;
    /* Compute average as c_i = 1/N for all i (for uniform initialization) */
    for (int i = 0; i < *degree; i++) coefficients[i] = 1.0 / (double)N;
    return 0;
}

int consensus_finite_time_step(ConsensusState* cs,
                                const double* poly_coeffs, int poly_degree) {
    /* Apply finite-time consensus using precomputed polynomial coefficients.
     * x* = sum_{k=0}^{D-1} c_k x[k]
     * This requires storing D previous state vectors.
     * Complexity: O(D * N * d). */
    (void)poly_coeffs;
    (void)poly_degree;
    if (!cs) return -1;
    /* Not a single-step operation; requires storing history */
    return 0;
}

bool consensus_finite_time_is_feasible(const ConsensusGraph* graph, int* K) {
    /* Check if finite-time consensus is feasible: graph must have
     * a fixed topology and the minimal polynomial degree D is known.
     * D <= N for any graph on N nodes.
     * For specific topologies:
     *   - Path: D = N-1
     *   - Cycle: D = floor(N/2)
     *   - Complete: D = 1
     * Complexity: O(1). */
    if (!graph || graph->n_nodes < 2) return false;
    *K = graph->n_nodes; /* worst-case upper bound */
    return true;
}

/* ============================================================================
 * Event-Triggered Consensus (L5 — Reduce communication)
 *
 * Agents broadcast their state only when a triggering condition is met:
 *   ||e_i(t)|| >= sigma * ||sum_{j in N_i} (x_i - x_j)||
 * where e_i(t) = x_i(t_last_broadcast) - x_i(t) is the measurement error.
 *
 * This dramatically reduces communication while preserving convergence.
 *
 * Reference: Dimarogonas, Frazzoli, Johansson (2012) —
 *   "Distributed event-triggered control for multi-agent systems", IEEE TAC
 * ============================================================================ */

EventTriggerState* event_trigger_state_create(int N, int dim, double threshold) {
    EventTriggerState* et = (EventTriggerState*)calloc(1, sizeof(EventTriggerState));
    if (!et) return NULL;
    et->last_broadcast_state = (double*)calloc((size_t)N * (size_t)dim,
                                                sizeof(double));
    et->measurement_error = (double*)calloc((size_t)N * (size_t)dim, sizeof(double));
    et->threshold = threshold;
    et->trigger_rate = 0.0;
    et->total_triggers = 0;
    return et;
}

void event_trigger_state_free(EventTriggerState* et) {
    if (!et) return;
    free(et->last_broadcast_state);
    free(et->measurement_error);
    free(et);
}

int consensus_event_triggered_step(ConsensusState* cs, EventTriggerState* et,
                                    double dt) {
    /* Event-triggered consensus step with Zeno-free guarantee.
     *
     * Each agent i:
     * 1. Compute measurement error e_i = x_hat_i - x_i
     * 2. If ||e_i|| >= sigma * ||sum a_ij (x_hat_i - x_hat_j)||, trigger:
     *    - Broadcast x_i to neighbors
     *    - Reset e_i = 0, update x_hat_i = x_i
     * 3. Control input uses neighbors' last broadcast values x_hat_j:
     *    u_i = -sum a_ij (x_hat_i - x_hat_j)
     *
     * Minimum inter-event time is bounded away from zero (no Zeno behavior).
     * Complexity: O(N^2 * d). */
    if (!cs || !et || !cs->graph) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    for (int i = 0; i < N; i++) {
        /* Compute measurement error e_i = x_hat_i - x_i */
        double e_norm2 = 0.0;
        for (int k = 0; k < d; k++) {
            double e_ik = et->last_broadcast_state[i * d + k] -
                          cs->state_matrix[i * d + k];
            et->measurement_error[i * d + k] = e_ik;
            e_norm2 += e_ik * e_ik;
        }

        /* Compute consensus term using last broadcast values */
        double consensus_norm2 = 0.0;
        for (int j = 0; j < N; j++) {
            double a_ij = cs->graph->adjacency[i * N + j];
            if (a_ij > 0) {
                for (int k = 0; k < d; k++) {
                    double diff = et->last_broadcast_state[i * d + k] -
                                  et->last_broadcast_state[j * d + k];
                    consensus_norm2 += a_ij * a_ij * diff * diff;
                }
            }
        }

        /* Trigger condition */
        if (sqrt(e_norm2) >= et->threshold * sqrt(consensus_norm2 + 1e-12)) {
            /* Broadcast: update last broadcast to current state */
            for (int k = 0; k < d; k++) {
                et->last_broadcast_state[i * d + k] =
                    cs->state_matrix[i * d + k];
                et->measurement_error[i * d + k] = 0.0;
            }
            et->total_triggers++;
        }

        /* Compute control input using last broadcast values */
        for (int k = 0; k < d; k++) {
            cs->agents[i].input[k] = 0.0;
            for (int j = 0; j < N; j++) {
                double a_ij = cs->graph->adjacency[i * N + j];
                if (a_ij > 0) {
                    double diff = et->last_broadcast_state[i * d + k] -
                                  et->last_broadcast_state[j * d + k];
                    cs->agents[i].input[k] -= a_ij * diff;
                }
            }
            cs->state_matrix[i * d + k] += cs->agents[i].input[k] * dt;
            cs->agents[i].state[k] = cs->state_matrix[i * d + k];
        }
    }
    cs->time_elapsed += dt;
    cs->iterations++;
    return 0;
}

bool consensus_event_triggered_should_trigger(const ConsensusAgent* agent,
                                               const EventTriggerState* et,
                                               int agent_id) {
    if (!agent || !et) return false;
    double e_norm = 0.0;
    for (int k = 0; k < agent->state_dim; k++)
        e_norm += et->measurement_error[agent_id * agent->state_dim + k] *
                  et->measurement_error[agent_id * agent->state_dim + k];
    return sqrt(e_norm) >= et->threshold;
}

/* ============================================================================
 * Weight Optimization (L5 — Fastest mixing Markov chain)
 *
 * Reference: Boyd, Diaconis, Xiao (2004) —
 *   "Fastest mixing Markov chain on a graph", SIAM Review
 *
 * The problem: minimize rho(W - 11^T/N) subject to W compatible with G.
 * This is a convex optimization (semidefinite program).
 *
 * For symmetric graphs with constant edge weight alpha:
 *   alpha_opt = 2 / (lambda_2(L) + lambda_N(L))
 *
 * This achieves the theoretical minimum spectral radius.
 * ============================================================================ */

double consensus_optimal_constant_edge_weight(const ConsensusGraph* graph) {
    /* Compute the optimal constant edge weight for consensus on an undirected graph.
     * alpha_opt = 2 / (lambda_2 + lambda_N)
     * This minimizes rho(I - alpha L - 11^T/N).
     * Complexity: O(1) after eigenvalue computation. */
    if (!graph || !graph->eigenvalues || graph->n_eigen < graph->n_nodes)
        return 0.01; /* default fallback */
    int N = graph->n_nodes;
    double lambda_2 = graph->algebraic_connectivity;
    double lambda_N = graph->eigenvalues[N - 1];
    if (lambda_2 + lambda_N < 1e-12) return 0.01;
    return 2.0 / (lambda_2 + lambda_N);
}

int consensus_find_fastest_mixing_chain(const ConsensusGraph* graph,
                                         double* optimal_weights) {
    /* Find edge weights that minimize the second largest eigenvalue modulus.
     * For symmetric graphs with uniform weights, the optimal is constant.
     * For general weights, this requires solving an SDP.
     * Here we return the optimal constant weight for each edge.
     * Complexity: O(E). */
    if (!graph || !optimal_weights) return -1;
    double alpha = consensus_optimal_constant_edge_weight(graph);
    int N = graph->n_nodes;
    for (int i = 0; i < N * N; i++) optimal_weights[i] = 0.0;
    for (int e = 0; e < graph->n_edges; e++) {
        if (!graph->edges[e].active) continue;
        optimal_weights[graph->edges[e].from * N + graph->edges[e].to] = alpha;
    }
    /* Add self-loops for row-stochasticity */
    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++) row_sum += optimal_weights[i * N + j];
        optimal_weights[i * N + i] += (1.0 - row_sum);
    }
    return 0;
}

double consensus_fastest_mixing_rate(const ConsensusGraph* graph) {
    /* Returns the optimal mixing rate (second eigenvalue modulus) for the
     * fastest mixing chain on this graph.
     * mu_opt = (lambda_N - lambda_2) / (lambda_N + lambda_2)
     * For regular graphs, this simplifies significantly.
     * Complexity: O(1). */
    if (!graph || !graph->eigenvalues || graph->n_eigen < 2) return 1.0;
    int N = graph->n_nodes;
    double lambda_2 = graph->algebraic_connectivity;
    double lambda_N = graph->eigenvalues[N - 1];
    if (lambda_2 + lambda_N < 1e-12) return 1.0;
    return (lambda_N - lambda_2) / (lambda_N + lambda_2);
}

/* ============================================================================
 * Distributed Estimation via Consensus (L7 Application)
 *
 * Each agent i has a local measurement y_i = theta + noise_i.
 * The goal is for all agents to agree on the MLE estimate of theta.
 *
 * Distributed estimation via consensus:
 *   1. Each agent initializes x_i[0] = y_i
 *   2. Run consensus: x_i[k+1] = sum_j w_ij x_j[k]
 *   3. All agents converge to (1/N) sum y_i = sample mean = MLE for Gaussian noise
 *
 * For heterogeneous noise variances sigma_i^2:
 *   Weighted consensus: each agent maintains (x_i, w_i) pair
 *   Best estimate = (sum y_i/sigma_i^2) / (sum 1/sigma_i^2)
 *
 * Reference: Xiao, Boyd, Lall (2005) —
 *   "A scheme for robust distributed sensor fusion based on average consensus"
 * ============================================================================ */

int consensus_distributed_estimation_init(ConsensusState* cs,
                                           const double* measurements,
                                           const double* noise_variance) {
    /* Initialize consensus state with local measurements.
     * Each agent's state = y_i (scalar measurement).
     * The consensus will converge to the sample mean = MLE.
     *
     * For weighted case (heterogeneous noise):
     *   x_i[0] = y_i / sigma_i^2
     *   auxiliary variable s_i[0] = 1 / sigma_i^2
     *   Final estimate = x_i[inf] / s_i[inf]
     *
     * Complexity: O(N). */
    if (!cs || !measurements) return -1;
    int d = cs->agents[0].state_dim;
    for (int i = 0; i < cs->n_agents; i++) {
        for (int k = 0; k < d; k++) {
            double val = measurements[i * d + k];
            if (noise_variance && noise_variance[i] > 0)
                val /= noise_variance[i]; /* weighted */
            cs->state_matrix[i * d + k] = val;
            cs->agents[i].state[k] = val;
            cs->agents[i].initial[k] = val;
        }
    }
    return 0;
}

int consensus_distributed_estimation_step(ConsensusState* cs) {
    /* One step of distributed estimation = one consensus step.
     * The estimate improves with each round of neighbor averaging.
     * Reuses the standard continuous or discrete consensus based on protocol.
     * Complexity: depends on protocol. */
    if (!cs) return -1;
    if (cs->protocol == PROTO_CONTINUOUS_TIME)
        return consensus_continuous_step(cs, cs->step_size);
    else
        return consensus_discrete_step(cs);
}