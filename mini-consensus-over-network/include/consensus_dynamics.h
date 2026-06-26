#ifndef CONSENSUS_DYNAMICS_H
#define CONSENSUS_DYNAMICS_H
#include "consensus_types.h"
/* ============================================================================
 * Consensus Dynamics (L4 - Fundamental Laws, L5 - Algorithms)
 *
 * Continuous-time consensus:  dx/dt = -L x
 *   where L is the graph Laplacian. Convergence to average:
 *     lim_{t->inf} x(t) = (1/N) 1 1^T x(0)  iff G connected
 *   Convergence rate: ||x(t) - x*|| <= exp(-lambda_2 t) ||x(0) - x*||
 *
 * Discrete-time consensus:  x[k+1] = W x[k]
 *   where W is a stochastic matrix compatible with G.
 *   Convergence to average iff W is doubly stochastic and primitive.
 *   Convergence rate: ||x[k] - x*|| <= rho(W - 11^T/N)^k
 *
 * Disagreement energy: Phi(x) = (1/2) x^T L x = (1/4) sum a_ij (x_i - x_j)^2
 *   dPhi/dt = -x^T L^2 x <= -lambda_2 Phi (continuous-time)
 *   Providing a Lyapunov function proof of convergence.
 *
 * Key Theorems (L4):
 *   T1: Consensus Theorem (Olfati-Saber 2004)
 *   T2: Convergence rate = O(exp(-lambda_2 t)) (undirected, continuous)
 *   T3: Discrete average consensus: W doubly stochastic + primitive => convergence
 *   T4: Lyapunov stability: V(x) = max_i x_i - min_i x_i non-increasing
 *   T5: Invariance principle: limit set = agreement subspace
 *   T6: Switching topology: consensus if union graph jointly connected
 * ============================================================================ */

/* Continuous-time consensus (L5 Algorithm) */
int consensus_continuous_step(ConsensusState* cs, double dt);
int consensus_continuous_simulate(ConsensusState* cs, double duration, double dt);
void consensus_continuous_compute_derivative(const ConsensusState* cs,
                                              double* dxdt);
double consensus_continuous_convergence_rate(const ConsensusState* cs);

/* Discrete-time consensus (L5 Algorithm) */
int consensus_discrete_step(ConsensusState* cs);
int consensus_discrete_simulate(ConsensusState* cs, int max_steps);
void consensus_discrete_compute_update(const ConsensusState* cs, double* x_next);
double consensus_discrete_convergence_rate(const ConsensusState* cs);

/* Disagreement energy and Lyapunov analysis (L4) */
double consensus_disagreement_energy(const ConsensusState* cs);
double consensus_disagreement_energy_derivative(const ConsensusState* cs);
double consensus_lyapunov_function(const ConsensusState* cs);
double consensus_lyapunov_derivative(const ConsensusState* cs);
double consensus_max_min_gap(const ConsensusState* cs);

/* Consensus value prediction (L4) */
void consensus_predict_final_value(const ConsensusState* cs, double* predicted);
double consensus_average_preservation_error(const ConsensusState* cs);

/* Convergence detection (L5) */
bool consensus_has_converged_epsilon(const ConsensusState* cs, double epsilon);
int consensus_estimate_convergence_time(const ConsensusState* cs, double epsilon);
double consensus_estimate_remaining_error(const ConsensusState* cs, int future_steps);

/* High-order consensus (L2 extended concept) */
typedef struct {
    double* position;
    double* velocity;
    double* acceleration;
    int dim;
    double* pos_matrix;  /* N x dim for all agents */
    double* vel_matrix;
} ConsensusSecondOrder;

ConsensusSecondOrder* consensus_second_order_create(int n_agents, int dim);
void consensus_second_order_free(ConsensusSecondOrder* so);
void consensus_second_order_set_state(ConsensusSecondOrder* so,
                                       const double* pos, const double* vel);
int consensus_second_order_step(ConsensusSecondOrder* so,
                                 const ConsensusGraph* graph, double dt, double gamma);
bool consensus_second_order_has_converged(const ConsensusSecondOrder* so,
                                           double pos_tol, double vel_tol);

/* Leader-follower consensus (L2) */
int consensus_leader_follower_step(ConsensusState* cs, int leader_id);
void consensus_leader_follower_set_leader(ConsensusState* cs, int leader_id,
                                           const double* leader_state);
double consensus_leader_follower_tracking_error(const ConsensusState* cs,
                                                  int leader_id);

/* Bipartite consensus / signed graphs (L8) */
int consensus_bipartite_step(ConsensusState* cs, double dt);
bool consensus_graph_is_structural_balance(const ConsensusGraph* g);

#endif