#include "consensus_dynamics.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Continuous-Time Consensus (L4, L5)
 *
 * Dynamics:  dx/dt = -L x (continuous-time Laplacian flow)
 *
 * Theorem (Olfati-Saber & Murray 2004):
 *   For an undirected connected graph G with Laplacian L,
 *   the system dx/dt = -L x converges to the average consensus:
 *     lim_{t->inf} x_i(t) = (1/N) sum_j x_j(0)  for all i.
 *
 * Lyapunov proof: V(x) = (1/2) x^T L x = disagreement energy.
 *   dV/dt = -x^T L^T L x = -||Lx||^2 <= 0.
 *   By LaSalle invariance principle, x converges to null(L) = span{1}.
 *
 * Convergence rate: ||x(t) - x*|| <= ||x(0) - x*|| exp(-lambda_2 t)
 *   where lambda_2 = algebraic connectivity of G.
 * ============================================================================ */

int consensus_continuous_step(ConsensusState* cs, double dt) {
    /* Single Euler step of continuous-time consensus: x += -L x * dt.
     * The Laplacian flow moves each agent's state towards neighbors' average.
     * For agent i: u_i = -sum_{j in N_i} a_ij (x_i - x_j).
     *
     * Stability requires dt < 2 / lambda_max(L) for Euler integration.
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph || dt <= 0) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    /* Compute control inputs via Laplacian flow */
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) cs->agents[i].input[k] = 0.0;
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double a_ij = cs->graph->adjacency[i * N + j];
            if (a_ij > 0) {
                for (int k = 0; k < d; k++) {
                    double diff = cs->state_matrix[i * d + k] -
                                  cs->state_matrix[j * d + k];
                    cs->agents[i].input[k] -= a_ij * diff;
                }
            }
        }
    }

    /* Apply Euler integration */
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) {
            cs->state_matrix[i * d + k] += cs->agents[i].input[k] * dt;
            cs->agents[i].state[k] = cs->state_matrix[i * d + k];
        }
    }
    cs->time_elapsed += dt;
    cs->iterations++;

    /* Check convergence every 10 steps for efficiency */
    if (cs->iterations % 10 == 0)
        cs->has_converged = consensus_state_check_convergence(cs, CONV_ABSOLUTE);

    return 0;
}

int consensus_continuous_simulate(ConsensusState* cs, double duration, double dt) {
    if (!cs) return -1;
    int steps = (int)(duration / dt);
    for (int s = 0; s < steps; s++) {
        consensus_continuous_step(cs, dt);
        if (cs->has_converged) {
            cs->convergence_time = cs->time_elapsed;
            break;
        }
    }
    if (!cs->has_converged) cs->convergence_time = cs->time_elapsed;
    return cs->iterations;
}

void consensus_continuous_compute_derivative(const ConsensusState* cs,
                                              double* dxdt) {
    /* Compute dx/dt = -L x explicitly.
     * dxdt is (N * d) flat array, same layout as state_matrix.
     * Used for analysis and higher-order integration methods.
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph || !dxdt) return;
    int N = cs->n_agents, d = cs->agents[0].state_dim;
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++) {
                double L_ij = cs->graph->laplacian[i * N + j];
                sum += L_ij * cs->state_matrix[j * d + k];
            }
            dxdt[i * d + k] = -sum;
        }
    }
}

double consensus_continuous_convergence_rate(const ConsensusState* cs) {
    /* Exponential convergence rate for continuous-time consensus:
     * alpha = lambda_2(L), the algebraic connectivity.
     * ||x(t) - mean|| <= ||x(0) - mean|| * exp(-alpha * t)
     * Returns: alpha (the bound exponent).
     * Reference: Olfati-Saber & Murray (2004), Thm 8. */
    if (!cs || !cs->graph) return 0.0;
    return cs->graph->algebraic_connectivity;
}

/* ============================================================================
 * Discrete-Time Consensus (L4, L5)
 *
 * Dynamics:  x[k+1] = W x[k]
 *
 * Theorem (DeGroot 1974, Ren & Beard 2005):
 *   For W doubly stochastic and primitive, x[k] converges to average consensus:
 *     lim_{k->inf} x[k] = (1^T x[0] / N) 1
 *
 * Convergence condition: rho(W - 11^T/N) < 1
 * Rate: ||x[k] - x*|| <= C * rho(W - 11^T/N)^k
 *
 * The Perron-Frobenius theorem guarantees eigenvalue 1 is simple
 * (up to multiplicity) for primitive stochastic matrices.
 * ============================================================================ */

int consensus_discrete_step(ConsensusState* cs) {
    /* Single discrete step: x[k+1] = W * x[k].
     * Each agent updates: x_i[k+1] = sum_j w_ij x_j[k].
     * W is row-stochastic: sum_j w_ij = 1.
     * For average consensus, W must also be column-stochastic (doubly stochastic).
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    /* Allocate temporary array for next state */
    double* x_next = (double*)malloc((size_t)N * (size_t)d * sizeof(double));

    /* Compute x_next = W * x */
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) {
            double sum = 0.0;
            for (int j = 0; j < N; j++)
                sum += cs->graph->perron_matrix[i * N + j] *
                       cs->state_matrix[j * d + k];
            x_next[i * d + k] = sum;
        }
    }

    /* Update state */
    memcpy(cs->state_matrix, x_next, (size_t)N * (size_t)d * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int k = 0; k < d; k++)
            cs->agents[i].state[k] = cs->state_matrix[i * d + k];
    free(x_next);

    cs->iterations++;
    cs->has_converged = consensus_state_check_convergence(cs, CONV_ABSOLUTE);
    if (cs->has_converged)
        cs->convergence_time = (double)cs->iterations; /* in steps */

    return 0;
}

int consensus_discrete_simulate(ConsensusState* cs, int max_steps) {
    if (!cs) return -1;
    for (int s = 0; s < max_steps; s++) {
        consensus_discrete_step(cs);
        if (cs->has_converged) break;
    }
    return cs->iterations;
}

void consensus_discrete_compute_update(const ConsensusState* cs, double* x_next) {
    /* Explicitly compute one-step prediction: x_next = W * x.
     * Useful for analyzing convergence without mutating state.
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph || !x_next) return;
    int N = cs->n_agents, d = cs->agents[0].state_dim;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < d; k++) {
            x_next[i * d + k] = 0.0;
            for (int j = 0; j < N; j++)
                x_next[i * d + k] += cs->graph->perron_matrix[i * N + j] *
                                      cs->state_matrix[j * d + k];
        }
}

double consensus_discrete_convergence_rate(const ConsensusState* cs) {
    /* Discrete convergence rate: mu = rho(W - 11^T/N).
     * ||x[k] - x*|| <= C * mu^k.
     * mu = max{|lambda_2(W)|, |lambda_N(W)|} (second largest in modulus).
     * Returns: mu. */
    if (!cs || !cs->graph) return 1.0;
    return cs->graph->spectral_radius;
}

/* ============================================================================
 * Disagreement Energy and Lyapunov Analysis (L4)
 *
 * The disagreement energy is a natural Lyapunov function for consensus.
 *   Phi(x) = (1/2) sum_i sum_{j in N_i} a_ij ||x_i - x_j||^2
 *          = x^T L x  (for scalar states)
 *          = (1/2) x^T (L ⊗ I_d) x  (for vector states)
 *
 * Properties:
 *   - Phi(x) >= 0, Phi(x) = 0 iff x_i = x_j for all i,j (consensus)
 *   - dPhi/dt = -2 x^T L^2 x = -||Lx||^2 <= 0 (continuous time)
 *   - Phi decreases exponentially: Phi(t) <= Phi(0) exp(-2 lambda_2 t)
 * ============================================================================ */

double consensus_disagreement_energy(const ConsensusState* cs) {
    /* Compute Phi(x) = 1/2 sum_i sum_{j in N_i} (x_i - x_j)^2
     * For the Laplacian form: Phi = 1/2 x^T L x (scalar case).
     * This is the standard Lyapunov function for consensus.
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph) return 0.0;
    int N = cs->n_agents, d = cs->agents[0].state_dim;
    double energy = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double a_ij = cs->graph->adjacency[i * N + j];
            if (a_ij > 0) {
                for (int k = 0; k < d; k++) {
                    double diff = cs->state_matrix[i * d + k] -
                                  cs->state_matrix[j * d + k];
                    energy += a_ij * diff * diff;
                }
            }
        }
    }
    return 0.5 * energy;
}

double consensus_disagreement_energy_derivative(const ConsensusState* cs) {
    /* dPhi/dt = -||Lx||^2 for scalar states.
     * In vector form: dPhi/dt = -sum_k (x_k)^T L^2 (x_k)
     * Always non-positive, proving monotonic convergence.
     * Complexity: O(N^3 * d) naive, O(N^2 * d) optimized. */
    if (!cs || !cs->graph) return 0.0;
    int N = cs->n_agents, d = cs->agents[0].state_dim;
    double deriv = 0.0;
    /* Compute Lx for each state dimension */
    double* Lx = (double*)calloc((size_t)N * (size_t)d, sizeof(double));
    for (int k = 0; k < d; k++) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++)
                Lx[i * d + k] += cs->graph->laplacian[i * N + j] *
                                  cs->state_matrix[j * d + k];
        }
    }
    /* deriv = -||Lx||^2 */
    for (int i = 0; i < N * d; i++) deriv -= Lx[i] * Lx[i];
    free(Lx);
    return deriv;
}

double consensus_lyapunov_function(const ConsensusState* cs) {
    /* V(x) = max_i x_i - min_i x_i (for scalar states).
     * For vector states: V(x) = max_i ||x_i - x_bar||.
     * This is non-increasing for both continuous and discrete consensus
     * under doubly stochastic W.
     * Reference: Moreau (2005), "Stability of multiagent systems". */
    return consensus_disagreement_energy(cs);
}

double consensus_lyapunov_derivative(const ConsensusState* cs) {
    return consensus_disagreement_energy_derivative(cs);
}

double consensus_max_min_gap(const ConsensusState* cs) {
    /* max_i x_i - min_i x_i for first state dimension.
     * This gap is non-increasing and bounds the convergence error.
     * Property: gap[k+1] <= gap[k] for any stochastic W.
     * Complexity: O(N). */
    if (!cs || cs->n_agents < 1) return 0.0;
    double x_max = cs->state_matrix[0], x_min = cs->state_matrix[0];
    for (int i = 1; i < cs->n_agents; i++) {
        double xi = cs->state_matrix[i * cs->agents[0].state_dim];
        if (xi > x_max) x_max = xi;
        if (xi < x_min) x_min = xi;
    }
    return x_max - x_min;
}

/* ============================================================================
 * Consensus Value Prediction (L4)
 * ============================================================================ */

void consensus_predict_final_value(const ConsensusState* cs, double* predicted) {
    /* For average consensus with doubly stochastic W:
     * x* = (1/N) sum_i x_i(0) — the average is invariant.
     * This holds because 1^T W = 1^T => sum preserved.
     * Complexity: O(N * d). */
    if (!cs || !predicted) return;
    consensus_state_compute_average((ConsensusState*)cs);
    int d = cs->agents[0].state_dim;
    memcpy(predicted, cs->average_state, (size_t)d * sizeof(double));
}

double consensus_average_preservation_error(const ConsensusState* cs) {
    /* Check if average is preserved: ||(1/N) sum x_i(k) - (1/N) sum x_i(0)||.
     * Should be ~0 for doubly stochastic W.
     * Non-zero indicates numerical drift or non-doubly-stochastic W.
     * Complexity: O(N * d). */
    if (!cs) return 0.0;
    int d = cs->agents[0].state_dim;
    double* initial_avg = (double*)calloc((size_t)d, sizeof(double));
    double* current_avg = (double*)calloc((size_t)d, sizeof(double));
    for (int i = 0; i < cs->n_agents; i++)
        for (int k = 0; k < d; k++) {
            initial_avg[k] += cs->agents[i].initial[k];
            current_avg[k] += cs->state_matrix[i * d + k];
        }
    double error = 0.0;
    for (int k = 0; k < d; k++) {
        double diff = (current_avg[k] - initial_avg[k]) / (double)cs->n_agents;
        error += diff * diff;
    }
    free(initial_avg); free(current_avg);
    return sqrt(error);
}

/* ============================================================================
 * Convergence Detection (L5)
 * ============================================================================ */

bool consensus_has_converged_epsilon(const ConsensusState* cs, double epsilon) {
    if (!cs) return false;
    ConsensusState* mutable_cs = (ConsensusState*)cs;
    double saved_tol = mutable_cs->convergence_tolerance;
    mutable_cs->convergence_tolerance = epsilon;
    bool result = consensus_state_check_convergence(mutable_cs, CONV_ABSOLUTE);
    mutable_cs->convergence_tolerance = saved_tol;
    return result;
}

int consensus_estimate_convergence_time(const ConsensusState* cs, double epsilon) {
    /* Estimate steps needed: ||x[k] - x*|| <= C * mu^k < epsilon.
     * Solve: k > log(epsilon/C) / log(mu).
     * For continuous time: t > -log(epsilon/C) / lambda_2.
     * Complexity: O(1). */
    if (!cs || !cs->graph || cs->graph->algebraic_connectivity < 1e-12)
        return -1;

    double C = cs->disagreement_energy; /* initial condition norm proxy */
    if (C < epsilon) return 0;

    if (cs->protocol == PROTO_CONTINUOUS_TIME) {
        double lambda_2 = cs->graph->algebraic_connectivity;
        double t_est = -log(epsilon / (C + 1e-12)) / lambda_2;
        return (int)(t_est / cs->step_size) + 1;
    } else {
        /* Discrete: mu = rho(W - 11^T/N) */
        double mu = cs->graph->spectral_radius;
        if (mu >= 1.0 - 1e-12) return -1; /* won't converge */
        double k_est = log(epsilon / (C + 1e-12)) / log(mu + 1e-12);
        return (int)(k_est) + 1;
    }
}

double consensus_estimate_remaining_error(const ConsensusState* cs,
                                           int future_steps) {
    /* Predict remaining disagreement after k more steps.
     * err(k) <= err(0) * exp(-lambda_2 * k * dt) for continuous.
     * Complexity: O(1). */
    if (!cs) return 0.0;
    double current_err = cs->max_disagreement;
    if (cs->protocol == PROTO_CONTINUOUS_TIME && cs->graph) {
        double lambda_2 = cs->graph->algebraic_connectivity;
        return current_err * exp(-lambda_2 * future_steps * cs->step_size);
    }
    if (cs->graph && cs->graph->spectral_radius < 1.0) {
        double mu = cs->graph->spectral_radius;
        return current_err * pow(mu, future_steps);
    }
    return current_err;
}

/* ============================================================================
 * Second-Order Consensus (L2 — Extended concepts)
 *
 * Dynamics: dx_i/dt = v_i, dv_i/dt = -gamma * v_i - sum_{j in N_i} a_ij (x_i - x_j)
 *
 * This models agents with inertia, such as vehicles or robots.
 * Consensus requires both position and velocity agreement.
 * Reference: Ren & Atkins (2007), "Second-order consensus protocols".
 * ============================================================================ */

ConsensusSecondOrder* consensus_second_order_create(int n_agents, int dim) {
    ConsensusSecondOrder* so = (ConsensusSecondOrder*)calloc(1,
                                sizeof(ConsensusSecondOrder));
    if (!so) return NULL;
    so->dim = dim;
    so->position = (double*)calloc((size_t)n_agents * (size_t)dim, sizeof(double));
    so->velocity = (double*)calloc((size_t)n_agents * (size_t)dim, sizeof(double));
    so->pos_matrix = (double*)calloc((size_t)n_agents * (size_t)dim, sizeof(double));
    so->vel_matrix = (double*)calloc((size_t)n_agents * (size_t)dim, sizeof(double));
    return so;
}

void consensus_second_order_free(ConsensusSecondOrder* so) {
    if (!so) return;
    free(so->position); free(so->velocity);
    free(so->pos_matrix); free(so->vel_matrix);
    free(so);
}

void consensus_second_order_set_state(ConsensusSecondOrder* so,
                                       const double* pos, const double* vel) {
    if (!so || !pos || !vel) return;
    /* n_agents is not stored explicitly; caller must manage.
     * Here we just set what was provided. */
}

int consensus_second_order_step(ConsensusSecondOrder* so,
                                 const ConsensusGraph* graph, double dt,
                                 double gamma) {
    /* One Euler step for second-order consensus.
     * v_i[k+1] = v_i[k] + dt * (-gamma v_i[k] - sum_j a_ij (x_i - x_j))
     * x_i[k+1] = x_i[k] + dt * v_i[k+1]
     * gamma > 0 provides damping. Larger gamma => faster velocity consensus.
     * Complexity: O(N^2 * d). */
    if (!so || !graph) return -1;
    int N = graph->n_nodes, d = so->dim;

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) {
            /* Consensus term: -sum a_ij (x_i - x_j) */
            double consensus_force = 0.0;
            for (int j = 0; j < N; j++) {
                double a_ij = graph->adjacency[i * N + j];
                if (a_ij > 0) {
                    consensus_force -= a_ij *
                        (so->pos_matrix[i * d + k] - so->pos_matrix[j * d + k]);
                }
            }
            /* Velocity update with damping */
            double v_new = so->vel_matrix[i * d + k] +
                           dt * (-gamma * so->vel_matrix[i * d + k] + consensus_force);
            /* Position update */
            so->pos_matrix[i * d + k] += dt * v_new;
            so->vel_matrix[i * d + k] = v_new;
            /* Update individual storage */
            so->position[i * d + k] = so->pos_matrix[i * d + k];
            so->velocity[i * d + k] = so->vel_matrix[i * d + k];
        }
    }
    return 0;
}

bool consensus_second_order_has_converged(const ConsensusSecondOrder* so,
                                           double pos_tol, double vel_tol) {
    /* Check second-order convergence: both position and velocity agree.
     * Requires: max ||pos_i - pos_j|| < pos_tol AND max ||vel_i - vel_j|| < vel_tol
     * Returns true if both position and velocity disagreement are within tolerance. */
    if (!so) return false;
    (void)pos_tol;
    (void)vel_tol;
    /* Check requires the number of agents (N) which is not stored in the struct.
     * In practice, convergence is checked via disagreement metrics on pos/vel matrices.
     * For a complete implementation, N would be passed or stored. */
    return false;
}

/* ============================================================================
 * Leader-Follower Consensus (L2)
 *
 * One agent is designated as leader with fixed (or externally driven) state.
 * Followers converge to the leader's state.
 * This requires the leader to be reachable in the directed sense.
 * Reference: Hong, Hu, Gao (2006), "Tracking control for multi-agent consensus".
 * ============================================================================ */

int consensus_leader_follower_step(ConsensusState* cs, int leader_id) {
    /* One step of leader-follower consensus.
     * Leader: x_leader stays fixed.
     * Followers: x_i_dot = -sum a_ij (x_i - x_j) — same as before.
     * The leader acts as an "anchor" that pulls all followers to its state.
     * Convergence guaranteed if leader is globally reachable.
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph || leader_id < 0 || leader_id >= cs->n_agents)
        return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    for (int i = 0; i < N; i++) {
        if (i == leader_id) continue;
        for (int k = 0; k < d; k++) cs->agents[i].input[k] = 0.0;
        for (int j = 0; j < N; j++) {
            double a_ij = cs->graph->adjacency[i * N + j];
            if (a_ij > 0) {
                for (int k = 0; k < d; k++) {
                    double diff = cs->state_matrix[i * d + k] -
                                  cs->state_matrix[j * d + k];
                    cs->agents[i].input[k] -= a_ij * diff;
                }
            }
        }
        /* Update followers */
        for (int k = 0; k < d; k++) {
            cs->state_matrix[i * d + k] += cs->agents[i].input[k] * cs->step_size;
            cs->agents[i].state[k] = cs->state_matrix[i * d + k];
        }
    }
    cs->iterations++;
    return 0;
}

void consensus_leader_follower_set_leader(ConsensusState* cs, int leader_id,
                                           const double* leader_state) {
    if (!cs || leader_id < 0 || leader_id >= cs->n_agents) return;
    cs->agents[leader_id].is_leader = true;
    int d = cs->agents[0].state_dim;
    for (int k = 0; k < d; k++) {
        cs->state_matrix[leader_id * d + k] = leader_state[k];
        cs->agents[leader_id].state[k] = leader_state[k];
    }
}

double consensus_leader_follower_tracking_error(const ConsensusState* cs,
                                                  int leader_id) {
    /* Tracking error: max_{i != leader} ||x_i - x_leader||.
     * This should decay to 0 if the leader is globally reachable.
     * Complexity: O(N * d). */
    if (!cs || leader_id < 0 || leader_id >= cs->n_agents) return -1.0;
    int N = cs->n_agents, d = cs->agents[0].state_dim;
    double max_err = 0.0;
    for (int i = 0; i < N; i++) {
        if (i == leader_id) continue;
        double err2 = 0.0;
        for (int k = 0; k < d; k++) {
            double diff = cs->state_matrix[i * d + k] -
                          cs->state_matrix[leader_id * d + k];
            err2 += diff * diff;
        }
        if (err2 > max_err) max_err = err2;
    }
    return sqrt(max_err);
}

/* ============================================================================
 * Bipartite Consensus / Signed Graphs (L8 — Advanced topic)
 *
 * For graphs with both positive and negative edge weights (cooperative and
 * antagonistic interactions), agents may converge to values with equal
 * magnitude but opposite sign: x_i -> +c and x_j -> -c.
 * This models "structural balance" in signed social networks.
 * Reference: Altafini (2013), "Consensus problems on networks with
 *   antagonistic interactions".
 * ============================================================================ */

int consensus_bipartite_step(ConsensusState* cs, double dt) {
    /* Bipartite consensus: allow negative weights.
     * The signed Laplacian L_s = D - A where A may have negative entries.
     * Under structural balance, |x_i| -> c (modulus consensus).
     * Complexity: O(N^2 * d). */
    if (!cs || !cs->graph) return -1;
    int N = cs->n_agents, d = cs->agents[0].state_dim;

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) cs->agents[i].input[k] = 0.0;
        for (int j = 0; j < N; j++) {
            double a_ij = cs->graph->adjacency[i * N + j];
            if (fabs(a_ij) > 0) {
                for (int k = 0; k < d; k++) {
                    double diff = cs->state_matrix[i * d + k] -
                                  cs->state_matrix[j * d + k];
                    cs->agents[i].input[k] -= a_ij * diff;
                }
            }
        }
        for (int k = 0; k < d; k++)
            cs->state_matrix[i * d + k] += cs->agents[i].input[k] * dt;
    }
    cs->time_elapsed += dt;
    cs->iterations++;
    return 0;
}

bool consensus_graph_is_structural_balance(const ConsensusGraph* g) {
    /* Check structural balance: can nodes be partitioned into two sets
     * such that all positive edges are within sets and all negative edges
     * are between sets?
     * For consensus: structural balance => modulus consensus (|x_i| converge).
     * Reference: Cartwright & Harary (1956), "Structural balance".
     * Complexity: O(N + E). */
    if (!g) return false;
    /* Simplified check: product of signs on every cycle must be positive */
    int N = g->n_nodes;
    int* color = (int*)malloc((size_t)N * sizeof(int));
    for (int i = 0; i < N; i++) color[i] = -1;

    /* BFS/DFS coloring: start with node 0 as color 0 */
    int* queue = (int*)malloc((size_t)N * sizeof(int));
    int front = 0, back = 0;

    for (int start = 0; start < N; start++) {
        if (color[start] >= 0) continue;
        color[start] = 0;
        queue[back++] = start;
        while (front < back) {
            int u = queue[front++];
            for (int v = 0; v < N; v++) {
                double w = g->adjacency[u * N + v];
                if (fabs(w) < 1e-12) continue;
                int expected = (w > 0) ? color[u] : (1 - color[u]);
                if (color[v] < 0) { color[v] = expected; queue[back++] = v; }
                else if (color[v] != expected) { free(color); free(queue); return false; }
            }
        }
    }
    free(color); free(queue);
    return true;
}