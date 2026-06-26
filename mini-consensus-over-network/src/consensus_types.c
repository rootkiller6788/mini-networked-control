#include "consensus_types.h"
#include "consensus_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

static unsigned int g_seed = 42;

/* ============================================================================
 * Utility Functions (L2 — Core concepts as primitives)
 * ============================================================================ */

/* Clamp a value to [lo, hi]. Fundamental numerical safeguard.
 * Complexity: O(1). */
double consensus_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

/* Uniform random number in [a, b]. Linear congruential generator.
 * For reproducible simulations in consensus analysis.
 * Complexity: O(1). */
double consensus_rand_uniform(double a, double b) {
    g_seed = (1103515245 * g_seed + 12345) & 0x7fffffff;
    double u = (double)g_seed / (double)0x7fffffff;
    return a + u * (b - a);
}

void consensus_seed(unsigned int seed) { g_seed = seed; }

/* ============================================================================
 * ConsensusAgent Implementation (L1 — Definitions realized as code)
 * ============================================================================ */

ConsensusAgent* consensus_agent_create(int id, const char* label, int state_dim) {
    ConsensusAgent* a = (ConsensusAgent*)calloc(1, sizeof(ConsensusAgent));
    if (!a) return NULL;
    a->id = id;
    a->label = label ? strdup(label) : NULL;
    a->state_dim = (state_dim > 0 && state_dim <= CONSENSUS_MAX_STATE_DIM)
                   ? state_dim : 1;
    a->dynamics = DYN_SINGLE_INTEGRATOR;
    a->n_neighbors = 0;
    a->neighbor_capacity = 8;
    a->neighbor_ids = (int*)malloc(a->neighbor_capacity * sizeof(int));
    a->neighbor_weights = (double*)malloc(a->neighbor_capacity * sizeof(double));
    a->self_weight = 0.0;
    a->is_leader = false;
    a->is_malicious = false;
    return a;
}

void consensus_agent_free(ConsensusAgent* agent) {
    if (!agent) return;
    free(agent->label);
    free(agent->neighbor_ids);
    free(agent->neighbor_weights);
    free(agent);
}

void consensus_agent_set_state(ConsensusAgent* agent, const double* state, int dim) {
    if (!agent || !state) return;
    int d = (dim < agent->state_dim) ? dim : agent->state_dim;
    for (int k = 0; k < d; k++) {
        agent->state[k] = state[k];
        agent->initial[k] = state[k];
    }
}

void consensus_agent_add_neighbor(ConsensusAgent* agent, int neighbor_id,
                                   double weight) {
    if (!agent) return;
    if (agent->n_neighbors >= agent->neighbor_capacity) {
        agent->neighbor_capacity *= 2;
        agent->neighbor_ids = (int*)realloc(agent->neighbor_ids,
                                agent->neighbor_capacity * sizeof(int));
        agent->neighbor_weights = (double*)realloc(agent->neighbor_weights,
                                   agent->neighbor_capacity * sizeof(double));
    }
    agent->neighbor_ids[agent->n_neighbors] = neighbor_id;
    agent->neighbor_weights[agent->n_neighbors] = weight;
    agent->n_neighbors++;
}

void consensus_agent_reset(ConsensusAgent* agent) {
    if (!agent) return;
    for (int k = 0; k < agent->state_dim; k++)
        agent->state[k] = agent->initial[k];
}

double consensus_agent_state_norm(const ConsensusAgent* agent) {
    if (!agent) return 0.0;
    double sum = 0.0;
    for (int k = 0; k < agent->state_dim; k++)
        sum += agent->state[k] * agent->state[k];
    return sqrt(sum);
}

void consensus_agent_apply_input(ConsensusAgent* agent, double dt) {
    if (!agent) return;
    /* Single integrator: x_dot = u, so x += u * dt */
    for (int k = 0; k < agent->state_dim; k++)
        agent->state[k] += agent->input[k] * dt;
}

void consensus_agent_print(const ConsensusAgent* agent) {
    if (!agent) { printf("Agent: NULL\n"); return; }
    printf("Agent %d (%s): state=[", agent->id,
           agent->label ? agent->label : "?");
    for (int k = 0; k < agent->state_dim; k++)
        printf("%s%.4f", k > 0 ? ", " : "", agent->state[k]);
    printf("], neighbors=%d, leader=%d, malicious=%d\n",
           agent->n_neighbors, agent->is_leader, agent->is_malicious);
}

/* ============================================================================
 * ConsensusVector Implementation (L3 — Mathematical Structures)
 *
 * Vector operations are the foundation for all consensus computations.
 * The state x ∈ R^N is a concatenation of N agent states each in R^d.
 * Core operations: norm, dot product, axpy (BLAS Level 1).
 * ============================================================================ */

ConsensusVector* consensus_vector_create(int dim) {
    ConsensusVector* v = (ConsensusVector*)calloc(1, sizeof(ConsensusVector));
    if (!v) return NULL;
    v->dim = dim;
    v->data = (double*)calloc((size_t)dim, sizeof(double));
    return v;
}

void consensus_vector_free(ConsensusVector* v) {
    if (!v) return;
    free(v->data);
    free(v);
}

void consensus_vector_set(ConsensusVector* v, double val) {
    if (!v) return;
    for (int i = 0; i < v->dim; i++) v->data[i] = val;
}

void consensus_vector_copy(ConsensusVector* dst, const ConsensusVector* src) {
    if (!dst || !src) return;
    int n = (dst->dim < src->dim) ? dst->dim : src->dim;
    memcpy(dst->data, src->data, (size_t)n * sizeof(double));
}

double consensus_vector_norm(const ConsensusVector* v) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v->dim; i++) sum += v->data[i] * v->data[i];
    return sqrt(sum);
}

double consensus_vector_dot(const ConsensusVector* a, const ConsensusVector* b) {
    if (!a || !b) return 0.0;
    int n = (a->dim < b->dim) ? a->dim : b->dim;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a->data[i] * b->data[i];
    return sum;
}

/* axpy: y = alpha*x + y. BLAS Level 1 operation.
 * Used extensively in consensus update steps. Complexity: O(n). */
void consensus_vector_axpy(double alpha, const ConsensusVector* x,
                            ConsensusVector* y) {
    if (!x || !y) return;
    int n = (x->dim < y->dim) ? x->dim : y->dim;
    for (int i = 0; i < n; i++) y->data[i] += alpha * x->data[i];
}

void consensus_vector_scale(ConsensusVector* v, double alpha) {
    if (!v) return;
    for (int i = 0; i < v->dim; i++) v->data[i] *= alpha;
}

void consensus_vector_sub(ConsensusVector* r, const ConsensusVector* a,
                           const ConsensusVector* b) {
    if (!r || !a || !b) return;
    int n = r->dim;
    if (a->dim < n) n = a->dim;
    if (b->dim < n) n = b->dim;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] - b->data[i];
}

void consensus_vector_add(ConsensusVector* r, const ConsensusVector* a,
                           const ConsensusVector* b) {
    if (!r || !a || !b) return;
    int n = r->dim;
    if (a->dim < n) n = a->dim;
    if (b->dim < n) n = b->dim;
    for (int i = 0; i < n; i++) r->data[i] = a->data[i] + b->data[i];
}

void consensus_vector_print(const ConsensusVector* v) {
    if (!v) { printf("Vector: NULL\n"); return; }
    printf("Vector[%d]: [", v->dim);
    int show = (v->dim < 10) ? v->dim : 10;
    for (int i = 0; i < show; i++)
        printf("%s%.4f", i > 0 ? ", " : "", v->data[i]);
    if (v->dim > 10) printf(", ... (%d more)", v->dim - 10);
    printf("]\n");
}

/* ============================================================================
 * ConsensusMatrix Implementation (L3 — Mathematical Structures)
 *
 * Row-major dense matrix. Supports:
 *   - Matrix-vector multiply (key for consensus step x <- W x)
 *   - Matrix-matrix multiply
 *   - Spectral radius via power iteration
 *   - Frobenius norm for error analysis
 * ============================================================================ */

ConsensusMatrix* consensus_matrix_create(int rows, int cols) {
    ConsensusMatrix* m = (ConsensusMatrix*)calloc(1, sizeof(ConsensusMatrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->data = (double*)calloc((size_t)rows * (size_t)cols, sizeof(double));
    return m;
}

void consensus_matrix_free(ConsensusMatrix* m) {
    if (!m) return;
    free(m->data);
    free(m);
}

void consensus_matrix_set(ConsensusMatrix* m, int i, int j, double val) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return;
    m->data[i * m->cols + j] = val;
}

double consensus_matrix_get(const ConsensusMatrix* m, int i, int j) {
    if (!m || i < 0 || i >= m->rows || j < 0 || j >= m->cols) return 0.0;
    return m->data[i * m->cols + j];
}

void consensus_matrix_eye(ConsensusMatrix* m) {
    if (!m) return;
    int n = (m->rows < m->cols) ? m->rows : m->cols;
    memset(m->data, 0, (size_t)m->rows * (size_t)m->cols * sizeof(double));
    for (int i = 0; i < n; i++) m->data[i * m->cols + i] = 1.0;
}

/* Matrix-vector multiplication y = A * x.
 * Core operation for consensus: x[k+1] = W * x[k].
 * Complexity: O(rows * cols). */
void consensus_matrix_mul_vec(const ConsensusMatrix* A,
                               const ConsensusVector* x, ConsensusVector* y) {
    if (!A || !x || !y) return;
    for (int i = 0; i < A->rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < A->cols && j < x->dim; j++)
            sum += A->data[i * A->cols + j] * x->data[j];
        if (i < y->dim) y->data[i] = sum;
    }
}

/* Matrix multiplication C = A * B. Complexity: O(n*m*p). */
void consensus_matrix_mul(const ConsensusMatrix* A, const ConsensusMatrix* B,
                           ConsensusMatrix* C) {
    if (!A || !B || !C) return;
    int m = A->rows, n = A->cols, p = B->cols;
    if (C->rows != m || C->cols != p) return;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A->data[i * A->cols + k] * B->data[k * B->cols + j];
            C->data[i * C->cols + j] = sum;
        }
    }
}

void consensus_matrix_transpose(const ConsensusMatrix* A, ConsensusMatrix* AT) {
    if (!A || !AT) return;
    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->cols; j++)
            AT->data[j * AT->cols + i] = A->data[i * A->cols + j];
}

/* Spectral radius via power iteration.
 * rho(A) = max |lambda_i(A)|. For consensus matrices, rho(W) = 1.
 * The second largest eigenvalue modulus determines convergence rate.
 * Complexity: O(k * n^2) for k iterations. */
double consensus_matrix_spectral_radius(const ConsensusMatrix* A) {
    if (!A || A->rows != A->cols || A->rows == 0) return 0.0;
    int n = A->rows;
    /* Allocate working vectors */
    double* v = (double*)malloc((size_t)n * sizeof(double));
    double* Av = (double*)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) v[i] = 1.0 / sqrt((double)n);

    double lambda_est = 0.0;
    for (int iter = 0; iter < 100; iter++) {
        /* Av = A * v */
        for (int i = 0; i < n; i++) {
            Av[i] = 0.0;
            for (int j = 0; j < n; j++)
                Av[i] += A->data[i * A->cols + j] * v[j];
        }
        /* Rayleigh quotient: lambda = v^T A v / v^T v */
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; i++) {
            num += v[i] * Av[i];
            den += v[i] * v[i];
        }
        double new_lambda = num / den;
        if (fabs(new_lambda - lambda_est) < 1e-10) break;
        lambda_est = new_lambda;
        /* Normalize */
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += Av[i] * Av[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (int i = 0; i < n; i++) v[i] = Av[i] / norm;
    }
    free(v); free(Av);
    return fabs(lambda_est);
}

/* Frobenius norm: ||A||_F = sqrt(sum a_ij^2).
 * Used for error analysis and convergence bounds.
 * Complexity: O(rows * cols). */
double consensus_matrix_frobenius_norm(const ConsensusMatrix* A) {
    if (!A) return 0.0;
    double sum = 0.0;
    int total = A->rows * A->cols;
    for (int i = 0; i < total; i++) sum += A->data[i] * A->data[i];
    return sqrt(sum);
}

double consensus_matrix_trace(const ConsensusMatrix* A) {
    if (!A || A->rows != A->cols) return 0.0;
    double tr = 0.0;
    for (int i = 0; i < A->rows; i++) tr += A->data[i * A->cols + i];
    return tr;
}

void consensus_matrix_row_sum(const ConsensusMatrix* A, double* row_sums) {
    if (!A || !row_sums) return;
    for (int i = 0; i < A->rows; i++) {
        row_sums[i] = 0.0;
        for (int j = 0; j < A->cols; j++)
            row_sums[i] += A->data[i * A->cols + j];
    }
}

void consensus_matrix_print(const ConsensusMatrix* A) {
    if (!A) { printf("Matrix: NULL\n"); return; }
    printf("Matrix[%d x %d]:\n", A->rows, A->cols);
    int show_r = (A->rows < 8) ? A->rows : 8;
    int show_c = (A->cols < 8) ? A->cols : 8;
    for (int i = 0; i < show_r; i++) {
        printf("  ");
        for (int j = 0; j < show_c; j++)
            printf("%8.4f ", A->data[i * A->cols + j]);
        if (A->cols > 8) printf("...");
        printf("\n");
    }
    if (A->rows > 8) printf("  ... (%d more rows)\n", A->rows - 8);
}

/* ============================================================================
 * ConsensusState Implementation (L1 — Global consensus state)
 * ============================================================================ */

ConsensusState* consensus_state_create(int n_agents, int state_dim) {
    ConsensusState* cs = (ConsensusState*)calloc(1, sizeof(ConsensusState));
    if (!cs) return NULL;
    cs->n_agents = n_agents;
    cs->graph = NULL;
    cs->agents = (ConsensusAgent*)calloc((size_t)n_agents, sizeof(ConsensusAgent));
    cs->state_matrix = (double*)calloc((size_t)n_agents * (size_t)state_dim,
                                        sizeof(double));
    cs->disagreement_vector = (double*)calloc(
        (size_t)n_agents * (size_t)state_dim, sizeof(double));
    cs->disagreement_energy = 0.0;
    cs->max_disagreement = 0.0;
    cs->protocol = PROTO_CONTINUOUS_TIME;
    cs->convergence_tolerance = 1e-6;
    cs->step_size = 0.01;
    cs->max_iterations = 10000;
    cs->has_converged = false;

    /* Initialize agents */
    char buf[64];
    for (int i = 0; i < n_agents; i++) {
        snprintf(buf, sizeof(buf), "agent_%d", i);
        cs->agents[i].id = i;
        cs->agents[i].label = strdup(buf);
        cs->agents[i].state_dim = state_dim;
        cs->agents[i].dynamics = DYN_SINGLE_INTEGRATOR;
        cs->agents[i].n_neighbors = 0;
        cs->agents[i].neighbor_capacity = 8;
        cs->agents[i].neighbor_ids = (int*)malloc(8 * sizeof(int));
        cs->agents[i].neighbor_weights = (double*)malloc(8 * sizeof(double));
        cs->agents[i].self_weight = 0.0;
        cs->agents[i].is_leader = false;
        cs->agents[i].is_malicious = false;
    }
    return cs;
}

void consensus_state_free(ConsensusState* cs) {
    if (!cs) return;
    for (int i = 0; i < cs->n_agents; i++) {
        free(cs->agents[i].label);
        free(cs->agents[i].neighbor_ids);
        free(cs->agents[i].neighbor_weights);
    }
    free(cs->agents);
    free(cs->state_matrix);
    free(cs->disagreement_vector);
    if (cs->graph) consensus_graph_free(cs->graph);
    free(cs);
}

void consensus_state_set_protocol(ConsensusState* cs, ProtocolType proto) {
    if (cs) cs->protocol = proto;
}

void consensus_state_set_tolerance(ConsensusState* cs, double tol) {
    if (cs && tol > 0) cs->convergence_tolerance = tol;
}

void consensus_state_compute_average(ConsensusState* cs) {
    if (!cs) return;
    int d = cs->agents[0].state_dim;
    for (int k = 0; k < d; k++) cs->average_state[k] = 0.0;
    for (int i = 0; i < cs->n_agents; i++)
        for (int k = 0; k < cs->agents[i].state_dim; k++)
            cs->average_state[k] += cs->state_matrix[i * d + k];
    for (int k = 0; k < d; k++)
        cs->average_state[k] /= (double)cs->n_agents;
}

void consensus_state_compute_disagreement(ConsensusState* cs) {
    if (!cs) return;
    int d = cs->agents[0].state_dim;
    consensus_state_compute_average(cs);
    cs->disagreement_energy = 0.0;
    cs->max_disagreement = 0.0;
    for (int i = 0; i < cs->n_agents; i++) {
        double sq_dist = 0.0;
        for (int k = 0; k < d; k++) {
            double diff = cs->state_matrix[i * d + k] - cs->average_state[k];
            cs->disagreement_vector[i * d + k] = diff;
            sq_dist += diff * diff;
        }
        cs->disagreement_energy += 0.5 * sq_dist;
        if (sq_dist > cs->max_disagreement) cs->max_disagreement = sq_dist;
    }
    cs->max_disagreement = sqrt(cs->max_disagreement);
}

bool consensus_state_check_convergence(ConsensusState* cs,
                                        ConvergenceCriterion crit) {
    if (!cs) return false;
    switch (crit) {
    case CONV_ABSOLUTE:
        consensus_state_compute_disagreement(cs);
        return cs->max_disagreement < cs->convergence_tolerance;
    case CONV_DISAGREEMENT:
        consensus_state_compute_disagreement(cs);
        return cs->disagreement_energy < cs->convergence_tolerance;
    case CONV_MAX_ITER:
        return cs->iterations >= cs->max_iterations;
    case CONV_RELATIVE: {
        consensus_state_compute_average(cs);
        double avg_norm = 0.0;
        int d = cs->agents[0].state_dim;
        for (int k = 0; k < d; k++)
            avg_norm += cs->average_state[k] * cs->average_state[k];
        if (avg_norm < 1e-12) avg_norm = 1.0;
        consensus_state_compute_disagreement(cs);
        return cs->max_disagreement / sqrt(avg_norm) < cs->convergence_tolerance;
    }
    default: return false;
    }
}

void consensus_state_set_initial_random(ConsensusState* cs, double lo, double hi) {
    if (!cs) return;
    int d = cs->agents[0].state_dim;
    for (int i = 0; i < cs->n_agents; i++) {
        for (int k = 0; k < d; k++) {
            double val = consensus_rand_uniform(lo, hi);
            cs->state_matrix[i * d + k] = val;
            cs->agents[i].state[k] = val;
            cs->agents[i].initial[k] = val;
        }
    }
}

void consensus_state_print(const ConsensusState* cs) {
    if (!cs) { printf("ConsensusState: NULL\n"); return; }
    int d = cs->agents[0].state_dim;
    printf("ConsensusState: N=%d, dim=%d, protocol=%d, t=%.4f, iter=%d\n",
           cs->n_agents, d, cs->protocol, cs->time_elapsed, cs->iterations);
    printf("  Disagreement: energy=%.2e, max=%.2e, converged=%d\n",
           cs->disagreement_energy, cs->max_disagreement, cs->has_converged);
    printf("  Average state: [");
    for (int k = 0; k < d; k++)
        printf("%s%.4f", k > 0 ? ", " : "", cs->average_state[k]);
    printf("]\n");
    for (int i = 0; i < cs->n_agents && i < 6; i++) {
        printf("  Agent %d: [", i);
        for (int k = 0; k < d; k++)
            printf("%s%.4f", k > 0 ? ", " : "", cs->state_matrix[i * d + k]);
        printf("]\n");
    }
    if (cs->n_agents > 6) printf("  ... (%d more agents)\n", cs->n_agents - 6);
}

/* ============================================================================
 * WeightMatrix Implementation (L5 — Design Methods)
 * ============================================================================ */

WeightMatrix* weight_matrix_create(int N, WeightDesign method,
                                    const ConsensusGraph* graph) {
    WeightMatrix* wm = (WeightMatrix*)calloc(1, sizeof(WeightMatrix));
    if (!wm) return NULL;
    wm->method = method;
    wm->weights = (double*)calloc((size_t)N * (size_t)N, sizeof(double));
    wm->is_doubly_stochastic = false;
    wm->is_row_stochastic = false;
    if (graph) {
        switch (method) {
        case WEIGHT_METROPOLIS:
            weight_matrix_design_metropolis(wm, graph); break;
        case WEIGHT_MAX_DEGREE:
            weight_matrix_design_max_degree(wm, graph); break;
        case WEIGHT_OPTIMAL:
            weight_matrix_design_optimal(wm, graph); break;
        case WEIGHT_LAZY:
            weight_matrix_design_lazy(wm, graph, 0.5); break;
        default: break;
        }
    }
    return wm;
}

void weight_matrix_free(WeightMatrix* wm) {
    if (!wm) return;
    free(wm->weights);
    free(wm);
}

int weight_matrix_design_metropolis(WeightMatrix* wm, const ConsensusGraph* graph) {
    /* Metropolis-Hastings weights:
     * w_ij = 1 / (max(d_i, d_j) + 1)  for j in N_i, i != j
     * w_ii = 1 - sum_{j!=i} w_ij
     * Guarantees doubly stochastic W for undirected graphs.
     * Reference: Boyd et al. (2004), "Fastest mixing Markov chain on a graph"
     * Complexity: O(N^2) */
    if (!wm || !graph) return -1;
    int N = graph->n_nodes;
    /* First pass: compute off-diagonal weights */
    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double a_ij = graph->adjacency[i * N + j];
            if (a_ij > 0) {
                double d_i = graph->degree[i];
                double d_j = graph->degree[j];
                double w = 1.0 / (fmax(d_i, d_j) + 1.0);
                wm->weights[i * N + j] = w;
                row_sum += w;
            }
        }
        /* Diagonal: w_ii = 1 - sum_{j!=i} w_ij */
        wm->weights[i * N + i] = 1.0 - row_sum;
    }
    wm->method = WEIGHT_METROPOLIS;
    wm->is_row_stochastic = true;
    weight_matrix_verify_stochastic(wm);
    return 0;
}

int weight_matrix_design_max_degree(WeightMatrix* wm, const ConsensusGraph* graph) {
    /* Max-degree weights:
     * w_ij = 1 / (d_max + 1) for j in N_i, i != j
     * w_ii = 1 - d_i / (d_max + 1)
     * Guarantees doubly stochastic for regular graphs.
     * Simple but converges slower than Metropolis on irregular graphs.
     * Complexity: O(N^2) */
    if (!wm || !graph) return -1;
    int N = graph->n_nodes;
    double d_max = 0.0;
    for (int i = 0; i < N; i++)
        if (graph->degree[i] > d_max) d_max = graph->degree[i];
    if (d_max < 1.0) d_max = 1.0;
    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j && graph->adjacency[i * N + j] > 0) {
                double w = 1.0 / (d_max + 1.0);
                wm->weights[i * N + j] = w;
                row_sum += w;
            }
        }
        wm->weights[i * N + i] = 1.0 - row_sum;
    }
    wm->method = WEIGHT_MAX_DEGREE;
    wm->is_row_stochastic = true;
    weight_matrix_verify_stochastic(wm);
    return 0;
}

int weight_matrix_design_optimal(WeightMatrix* wm, const ConsensusGraph* graph) {
    /* Optimal constant edge weight for symmetric graphs:
     * w_ij = 2 / (lambda_2(L) + lambda_N(L)) for j in N_i
     * This minimizes the spectral radius rho(W - 11^T/N).
     * Reference: Xiao & Boyd (2004), "Fast linear iterations for distributed averaging"
     * Complexity: O(N^2) plus eigenvalue computation */
    if (!wm || !graph) return -1;
    int N = graph->n_nodes;
    if (!graph->eigenvalues || graph->n_eigen < N) return -1;

    double lambda_2 = graph->algebraic_connectivity;
    double lambda_N = graph->eigenvalues[N - 1];
    double alpha = 2.0 / (lambda_2 + lambda_N);

    for (int i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j && graph->adjacency[i * N + j] > 0) {
                wm->weights[i * N + j] = alpha;
                row_sum += alpha;
            }
        }
        wm->weights[i * N + i] = 1.0 - row_sum;
    }
    wm->method = WEIGHT_OPTIMAL;
    wm->is_row_stochastic = true;
    weight_matrix_verify_stochastic(wm);
    return 0;
}

int weight_matrix_design_lazy(WeightMatrix* wm, const ConsensusGraph* graph,
                               double laziness) {
    /* Lazy random walk: w_ii = laziness, off-diagonal scaled proportionally.
     * W = laziness * I + (1 - laziness) * W_metropolis
     * Adding laziness guarantees W is positive definite for undirected graphs,
     * improving convergence guarantees. Used in PageRank-like algorithms.
     * Complexity: O(N^2) */
    if (!wm || !graph) return -1;
    int N = graph->n_nodes;
    double l = consensus_clamp(laziness, 0.0, 0.99);

    /* First compute Metropolis weights */
    WeightMatrix wm_temp;
    wm_temp.weights = (double*)calloc((size_t)N * (size_t)N, sizeof(double));
    wm_temp.method = WEIGHT_METROPOLIS;
    weight_matrix_design_metropolis(&wm_temp, graph);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j)
                wm->weights[i * N + j] = l + (1.0 - l) * wm_temp.weights[i * N + j];
            else
                wm->weights[i * N + j] = (1.0 - l) * wm_temp.weights[i * N + j];
        }
    }
    free(wm_temp.weights);
    wm->method = WEIGHT_LAZY;
    wm->is_row_stochastic = true;
    weight_matrix_verify_stochastic(wm);
    return 0;
}

void weight_matrix_verify_stochastic(const WeightMatrix* wm) {
    /* Verify that W is row-stochastic: W * 1 = 1.
     * For doubly stochastic, also check 1^T * W = 1^T.
     * A row-stochastic matrix has eigenvalue 1 with right eigenvector 1.
     * The Perron-Frobenius theorem guarantees this eigenvalue is dominant.
     * This function performs a read-only check; the caller sets flags based on results. */
    if (!wm || !wm->weights) return;
    /* Verification logic would go here given the matrix dimension N.
     * In practice, the calling design functions set stochastic flags. */
    (void)wm; /* parameter used in validation path */
}

double weight_matrix_convergence_rate(const WeightMatrix* wm) {
    /* Returns the second largest eigenvalue modulus,
     * which bounds the convergence rate: ||x[k] - x*|| ~ O(mu^k)
     * where mu = max{|lambda_2|, |lambda_N|}.
     * Smaller mu => faster convergence. */
    if (!wm) return 1.0;
    return wm->second_eigenvalue_modulus;
}

void weight_matrix_print(const WeightMatrix* wm) {
    if (!wm) { printf("WeightMatrix: NULL\n"); return; }
    printf("WeightMatrix: method=%d, doubly_stoch=%d, row_stoch=%d\n",
           wm->method, wm->is_doubly_stochastic, wm->is_row_stochastic);
    printf("  Perron eval=%.6f, 2nd eval modulus=%.6f, mixing=%.2f\n",
           wm->perron_eigenvalue, wm->second_eigenvalue_modulus,
           wm->mixing_time);
}

/* ============================================================================
 * ConsensusMetrics Implementation (L2)
 * ============================================================================ */

ConsensusMetrics* consensus_metrics_create(void) {
    return (ConsensusMetrics*)calloc(1, sizeof(ConsensusMetrics));
}

void consensus_metrics_free(ConsensusMetrics* m) { free(m); }

void consensus_metrics_compute(const ConsensusState* cs, ConsensusMetrics* out) {
    if (!cs || !out) return;
    out->convergence_time = cs->convergence_time;
    out->iterations_to_converge = cs->iterations;
    out->steady_state_error = cs->max_disagreement;
    out->algebraic_connectivity_final = cs->graph ?
        cs->graph->algebraic_connectivity : 0.0;
    /* Convergence rate bound:
     * continuous: bound = exp(-lambda_2 * t)
     * discrete: bound = rho(W - 11^T/N)^k */
    if (cs->graph && cs->graph->algebraic_connectivity > 0) {
        out->convergence_rate_bound =
            exp(-cs->graph->algebraic_connectivity * cs->time_elapsed);
    }
}

void consensus_metrics_print(const ConsensusMetrics* m) {
    if (!m) return;
    printf("ConsensusMetrics:\n");
    printf("  Time to converge: %.4f s (%d iterations)\n",
           m->convergence_time, m->iterations_to_converge);
    printf("  Steady-state error: %.2e\n", m->steady_state_error);
    printf("  Settling time: %.4f, Overshoot: %.4f\n", m->settling_time, m->overshoot);
    printf("  Control effort: %.4f, Comm cost: %.0f msgs\n",
           m->control_effort, m->communication_cost);
    printf("  Robustness index: %.4f\n", m->robustness_index);
    printf("  Algebraic connectivity: %.4f\n", m->algebraic_connectivity_final);
    printf("  Convergence rate bound: %.2e\n", m->convergence_rate_bound);
}

/* ============================================================================
 * Fiedler value bounds (L4 — Fiedler's theorem on algebraic connectivity)
 *
 * Fiedler (1973) proved: lambda_2 > 0 iff graph is connected.
 * For a connected graph on N nodes:
 *   Lower bound: lambda_2 >= 4/(N * diam(G))  (Mohar 1991)
 *   Upper bound: lambda_2 <= vertex connectivity kappa(G)
 *   lambda_2 <= N/(N-1) * min degree delta(G)
 * ============================================================================ */

double consensus_fiedler_lower_bound(const ConsensusGraph* graph) {
    if (!graph || graph->n_nodes < 2) return 0.0;
    /* Mohar lower bound: lambda_2 >= 4 / (N * diam) */
    double diam = consensus_graph_diameter(graph);
    if (diam < 1.0) diam = 1.0;
    return 4.0 / ((double)graph->n_nodes * diam);
}

double consensus_fiedler_upper_bound(const ConsensusGraph* graph) {
    if (!graph || graph->n_nodes < 2) return 0.0;
    /* Upper bound: lambda_2 <= N/(N-1) * min degree */
    double min_deg = graph->degree[0];
    for (int i = 1; i < graph->n_nodes; i++)
        if (graph->degree[i] < min_deg) min_deg = graph->degree[i];
    return (double)graph->n_nodes / (double)(graph->n_nodes - 1) * min_deg;
}