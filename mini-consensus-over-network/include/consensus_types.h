#ifndef CONSENSUS_TYPES_H
#define CONSENSUS_TYPES_H
#include <stdbool.h>
#include <stddef.h>
/* ============================================================================
 * Consensus over Network - Core Type Definitions (L1 Definitions)
 * References: Olfati-Saber & Murray (2004), Ren & Beard (2008),
 *   Mesbahi & Egerstedt (2010), Jadbabaie-Lin-Morse (2003),
 *   Fiedler (1973), DeGroot (1974), Tsitsiklis (1984)
 * ============================================================================ */
#define CONSENSUS_MAX_STATE_DIM  16
#define CONSENSUS_MAX_AGENTS     1024
#define CONSENSUS_MAX_NEIGHBORS  256
#define CONSENSUS_EPSILON        1e-12

typedef enum { TOPO_UNDIRECTED=0, TOPO_DIRECTED=1, TOPO_BIDIRECTIONAL=2, TOPO_TIME_VARYING=3 } TopologyType;
typedef enum { PROTO_CONTINUOUS_TIME=0, PROTO_DISCRETE_TIME=1, PROTO_GOSSIP=2, PROTO_EVENT_TRIGGERED=3, PROTO_FINITE_TIME=4 } ProtocolType;
typedef enum { WEIGHT_UNIFORM=0, WEIGHT_METROPOLIS=1, WEIGHT_MAX_DEGREE=2, WEIGHT_LAPLACIAN=3, WEIGHT_OPTIMAL=4, WEIGHT_LAZY=5 } WeightDesign;
typedef enum { CONV_ABSOLUTE=0, CONV_RELATIVE=1, CONV_DISAGREEMENT=2, CONV_MAX_ITER=3 } ConvergenceCriterion;
typedef enum { DYN_SINGLE_INTEGRATOR=0, DYN_DOUBLE_INTEGRATOR=1, DYN_GENERAL_LINEAR=2, DYN_NONLINEAR=3, DYN_UNICYCLE=4 } AgentDynamicsType;
typedef enum { NET_FIXED=0, NET_SWITCHING=1, NET_RANDOM=2, NET_PROXIMITY=3, NET_DELAYED=4, NET_PACKET_LOSS=5 } NetworkModel;
typedef enum { LAPLACIAN_COMBINATORIAL=0, LAPLACIAN_NORMALIZED=1, LAPLACIAN_SIGNLESS=2, LAPLACIAN_RANDOM_WALK=3 } LaplacianType;

typedef struct { double* data; int dim; } ConsensusVector;
typedef struct { double* data; int rows; int cols; } ConsensusMatrix;

typedef struct {
    int id; char* label;
    double state[16]; double input[16]; double initial[16];
    int state_dim; AgentDynamicsType dynamics;
    double* neighbor_weights; int* neighbor_ids; int n_neighbors; int neighbor_capacity;
    double self_weight; bool is_leader; bool is_malicious;
} ConsensusAgent;

typedef struct { int from, to; double weight; bool active; } GraphEdge;

typedef struct {
    int n_nodes, n_edges; GraphEdge* edges; int edge_capacity;
    TopologyType topology; NetworkModel net_model;
    double *adjacency, *degree, *laplacian, *normalized_laplacian;
    double *perron_matrix, *incidence;
    bool is_connected, is_balanced, has_spanning_tree;
    double algebraic_connectivity, spectral_radius;
    double* eigenvalues; int n_eigen;
} ConsensusGraph;

typedef struct {
    ConsensusAgent* agents; int n_agents; ConsensusGraph* graph;
    double* state_matrix; double* disagreement_vector;
    double average_state[16]; double disagreement_energy;
    double max_disagreement; double consensus_value[16];
    double time_elapsed; int iterations; bool has_converged;
    double convergence_time; ProtocolType protocol; LaplacianType lap_type;
    double convergence_tolerance, step_size; int max_iterations;
} ConsensusState;

typedef struct {
    WeightDesign method; double* weights;
    bool is_doubly_stochastic, is_row_stochastic;
    double perron_eigenvalue, second_eigenvalue_modulus, mixing_time;
} WeightMatrix;

typedef struct {
    double convergence_time; int iterations_to_converge;
    double steady_state_error, settling_time, overshoot;
    double control_effort, communication_cost, robustness_index;
    double algebraic_connectivity_final, convergence_rate_bound;
} ConsensusMetrics;

ConsensusAgent* consensus_agent_create(int id, const char* label, int state_dim);
void consensus_agent_free(ConsensusAgent* agent);
void consensus_agent_set_state(ConsensusAgent* agent, const double* state, int dim);
void consensus_agent_add_neighbor(ConsensusAgent* agent, int neighbor_id, double weight);
void consensus_agent_reset(ConsensusAgent* agent);
double consensus_agent_state_norm(const ConsensusAgent* agent);
void consensus_agent_apply_input(ConsensusAgent* agent, double dt);
void consensus_agent_print(const ConsensusAgent* agent);

ConsensusVector* consensus_vector_create(int dim);
void consensus_vector_free(ConsensusVector* v);
void consensus_vector_set(ConsensusVector* v, double val);
void consensus_vector_copy(ConsensusVector* dst, const ConsensusVector* src);
double consensus_vector_norm(const ConsensusVector* v);
double consensus_vector_dot(const ConsensusVector* a, const ConsensusVector* b);
void consensus_vector_axpy(double alpha, const ConsensusVector* x, ConsensusVector* y);
void consensus_vector_scale(ConsensusVector* v, double alpha);
void consensus_vector_sub(ConsensusVector* r, const ConsensusVector* a, const ConsensusVector* b);
void consensus_vector_add(ConsensusVector* r, const ConsensusVector* a, const ConsensusVector* b);
void consensus_vector_print(const ConsensusVector* v);

ConsensusMatrix* consensus_matrix_create(int rows, int cols);
void consensus_matrix_free(ConsensusMatrix* m);
void consensus_matrix_set(ConsensusMatrix* m, int i, int j, double val);
double consensus_matrix_get(const ConsensusMatrix* m, int i, int j);
void consensus_matrix_eye(ConsensusMatrix* m);
void consensus_matrix_mul_vec(const ConsensusMatrix* A, const ConsensusVector* x, ConsensusVector* y);
void consensus_matrix_mul(const ConsensusMatrix* A, const ConsensusMatrix* B, ConsensusMatrix* C);
void consensus_matrix_transpose(const ConsensusMatrix* A, ConsensusMatrix* AT);
double consensus_matrix_spectral_radius(const ConsensusMatrix* A);
double consensus_matrix_frobenius_norm(const ConsensusMatrix* A);
double consensus_matrix_trace(const ConsensusMatrix* A);
void consensus_matrix_row_sum(const ConsensusMatrix* A, double* row_sums);
void consensus_matrix_print(const ConsensusMatrix* A);

ConsensusState* consensus_state_create(int n_agents, int state_dim);
void consensus_state_free(ConsensusState* cs);
void consensus_state_set_protocol(ConsensusState* cs, ProtocolType proto);
void consensus_state_set_tolerance(ConsensusState* cs, double tol);
void consensus_state_compute_average(ConsensusState* cs);
void consensus_state_compute_disagreement(ConsensusState* cs);
bool consensus_state_check_convergence(ConsensusState* cs, ConvergenceCriterion crit);
void consensus_state_set_initial_random(ConsensusState* cs, double lo, double hi);
void consensus_state_print(const ConsensusState* cs);

WeightMatrix* weight_matrix_create(int N, WeightDesign method, const ConsensusGraph* graph);
void weight_matrix_free(WeightMatrix* wm);
int weight_matrix_design_metropolis(WeightMatrix* wm, const ConsensusGraph* graph);
int weight_matrix_design_max_degree(WeightMatrix* wm, const ConsensusGraph* graph);
int weight_matrix_design_optimal(WeightMatrix* wm, const ConsensusGraph* graph);
int weight_matrix_design_lazy(WeightMatrix* wm, const ConsensusGraph* graph, double laziness);
void weight_matrix_verify_stochastic(const WeightMatrix* wm);
double weight_matrix_convergence_rate(const WeightMatrix* wm);
void weight_matrix_print(const WeightMatrix* wm);

ConsensusMetrics* consensus_metrics_create(void);
void consensus_metrics_free(ConsensusMetrics* m);
void consensus_metrics_compute(const ConsensusState* cs, ConsensusMetrics* out);
void consensus_metrics_print(const ConsensusMetrics* m);

double consensus_clamp(double x, double lo, double hi);
double consensus_rand_uniform(double a, double b);
void consensus_seed(unsigned int seed);
double consensus_fiedler_lower_bound(const ConsensusGraph* graph);
double consensus_fiedler_upper_bound(const ConsensusGraph* graph);
#endif