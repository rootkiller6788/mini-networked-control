#ifndef CONSENSUS_GRAPH_H
#define CONSENSUS_GRAPH_H
#include "consensus_types.h"
/* ============================================================================
 * Consensus Graph Theory (L3 - Mathematical Structures, L4 - Theorems)
 *
 * Graph-theoretic foundations for consensus:
 *   G = (V,E), Laplacian L = D - A, algebraic connectivity lambda_2(L)
 *
 * Key theorems:
 *   Thm 1 (Olfati-Saber 2004): Continuous consensus xdot = -L x converges
 *        to average iff G is connected (undirected) or has spanning tree (directed).
 *   Thm 2 (Fiedler 1973): lambda_2(L) > 0 iff G is connected.
 *   Thm 3 (Perron-Frobenius): W stochastic => rho(W) = 1, limit exists.
 *   Thm 4 (Ren & Beard 2005): Discrete consensus x[k+1]=Wx[k] converges
 *        iff W is primitive and doubly stochastic for average consensus.
 *   Thm 5 (Convergence Rate): ||x(t)-x*|| <= ||x(0)-x*|| exp(-lambda_2 t)
 *        for continuous time. For discrete: ||x[k]-x*|| <= rho(W-11^T/N)^k.
 * ============================================================================ */

ConsensusGraph* consensus_graph_create(int n_nodes, TopologyType topo);
void consensus_graph_free(ConsensusGraph* g);
int consensus_graph_add_edge(ConsensusGraph* g, int from, int to, double weight);
int consensus_graph_add_undirected_edge(ConsensusGraph* g, int a, int b, double weight);
void consensus_graph_remove_edge(ConsensusGraph* g, int from, int to);
void consensus_graph_build_adjacency(ConsensusGraph* g);
void consensus_graph_build_degree(ConsensusGraph* g);
void consensus_graph_build_laplacian(ConsensusGraph* g, LaplacianType type);
void consensus_graph_build_perron_matrix(ConsensusGraph* g, WeightDesign method);
void consensus_graph_build_incidence(ConsensusGraph* g);
void consensus_graph_compute_eigenvalues(ConsensusGraph* g);
bool consensus_graph_is_connected(const ConsensusGraph* g);
bool consensus_graph_has_spanning_tree(const ConsensusGraph* g);
int consensus_graph_connected_components(const ConsensusGraph* g, int** comp_map);
double consensus_graph_algebraic_connectivity(ConsensusGraph* g);
double consensus_graph_spectral_gap(const ConsensusGraph* g);
double consensus_graph_average_path_length(const ConsensusGraph* g);
double consensus_graph_clustering_coefficient(const ConsensusGraph* g);
double consensus_graph_diameter(const ConsensusGraph* g);
void consensus_graph_set_proximity(ConsensusGraph* g, const double* positions, int dim, double radius);
ConsensusGraph* consensus_graph_create_path(int N);
ConsensusGraph* consensus_graph_create_cycle(int N);
ConsensusGraph* consensus_graph_create_star(int N);
ConsensusGraph* consensus_graph_create_complete(int N);
ConsensusGraph* consensus_graph_create_erdos_renyi(int N, double p);
ConsensusGraph* consensus_graph_create_small_world(int N, int k, double p_rewire);
ConsensusGraph* consensus_graph_create_barabasi_albert(int N, int m0, int m);
void consensus_graph_print(const ConsensusGraph* g);

/* Laplacian spectrum bounds */
double consensus_laplacian_spectral_radius_bound(const ConsensusGraph* g);
double consensus_fiedler_mohar_lower(const ConsensusGraph* g);
double consensus_fiedler_mohar_upper(const ConsensusGraph* g);
int consensus_graph_algebraic_connectivity_perturbation(const ConsensusGraph* g, int add_edge_u, int add_edge_v);

/* Peripheral node detection (bottleneck identification) */
int consensus_graph_find_bottleneck_edge(const ConsensusGraph* g, int* u, int* v);

#endif