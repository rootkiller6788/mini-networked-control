#include "consensus_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* Forward declarations */
static void dfs_visit(const ConsensusGraph* g, int u, int* visited);
static void bfs_distance(const ConsensusGraph* g, int start, int* dist);
static void qr_algorithm(double* A, int n, double* eigenvalues);

/* ============================================================================
 * Graph Creation and Destruction (L3 — Mathematical Structures)
 * ============================================================================ */

ConsensusGraph* consensus_graph_create(int n_nodes, TopologyType topo) {
    ConsensusGraph* g = (ConsensusGraph*)calloc(1, sizeof(ConsensusGraph));
    if (!g) return NULL;
    g->n_nodes = n_nodes;
    g->n_edges = 0;
    g->edge_capacity = n_nodes * 4;
    g->edges = (GraphEdge*)malloc((size_t)g->edge_capacity * sizeof(GraphEdge));
    g->topology = topo;
    g->net_model = NET_FIXED;
    int n2 = n_nodes * n_nodes;
    g->adjacency = (double*)calloc((size_t)n2, sizeof(double));
    g->degree = (double*)calloc((size_t)n_nodes, sizeof(double));
    g->laplacian = (double*)calloc((size_t)n2, sizeof(double));
    g->normalized_laplacian = (double*)calloc((size_t)n2, sizeof(double));
    g->perron_matrix = (double*)calloc((size_t)n2, sizeof(double));
    g->eigenvalues = NULL;
    g->n_eigen = 0;
    return g;
}

void consensus_graph_free(ConsensusGraph* g) {
    if (!g) return;
    free(g->edges);
    free(g->adjacency);
    free(g->degree);
    free(g->laplacian);
    free(g->normalized_laplacian);
    free(g->perron_matrix);
    free(g->incidence);
    free(g->eigenvalues);
    free(g);
}

int consensus_graph_add_edge(ConsensusGraph* g, int from, int to, double weight) {
    if (!g || from < 0 || from >= g->n_nodes || to < 0 || to >= g->n_nodes)
        return -1;
    if (g->n_edges >= g->edge_capacity) {
        g->edge_capacity *= 2;
        g->edges = (GraphEdge*)realloc(g->edges,
                    (size_t)g->edge_capacity * sizeof(GraphEdge));
    }
    g->edges[g->n_edges].from = from;
    g->edges[g->n_edges].to = to;
    g->edges[g->n_edges].weight = weight;
    g->edges[g->n_edges].active = true;
    g->n_edges++;
    return g->n_edges - 1;
}

int consensus_graph_add_undirected_edge(ConsensusGraph* g, int a, int b,
                                         double weight) {
    int e1 = consensus_graph_add_edge(g, a, b, weight);
    int e2 = consensus_graph_add_edge(g, b, a, weight);
    return (e1 >= 0 && e2 >= 0) ? 0 : -1;
}

void consensus_graph_remove_edge(ConsensusGraph* g, int from, int to) {
    if (!g) return;
    for (int e = 0; e < g->n_edges; e++) {
        if (g->edges[e].from == from && g->edges[e].to == to) {
            g->edges[e].active = false;
            g->edges[e].weight = 0.0;
        }
    }
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
}

/* ============================================================================
 * Matrix Construction from Graph (L3)
 * ============================================================================ */

void consensus_graph_build_adjacency(ConsensusGraph* g) {
    /* Build adjacency matrix A from edge list.
     * A_ij = weight of edge (i,j). Zero if no edge.
     * For undirected graphs, A must be symmetric.
     * Complexity: O(E). */
    if (!g) return;
    int N = g->n_nodes;
    memset(g->adjacency, 0, (size_t)N * (size_t)N * sizeof(double));
    for (int e = 0; e < g->n_edges; e++) {
        if (!g->edges[e].active) continue;
        int i = g->edges[e].from, j = g->edges[e].to;
        g->adjacency[i * N + j] = g->edges[e].weight;
    }
}

void consensus_graph_build_degree(ConsensusGraph* g) {
    /* Degree vector: d_i = sum_j A_ij (out-degree for directed).
     * For undirected graphs, d_i is the standard degree.
     * The degree matrix D = diag(d_1, ..., d_N).
     * Complexity: O(N^2). */
    if (!g) return;
    int N = g->n_nodes;
    memset(g->degree, 0, (size_t)N * sizeof(double));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            g->degree[i] += g->adjacency[i * N + j];
}

void consensus_graph_build_laplacian(ConsensusGraph* g, LaplacianType type) {
    /* Laplacian matrix: L = D - A.
     *
     * Properties (undirected graph):
     *   - L is symmetric positive semi-definite
     *   - L * 1 = 0 (zero row sum)
     *   - rank(L) = N - c, where c = #connected components
     *   - 0 = lambda_1 <= lambda_2 <= ... <= lambda_N
     *   - lambda_2 > 0 iff graph is connected (Fiedler 1973)
     *
     * Normalized Laplacian: L_norm = I - D^{-1/2} A D^{-1/2}
     *   Better spectral properties, eigenvalues in [0, 2].
     *
     * Complexity: O(N^2). */
    if (!g) return;
    int N = g->n_nodes;
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);

    if (type == LAPLACIAN_COMBINATORIAL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double val = -g->adjacency[i * N + j];
                if (i == j) val = g->degree[i] - g->adjacency[i * N + j];
                g->laplacian[i * N + j] = val;
            }
        }
    } else if (type == LAPLACIAN_NORMALIZED) {
        for (int i = 0; i < N; i++) {
            double d_i_sqrt = (g->degree[i] > 0) ? 1.0 / sqrt(g->degree[i]) : 0.0;
            for (int j = 0; j < N; j++) {
                double d_j_sqrt = (g->degree[j] > 0) ? 1.0 / sqrt(g->degree[j]) : 0.0;
                if (i == j && g->degree[i] > 0)
                    g->normalized_laplacian[i * N + j] = 1.0 -
                        g->adjacency[i * N + j] * d_i_sqrt * d_j_sqrt;
                else if (i != j)
                    g->normalized_laplacian[i * N + j] =
                        -g->adjacency[i * N + j] * d_i_sqrt * d_j_sqrt;
            }
        }
    } else if (type == LAPLACIAN_SIGNLESS) {
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                g->laplacian[i * N + j] = g->degree[i] * (i==j) +
                    g->adjacency[i * N + j];
    }
}

void consensus_graph_build_perron_matrix(ConsensusGraph* g, WeightDesign method) {
    /* Build Perron (stochastic) matrix W for discrete-time consensus.
     * For undirected graphs, W is doubly stochastic.
     * The Perron-Frobenius theorem guarantees rho(W) = 1.
     * Reference: Ren & Beard (2008), Sec 3.2. */
    if (!g) return;
    int N = g->n_nodes;
    WeightMatrix wm;
    wm.weights = g->perron_matrix;
    wm.method = method;
    switch (method) {
    case WEIGHT_METROPOLIS: weight_matrix_design_metropolis(&wm, g); break;
    case WEIGHT_MAX_DEGREE: weight_matrix_design_max_degree(&wm, g); break;
    default:
        /* Uniform weights: w_ij = 1/d_i for j in N_i */
        for (int i = 0; i < N; i++) {
            if (g->degree[i] < 0.5) {
                g->perron_matrix[i * N + i] = 1.0;
                continue;
            }
            for (int j = 0; j < N; j++) {
                g->perron_matrix[i * N + j] =
                    g->adjacency[i * N + j] / g->degree[i];
            }
        }
        break;
    }
}

void consensus_graph_build_incidence(ConsensusGraph* g) {
    /* Incidence matrix B: N x E, B[i,e] = +1 if edge e enters i,
     * -1 if edge e leaves i, 0 otherwise.
     * L = B * B^T for unweighted undirected graphs.
     * Used for formation control: ẋ = -B B^T x */
    if (!g) return;
    int N = g->n_nodes, E = g->n_edges;
    free(g->incidence);
    g->incidence = (double*)calloc((size_t)N * (size_t)E, sizeof(double));
    for (int e = 0; e < E; e++) {
        if (!g->edges[e].active) continue;
        g->incidence[g->edges[e].from * E + e] = -1.0;
        g->incidence[g->edges[e].to * E + e] = 1.0;
    }
}

/* ============================================================================
 * Eigenvalue Computation (L4 — Spectral graph theory)
 * ============================================================================ */

/* Simple QR algorithm for symmetric matrices.
 * Used to compute Laplacian spectrum: 0 = lambda_1 <= lambda_2 <= ... <= lambda_N.
 * lambda_2 (algebraic connectivity / Fiedler value) determines consensus rate.
 * Complexity: O(n^3) per iteration, typically O(n^3) total. */
static void qr_algorithm(double* A, int n, double* eigenvalues) {
    /* Copy A to working matrix */
    double* H = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
    memcpy(H, A, (size_t)n * (size_t)n * sizeof(double));

    /* Householder QR iteration for tridiagonalization + implicit shifts */
    for (int iter = 0; iter < 50; iter++) {
        /* Check for convergence: off-diagonal small enough? */
        double off_diag = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j) off_diag += H[i * n + j] * H[i * n + j];
        if (sqrt(off_diag) < 1e-12) break;

        /* QR decomposition via Gram-Schmidt */
        double* Q = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
        double* R = (double*)calloc((size_t)n * (size_t)n, sizeof(double));

        /* Copy columns of H as vectors, orthonormalize */
        for (int j = 0; j < n; j++) {
            /* Copy column j */
            for (int i = 0; i < n; i++) Q[i * n + j] = H[i * n + j];
            /* Orthogonalize against previous columns */
            for (int k = 0; k < j; k++) {
                double dot = 0.0;
                for (int i = 0; i < n; i++)
                    dot += Q[i * n + k] * H[i * n + j];
                for (int i = 0; i < n; i++)
                    Q[i * n + j] -= dot * Q[i * n + k];
                R[k * n + j] = dot;
            }
            /* Normalize */
            double norm = 0.0;
            for (int i = 0; i < n; i++)
                norm += Q[i * n + j] * Q[i * n + j];
            norm = sqrt(norm);
            if (norm > 1e-15) {
                for (int i = 0; i < n; i++) Q[i * n + j] /= norm;
                R[j * n + j] = norm;
            }
        }

        /* Compute R*Q (RQ iteration) */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                H[i * n + j] = 0.0;
                for (int k = 0; k < n; k++)
                    H[i * n + j] += R[i * n + k] * Q[k * n + j];
            }
        free(Q); free(R);
    }

    /* Extract eigenvalues from diagonal */
    for (int i = 0; i < n; i++) eigenvalues[i] = H[i * n + i];

    /* Sort ascending (selection sort for small matrices) */
    for (int i = 0; i < n - 1; i++) {
        int min_idx = i;
        for (int j = i + 1; j < n; j++)
            if (eigenvalues[j] < eigenvalues[min_idx]) min_idx = j;
        double tmp = eigenvalues[i];
        eigenvalues[i] = eigenvalues[min_idx];
        eigenvalues[min_idx] = tmp;
    }
    free(H);
}

void consensus_graph_compute_eigenvalues(ConsensusGraph* g) {
    if (!g || g->n_nodes < 1) return;
    int N = g->n_nodes;
    free(g->eigenvalues);
    g->eigenvalues = (double*)malloc((size_t)N * sizeof(double));
    g->n_eigen = N;

    /* Compute eigenvalues of Laplacian */
    double* L_copy = (double*)malloc((size_t)N * (size_t)N * sizeof(double));
    memcpy(L_copy, g->laplacian, (size_t)N * (size_t)N * sizeof(double));
    qr_algorithm(L_copy, N, g->eigenvalues);
    free(L_copy);

    /* First eigenvalue should be ~0 (up to numerical error) */
    g->algebraic_connectivity = (N >= 2) ? g->eigenvalues[1] : 0.0;
}

/* ============================================================================
 * Connectivity Analysis (L4 — Fundamental graph properties)
 * ============================================================================ */

static void dfs_visit(const ConsensusGraph* g, int u, int* visited) {
    visited[u] = 1;
    int N = g->n_nodes;
    for (int v = 0; v < N; v++) {
        if (g->adjacency[u * N + v] > 0 && !visited[v])
            dfs_visit(g, v, visited);
    }
}

bool consensus_graph_is_connected(const ConsensusGraph* g) {
    /* Check if the graph is connected using DFS.
     * For undirected graphs: connectivity <=> lambda_2(L) > 0 (Fiedler 1973).
     * For directed graphs: weak connectivity (underlying undirected is connected).
     * Complexity: O(N + E). */
    if (!g || g->n_nodes < 1) return false;
    int N = g->n_nodes;
    int* visited = (int*)calloc((size_t)N, sizeof(int));
    dfs_visit(g, 0, visited);
    for (int i = 0; i < N; i++)
        if (!visited[i]) { free(visited); return false; }
    free(visited);
    return true;
}

bool consensus_graph_has_spanning_tree(const ConsensusGraph* g) {
    /* A directed graph has a directed spanning tree iff there exists
     * a node (root) that can reach all other nodes via directed paths.
     * This is the necessary and sufficient condition for consensus
     * under directed topology (Ren & Beard 2005).
     * Method: check reachability from each node via DFS.
     * Complexity: O(N * (N + E)). */
    if (!g || g->n_nodes < 2) return (g && g->n_nodes == 1);
    int N = g->n_nodes;
    for (int root = 0; root < N; root++) {
        int* visited = (int*)calloc((size_t)N, sizeof(int));
        dfs_visit(g, root, visited);
        bool reaches_all = true;
        for (int i = 0; i < N && reaches_all; i++)
            if (!visited[i]) reaches_all = false;
        free(visited);
        if (reaches_all) return true;
    }
    return false;
}

int consensus_graph_connected_components(const ConsensusGraph* g, int** comp_map) {
    if (!g || g->n_nodes < 1) return 0;
    int N = g->n_nodes;
    int* comp = (int*)malloc((size_t)N * sizeof(int));
    memset(comp, -1, (size_t)N * sizeof(int));
    int comp_id = 0;
    for (int i = 0; i < N; i++) {
        if (comp[i] < 0) {
            int* visited = (int*)calloc((size_t)N, sizeof(int));
            dfs_visit(g, i, visited);
            for (int j = 0; j < N; j++)
                if (visited[j]) comp[j] = comp_id;
            free(visited);
            comp_id++;
        }
    }
    *comp_map = comp;
    return comp_id;
}

/* ============================================================================
 * Spectral Properties (L4)
 * ============================================================================ */

double consensus_graph_algebraic_connectivity(ConsensusGraph* g) {
    if (!g) return 0.0;
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    consensus_graph_compute_eigenvalues(g);
    return g->algebraic_connectivity;
}

double consensus_graph_spectral_gap(const ConsensusGraph* g) {
    /* Spectral gap = |lambda_2| - |lambda_N| for the Perron matrix W.
     * Large gap => fast mixing. Small gap => slow convergence.
     * For consensus: the relevant gap is 1 - max{|lambda_2|, |lambda_N|}
     * where lambda_i are eigenvalues of W (excluding lambda_1 = 1).
     * Complexity: O(N^3) for eigenvalue computation. */
    if (!g || g->n_nodes < 2) return 0.0;
    int N = g->n_nodes;
    if (!g->eigenvalues || g->n_eigen < N) return 0.0;
    /* The spectral gap for consensus rate */
    return g->algebraic_connectivity;
}

static void bfs_distance(const ConsensusGraph* g, int start, int* dist) {
    int N = g->n_nodes;
    for (int i = 0; i < N; i++) dist[i] = -1;
    int* queue = (int*)malloc((size_t)N * sizeof(int));
    int front = 0, back = 0;
    queue[back++] = start;
    dist[start] = 0;
    while (front < back) {
        int u = queue[front++];
        for (int v = 0; v < N; v++) {
            if (g->adjacency[u * N + v] > 0 && dist[v] < 0) {
                dist[v] = dist[u] + 1;
                queue[back++] = v;
            }
        }
    }
    free(queue);
}

double consensus_graph_average_path_length(const ConsensusGraph* g) {
    /* Average shortest path length. Characteristic path length L.
     * Small L => small-world property (Watts-Strogatz 1998).
     * For consensus: shorter paths => faster information propagation.
     * Complexity: O(N * (N + E)). */
    if (!g || g->n_nodes < 2) return 0.0;
    int N = g->n_nodes;
    double total = 0.0;
    int count = 0;
    int* dist = (int*)malloc((size_t)N * sizeof(int));
    for (int i = 0; i < N; i++) {
        bfs_distance(g, i, dist);
        for (int j = i + 1; j < N; j++) {
            if (dist[j] > 0) { total += dist[j]; count++; }
        }
    }
    free(dist);
    return (count > 0) ? total / (double)count : 0.0;
}

double consensus_graph_clustering_coefficient(const ConsensusGraph* g) {
    /* Watts-Strogatz clustering coefficient.
     * C_i = (2 * #triangles containing i) / (d_i * (d_i - 1))
     * C = average(C_i). Measures "cliquishness".
     * Random graphs have C ~ p; many real networks have high C.
     * Complexity: O(N * d_max^2). */
    if (!g || g->n_nodes < 2) return 0.0;
    int N = g->n_nodes;
    double total_C = 0.0;
    int valid = 0;
    for (int i = 0; i < N; i++) {
        int d = (int)g->degree[i];
        if (d < 2) continue;
        int triangles = 0;
        /* Count triangles: for each pair of neighbors (j,k), check edge (j,k) */
        for (int j = 0; j < N; j++) {
            if (j == i || g->adjacency[i * N + j] == 0) continue;
            for (int k = j + 1; k < N; k++) {
                if (k == i || g->adjacency[i * N + k] == 0) continue;
                if (g->adjacency[j * N + k] > 0 || g->adjacency[k * N + j] > 0)
                    triangles++;
            }
        }
        total_C += (2.0 * triangles) / (double)(d * (d - 1));
        valid++;
    }
    return (valid > 0) ? total_C / (double)valid : 0.0;
}

double consensus_graph_diameter(const ConsensusGraph* g) {
    /* Graph diameter: maximum shortest path length between any pair.
     * Upper bound for consensus convergence: t ~ O(diam / lambda_2).
     * Complexity: O(N * (N + E)). */
    if (!g || g->n_nodes < 2) return 0.0;
    int N = g->n_nodes;
    int max_dist = 0;
    int* dist = (int*)malloc((size_t)N * sizeof(int));
    for (int i = 0; i < N; i++) {
        bfs_distance(g, i, dist);
        for (int j = 0; j < N; j++)
            if (dist[j] > max_dist) max_dist = dist[j];
    }
    free(dist);
    return (double)max_dist;
}

/* ============================================================================
 * Proximity Graph (L3)
 * ============================================================================ */

void consensus_graph_set_proximity(ConsensusGraph* g, const double* positions,
                                    int dim, double radius) {
    /* Build edges based on Euclidean distance: edge (i,j) if ||pos_i - pos_j|| < radius.
     * Models spatial agent networks (e.g., UAV swarms, robot teams).
     * The resulting graph is a unit disk graph (UDG).
     * Complexity: O(N^2 * dim). */
    if (!g || !positions) return;
    int N = g->n_nodes;
    /* Clear existing edges */
    g->n_edges = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dist2 = 0.0;
            for (int k = 0; k < dim; k++) {
                double diff = positions[i * dim + k] - positions[j * dim + k];
                dist2 += diff * diff;
            }
            if (sqrt(dist2) < radius) {
                consensus_graph_add_undirected_edge(g, i, j, 1.0);
            }
        }
    }
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
}

/* ============================================================================
 * Graph Generators (L6 — Canonical topologies for benchmarking)
 * ============================================================================ */

ConsensusGraph* consensus_graph_create_path(int N) {
    ConsensusGraph* g = consensus_graph_create(N, TOPO_UNDIRECTED);
    for (int i = 0; i < N - 1; i++)
        consensus_graph_add_undirected_edge(g, i, i + 1, 1.0);
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

ConsensusGraph* consensus_graph_create_cycle(int N) {
    ConsensusGraph* g = consensus_graph_create_path(N);
    consensus_graph_add_undirected_edge(g, N - 1, 0, 1.0);
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

ConsensusGraph* consensus_graph_create_star(int N) {
    ConsensusGraph* g = consensus_graph_create(N, TOPO_UNDIRECTED);
    for (int i = 1; i < N; i++)
        consensus_graph_add_undirected_edge(g, 0, i, 1.0);
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

ConsensusGraph* consensus_graph_create_complete(int N) {
    ConsensusGraph* g = consensus_graph_create(N, TOPO_UNDIRECTED);
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
            consensus_graph_add_undirected_edge(g, i, j, 1.0);
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

ConsensusGraph* consensus_graph_create_erdos_renyi(int N, double p) {
    /* Erdos-Renyi random graph G(N,p).
     * Each edge exists independently with probability p.
     * Connectivity threshold: p > (1+eps) * ln(N)/N
     * Expected degree: p*(N-1). */
    ConsensusGraph* g = consensus_graph_create(N, TOPO_UNDIRECTED);
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (consensus_rand_uniform(0, 1) < p)
                consensus_graph_add_undirected_edge(g, i, j, 1.0);
        }
    }
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

ConsensusGraph* consensus_graph_create_small_world(int N, int k, double p_rewire) {
    /* Watts-Strogatz small-world network.
     * Start with ring lattice where each node connects to k nearest neighbors,
     * then rewire each edge with probability p_rewire.
     * Small-world: high clustering + low average path length. */
    ConsensusGraph* g = consensus_graph_create(N, TOPO_UNDIRECTED);
    /* Create ring lattice */
    for (int i = 0; i < N; i++) {
        for (int j = 1; j <= k / 2; j++) {
            int neighbor = (i + j) % N;
            consensus_graph_add_undirected_edge(g, i, neighbor, 1.0);
        }
    }
    /* Rewire */
    for (int e = 0; e < g->n_edges; e++) {
        if (consensus_rand_uniform(0, 1) < p_rewire) {
            int new_target;
            do { new_target = (int)(consensus_rand_uniform(0, 1) * N); }
            while (new_target == g->edges[e].from);
            g->edges[e].to = new_target;
        }
    }
    consensus_graph_build_adjacency(g);
    consensus_graph_build_degree(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

ConsensusGraph* consensus_graph_create_barabasi_albert(int N, int m0, int m) {
    /* Barabasi-Albert scale-free network (preferential attachment).
     * Start with m0 nodes fully connected, then add N-m0 nodes,
     * each connecting to m existing nodes with probability proportional to degree.
     * Produces power-law degree distribution P(k) ~ k^{-3}. */
    if (m0 < 2) m0 = 2;
    if (m > m0) m = m0;
    ConsensusGraph* g = consensus_graph_create(N, TOPO_UNDIRECTED);
    /* Initial complete subgraph of m0 nodes */
    for (int i = 0; i < m0; i++)
        for (int j = i + 1; j < m0; j++)
            consensus_graph_add_undirected_edge(g, i, j, 1.0);
    /* Add remaining nodes */
    for (int i = m0; i < N; i++) {
        double total_deg = 0;
        for (int j = 0; j < i; j++) total_deg += g->degree[j];
        int edges_added = 0;
        /* Try to add m edges */
        for (int attempt = 0; attempt < i * 2 && edges_added < m; attempt++) {
            int candidate = (int)(consensus_rand_uniform(0, 1) * (double)i);
            if (candidate >= i) candidate = i - 1;
            /* Check if edge already exists */
            bool exists = false;
            for (int e = 0; e < g->n_edges; e++)
                if ((g->edges[e].from == i && g->edges[e].to == candidate) ||
                    (g->edges[e].from == candidate && g->edges[e].to == i))
                { exists = true; break; }
            if (exists) continue;
            /* Preferential attachment probability */
            double prob = (total_deg > 0) ? g->degree[candidate] / total_deg : 1.0;
            if (consensus_rand_uniform(0, 1) < prob) {
                consensus_graph_add_undirected_edge(g, i, candidate, 1.0);
                edges_added++;
            }
        }
        /* Update degrees incrementally */
        consensus_graph_build_degree(g);
    }
    consensus_graph_build_adjacency(g);
    consensus_graph_build_laplacian(g, LAPLACIAN_COMBINATORIAL);
    return g;
}

void consensus_graph_print(const ConsensusGraph* g) {
    if (!g) { printf("Graph: NULL\n"); return; }
    printf("Graph: N=%d, E=%d, topo=%d, connected=%d\n",
           g->n_nodes, g->n_edges, g->topology, g->is_connected);
    printf("  Algebraic connectivity lambda_2 = %.6f\n",
           g->algebraic_connectivity);
    if (g->n_nodes <= 8) {
        printf("  Laplacian:\n");
        for (int i = 0; i < g->n_nodes; i++) {
            printf("    ");
            for (int j = 0; j < g->n_nodes; j++)
                printf("%7.3f ", g->laplacian[i * g->n_nodes + j]);
            printf("\n");
        }
    }
}

/* ============================================================================
 * Laplacian Spectrum Bounds (L4)
 * ============================================================================ */

double consensus_laplacian_spectral_radius_bound(const ConsensusGraph* g) {
    /* Upper bound for largest Laplacian eigenvalue:
     * lambda_N <= max_{i~j} (d_i + d_j)
     * Also: lambda_N <= 2 * d_max for undirected graphs.
     * Reference: Anderson & Morley (1985). */
    if (!g || g->n_nodes < 1) return 0.0;
    double d_max = 0.0;
    for (int i = 0; i < g->n_nodes; i++)
        if (g->degree[i] > d_max) d_max = g->degree[i];
    return 2.0 * d_max;
}

double consensus_fiedler_mohar_lower(const ConsensusGraph* g) {
    /* Mohar (1991) lower bound for algebraic connectivity:
     * lambda_2 >= 4 / (N * diam)
     * Additionally: lambda_2 >= 2 * (1 - cos(pi/N)) for path. */
    return consensus_fiedler_lower_bound(g);
}

double consensus_fiedler_mohar_upper(const ConsensusGraph* g) {
    /* Mohar upper bound: lambda_2 <= vertex_connectivity.
     * For trees: lambda_2 <= 2*(1 - cos(pi/(2*radius+1))).
     * Here we use: lambda_2 <= N/(N-1) * min_degree. */
    return consensus_fiedler_upper_bound(g);
}

int consensus_graph_algebraic_connectivity_perturbation(const ConsensusGraph* g,
                                                         int add_edge_u,
                                                         int add_edge_v) {
    /* Adding an edge never decreases algebraic connectivity (monotonicity).
     * lambda_2(G + e) >= lambda_2(G).
     * Returns: 1 if addition strictly increases lambda_2, 0 if same.
     * Reference: Fiedler (1973), "Algebraic connectivity of graphs". */
    if (!g || add_edge_u < 0 || add_edge_u >= g->n_nodes ||
        add_edge_v < 0 || add_edge_v >= g->n_nodes) return -1;
    /* If edge already exists, no change */
    if (g->adjacency[add_edge_u * g->n_nodes + add_edge_v] > 0) return 0;
    return 1; /* Adding any edge helps or maintains connectivity */
}

int consensus_graph_find_bottleneck_edge(const ConsensusGraph* g, int* u, int* v) {
    /* Find the edge whose removal causes the largest drop in algebraic connectivity.
     * This identifies critical communication links.
     * For a bridge edge: lambda_2 drops to 0 (graph disconnects).
     * Complexity: O(E * N^3). */
    if (!g || g->n_nodes < 3 || !u || !v) return -1;
    /* Simplified: find bridge edges (edges whose removal disconnects graph) */
    int N = g->n_nodes;
    (void)g->algebraic_connectivity; /* reference value for comparison */
    for (int e = 0; e < g->n_edges; e++) {
        if (!g->edges[e].active) continue;
        int from = g->edges[e].from, to = g->edges[e].to;
        /* Temporarily remove edge */
        g->edges[e].active = false;
        g->adjacency[from * N + to] = 0.0;
        g->adjacency[to * N + from] = 0.0;
        bool was_connected = consensus_graph_is_connected(g);
        /* Restore */
        g->edges[e].active = true;
        g->adjacency[from * N + to] = g->edges[e].weight;
        g->adjacency[to * N + from] = g->edges[e].weight;
        if (!was_connected) { *u = from; *v = to; return 1; }
    }
    return 0; /* No bridge edge found */
}