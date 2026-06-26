/**
 * blc_encoding.c — Encoding Strategies for Bandwidth-Limited Control
 *
 * Implementation of encoding and quantization strategies:
 *  - Logarithmic quantization (optimal for stabilization)
 *  - Lloyd-Max optimal quantization (minimum MSE for given distribution)
 *  - Vector quantization via LBG algorithm
 *  - Huffman entropy coding
 *  - Arithmetic coding
 *  - Run-length encoding for sparse control signals
 *
 * Key results:
 *  - Elia & Mitter (2001): logarithmic quantization is optimal
 *    for stabilizing linear systems under data rate constraints
 *  - Lloyd (1957), Max (1960): iterative algorithm for optimal
 *    scalar quantizer design
 *  - Linde, Buzo, Gray (1980): LBG algorithm for vector quantization
 *
 * Knowledge coverage: L5 (Algorithms), L7 (Applications)
 */

#include "blc_encoding.h"
#include "blc_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * Logarithmic Quantizer
 *
 * For stabilizing an unstable scalar system with pole λ > 0,
 * the coarsest stabilizing quantizer has density δ = e^{λ·T_s}.
 *
 * The quantization levels are: q_k = ρ^k · q_0, k = 0, ±1, ±2, ...
 * where ρ = δ is the optimal density.
 *
 * A logarithmic quantizer with N positive levels can stabilize
 * any system with:
 *   λ_max ≤ (1/T_s) · ln(ρ)  where ρ = (1 + Δ)/(1 - Δ), Δ = 1/(2N)
 *
 * @ref Elia & Mitter (2001), "Stabilization of linear systems with
 *      limited information", IEEE TAC
 * ================================================================ */

int blc_log_quant_init(BLCLogQuantizer* lq, double density,
                        int num_pos, double deadzone) {
    if (!lq || density <= 1.0 || num_pos < 1 || deadzone <= 0.0) return -1;

    lq->num_levels_pos = num_pos;
    lq->num_levels_neg = num_pos;
    lq->density        = density;
    lq->deadzone        = deadzone;

    int total_levels = num_pos * 2 + 1;  /** +1 for the deadzone (0) */
    lq->level_table = (double*)calloc((size_t)total_levels, sizeof(double));
    lq->recon_table = (double*)calloc((size_t)total_levels, sizeof(double));

    if (!lq->level_table || !lq->recon_table) {
        free(lq->level_table);
        free(lq->recon_table);
        return -2;
    }

    /** Build symmetric logarithmic levels:
     *  Level 0: deadzone [-d, d]
     *  Level +k: [d·ρ^{k-1}, d·ρ^k]  for k = 1..num_pos
     *  Level -k: [-d·ρ^k, -d·ρ^{k-1}] */

    /** Decision boundaries */
    double pos_boundaries[128];  /** num_pos + 1 entries */
    pos_boundaries[0] = deadzone;
    for (int k = 1; k <= num_pos; k++) {
        pos_boundaries[k] = deadzone * pow(density, (double)k);
    }
    lq->overload_value = pos_boundaries[num_pos];

    /** Reconstruction values: centroid of each interval.
     *  For logarithmic quantizer: r_k = (b_k · b_{k-1})^{1/2} ≈ geometric mean
     *  or more precisely: r_k = (b_k - b_{k-1}) / ln(b_k / b_{k-1})
     *  For simplicity: geometric mean */
    for (int k = 0; k < num_pos; k++) {
        double lo = (k == 0) ? 0.0 : pos_boundaries[k];
        double hi = pos_boundaries[k+1];
        /** Optimal reconstruction for uniform distribution in log-space:
         *  r_k = (hi - lo) / ln(hi/lo)  → approximates the centroid */
        if (lo < 1e-15) {
            lq->recon_table[num_pos + k] = hi / density;
        } else {
            lq->recon_table[num_pos + k] = (hi - lo) / log(hi / lo);
        }
        lq->recon_table[num_pos - 1 - k] = -lq->recon_table[num_pos + k];
    }
    /** Deadzone reconstruction: 0 */
    lq->recon_table[num_pos] = 0.0;

    /** Level boundaries for encoding */
    for (int k = 0; k < total_levels; k++) {
        lq->level_table[k] = (double)(k - num_pos);  /** Index: -N ... 0 ... +N */
    }

    return 0;
}

void blc_log_quant_free(BLCLogQuantizer* lq) {
    if (!lq) return;
    free(lq->level_table);
    free(lq->recon_table);
    lq->level_table = NULL;
    lq->recon_table = NULL;
}

int blc_log_quant_encode(const BLCLogQuantizer* lq, double value) {
    if (!lq) return 0;
    int N = lq->num_levels_pos;

    /** Determine which interval the value falls into */
    double abs_val = fabs(value);

    if (abs_val <= lq->deadzone) {
        return N;  /** Deadzone → middle index */
    }

    /** Find k such that d·ρ^k < abs_val ≤ d·ρ^{k+1}
     *  k = floor(log(abs_val / deadzone) / log(density)) */
    double k_d = log(abs_val / lq->deadzone) / log(lq->density);
    int k = (int)floor(k_d);

    /** Clamp to valid range */
    if (k >= N) {
        /** Overload: saturate to maximum level */
        return (value > 0) ? (2 * N) : 0;
    }
    if (k < 0) k = 0;

    if (value > 0) {
        return N + k;  /** Positive side: recon_table[N + k] */
    } else {
        return N - 1 - k;  /** Negative side: recon_table[N - 1 - k] */
    }
}

double blc_log_quant_decode(const BLCLogQuantizer* lq, int index) {
    if (!lq) return 0.0;
    int total = lq->num_levels_pos * 2 + 1;
    if (index < 0 || index >= total) return 0.0;
    return lq->recon_table[index];
}

double blc_log_quant_optimal_density(double lambda_max, double Ts) {
    /** δ_opt = e^{λ_max · T_s}
     *  This is the minimum density that stabilizes a system with
     *  unstable eigenvalue λ_max when sampled at period T_s. */
    if (lambda_max <= 0.0 || Ts <= 0.0) return 1.1;  /** Minimum practical */
    return exp(lambda_max * Ts);
}

int blc_log_quant_min_levels(double density) {
    /** N_min = ceil((density + 1) / (density - 1) / 2)
     *  This comes from the stability condition:
     *    Δ = 1/(2N+1) < (δ-1)/(δ+1) */
    if (density <= 1.0) return 1;
    double N_d = (density + 1.0) / (density - 1.0) / 2.0;
    return (int)ceil(N_d);
}

/* ================================================================
 * Lloyd-Max Optimal Quantizer
 *
 * Lloyd-Max algorithm: iteratively optimize decision boundaries
 * and reconstruction levels to minimize MSE for a given
 * probability density function.
 *
 * Algorithm:
 *   1. Initialize reconstruction levels {r_i}
 *   2. Set boundaries: b_i = (r_i + r_{i+1}) / 2
 *   3. Update reconstructions: r_i = E[x | x ∈ (b_{i-1}, b_i]]
 *   4. Repeat until convergence
 *
 * For a uniform source over [0, 1] with N levels:
 *   MSE_min = 1 / (12·N²)
 *
 * @ref Lloyd (1957), "Least squares quantization in PCM"
 * @ref Max (1960), "Quantizing for minimum distortion"
 * ================================================================ */

int blc_lloyd_max_design(BLCLloydQuantizer* q, const double* pdf_samples,
                          int n_samples, double range_lo, double range_hi,
                          int max_iter, double tol) {
    if (!q || !pdf_samples || n_samples < 1) return -1;
    int N = q->num_levels;
    if (N < 1) return -2;

    /** Allocate boundaries and reconstruction arrays */
    q->boundaries = (double*)calloc((size_t)(N + 1), sizeof(double));
    q->reconstructions = (double*)calloc((size_t)N, sizeof(double));
    if (!q->boundaries || !q->reconstructions) return -3;

    /** Initialize with uniform spacing */
    double width = (range_hi - range_lo) / (double)N;
    for (int i = 0; i <= N; i++) {
        q->boundaries[i] = range_lo + (double)i * width;
    }
    for (int i = 0; i < N; i++) {
        q->reconstructions[i] = (q->boundaries[i] + q->boundaries[i+1]) / 2.0;
    }

    double prev_mse = DBL_MAX;
    q->iterations = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        /** Step 1: Update reconstruction levels as centroids */
        double* sum_val  = (double*)calloc((size_t)N, sizeof(double));
        int*    count    = (int*)calloc((size_t)N, sizeof(int));
        if (!sum_val || !count) {
            free(sum_val); free(count);
            return -3;
        }

        for (int s = 0; s < n_samples; s++) {
            double x = pdf_samples[s];
            /** Find which interval x belongs to */
            int idx = 0;
            for (int i = 0; i < N; i++) {
                if (x >= q->boundaries[i] && x < q->boundaries[i+1]) {
                    idx = i;
                    break;
                }
            }
            /** Handle last boundary inclusively */
            if (x >= q->boundaries[N]) idx = N - 1;

            sum_val[idx] += x;
            count[idx]++;
        }

        /** Compute new reconstruction levels as centroids */
        double mse = 0.0;
        for (int i = 0; i < N; i++) {
            if (count[i] > 0) {
                q->reconstructions[i] = sum_val[i] / (double)count[i];
            }
            /** Accumulate MSE */
            for (int s = 0; s < n_samples; s++) {
                double x = pdf_samples[s];
                if (x >= q->boundaries[i] && x < q->boundaries[i+1]) {
                    double e = x - q->reconstructions[i];
                    mse += e * e;
                }
            }
        }
        mse /= (double)n_samples;
        free(sum_val);
        free(count);

        /** Step 2: Update boundaries as midpoints */
        for (int i = 1; i < N; i++) {
            q->boundaries[i] = (q->reconstructions[i-1] + q->reconstructions[i]) / 2.0;
        }

        /** Check convergence */
        if (fabs(prev_mse - mse) < tol) {
            q->mse = mse;
            q->iterations = iter + 1;
            return iter + 1;
        }
        prev_mse = mse;
    }

    q->mse = prev_mse;
    q->iterations = max_iter;
    return max_iter;
}

void blc_lloyd_max_free(BLCLloydQuantizer* q) {
    if (!q) return;
    free(q->boundaries);
    free(q->reconstructions);
    q->boundaries = NULL;
    q->reconstructions = NULL;
}

int blc_lloyd_max_encode(const BLCLloydQuantizer* q, double value) {
    if (!q || !q->boundaries) return 0;
    int N = q->num_levels;
    for (int i = 0; i < N; i++) {
        if (value >= q->boundaries[i] && value < q->boundaries[i+1]) {
            return i;
        }
    }
    return N - 1;  /** Saturation */
}

double blc_lloyd_max_decode(const BLCLloydQuantizer* q, int index) {
    if (!q || !q->reconstructions || index < 0 || index >= q->num_levels)
        return 0.0;
    return q->reconstructions[index];
}

double blc_lloyd_max_get_mse(const BLCLloydQuantizer* q) {
    return q ? q->mse : 0.0;
}

/* ================================================================
 * Vector Quantizer (LBG Algorithm)
 *
 * The LBG algorithm (Linde, Buzo, Gray, 1980) is the vector
 * generalization of Lloyd-Max. It partitions the input space
 * into Voronoi regions and computes their centroids.
 *
 * Algorithm:
 *   1. Start with 1 codeword: the centroid of all training data
 *   2. Split: double the codebook by perturbing each codeword
 *   3. Classify: assign each training vector to nearest codeword
 *   4. Update: recompute each codeword as centroid of its cell
 *   5. Repeat 3-4 until convergence
 *   6. Repeat 2-5 until desired codebook size
 * ================================================================ */

int blc_vq_init(BLCVectorQuantizer* vq, int dimension, int codebook_size) {
    if (!vq || dimension < 1 || dimension > BLC_MAX_STATES
        || codebook_size < 1) return -1;

    vq->dimension      = dimension;
    vq->codebook_size  = codebook_size;
    vq->mse            = 0.0;
    vq->bits_per_sample = log2((double)codebook_size) / (double)dimension;

    /** Allocate codebook: codebook_size × dimension */
    vq->codebook = (double**)calloc((size_t)codebook_size, sizeof(double*));
    vq->cell_counts = (int*)calloc((size_t)codebook_size, sizeof(int));
    if (!vq->codebook || !vq->cell_counts) {
        free(vq->codebook);
        free(vq->cell_counts);
        return -2;
    }
    for (int i = 0; i < codebook_size; i++) {
        vq->codebook[i] = (double*)calloc((size_t)dimension, sizeof(double));
        if (!vq->codebook[i]) return -2;
    }

    return 0;
}

void blc_vq_free(BLCVectorQuantizer* vq) {
    if (!vq) return;
    if (vq->codebook) {
        for (int i = 0; i < vq->codebook_size; i++) {
            free(vq->codebook[i]);
        }
        free(vq->codebook);
    }
    free(vq->cell_counts);
    vq->codebook = NULL;
    vq->cell_counts = NULL;
}

int blc_vq_train(BLCVectorQuantizer* vq, const double* training_data,
                  int n_vectors, double epsilon, int max_iter) {
    if (!vq || !training_data || n_vectors < 1) return -1;
    int dim = vq->dimension;
    int K   = vq->codebook_size;

    /** Initialize codebook with first K training vectors */
    for (int k = 0; k < K && k < n_vectors; k++) {
        for (int d = 0; d < dim; d++) {
            vq->codebook[k][d] = training_data[k * dim + d];
        }
    }

    double prev_distortion = DBL_MAX;
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /** Classification: assign each vector to nearest codeword */
        double* sum_vec = (double*)calloc((size_t)(K * dim), sizeof(double));
        int*    counts  = (int*)calloc((size_t)K, sizeof(int));
        if (!sum_vec || !counts) { free(sum_vec); free(counts); return -2; }

        double distortion = 0.0;

        for (int v = 0; v < n_vectors; v++) {
            /** Find nearest codeword */
            double min_dist = DBL_MAX;
            int nearest = 0;
            for (int k = 0; k < K; k++) {
                double dist = 0.0;
                for (int d = 0; d < dim; d++) {
                    double diff = training_data[v * dim + d] - vq->codebook[k][d];
                    dist += diff * diff;
                }
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest = k;
                }
            }
            distortion += min_dist;
            counts[nearest]++;
            for (int d = 0; d < dim; d++) {
                sum_vec[nearest * dim + d] += training_data[v * dim + d];
            }
        }
        distortion /= (double)n_vectors;

        /** Update codewords as centroids */
        for (int k = 0; k < K; k++) {
            if (counts[k] > 0) {
                for (int d = 0; d < dim; d++) {
                    vq->codebook[k][d] = sum_vec[k * dim + d] / (double)counts[k];
                }
            }
            vq->cell_counts[k] = counts[k];
        }

        free(sum_vec);
        free(counts);

        vq->mse = distortion;

        /** Check convergence */
        if (fabs(prev_distortion - distortion) / (prev_distortion + 1e-10) < epsilon) {
            break;
        }
        prev_distortion = distortion;
    }

    return iter + 1;
}

int blc_vq_encode(const BLCVectorQuantizer* vq, const double* x) {
    if (!vq || !x) return 0;
    int K = vq->codebook_size;
    int dim = vq->dimension;

    double min_dist = DBL_MAX;
    int nearest = 0;
    for (int k = 0; k < K; k++) {
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = x[d] - vq->codebook[k][d];
            dist += diff * diff;
        }
        if (dist < min_dist) {
            min_dist = dist;
            nearest = k;
        }
    }
    return nearest;
}

void blc_vq_decode(const BLCVectorQuantizer* vq, int index, double* x_hat) {
    if (!vq || !x_hat || index < 0 || index >= vq->codebook_size) return;
    for (int d = 0; d < vq->dimension; d++) {
        x_hat[d] = vq->codebook[index][d];
    }
}

double blc_vq_distortion(const BLCVectorQuantizer* vq,
                          const double* data, int n_vectors) {
    if (!vq || !data || n_vectors < 1) return 0.0;
    int dim = vq->dimension;
    double total = 0.0;

    for (int v = 0; v < n_vectors; v++) {
        int idx = blc_vq_encode(vq, &data[v * dim]);
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
            double diff = data[v * dim + d] - vq->codebook[idx][d];
            dist += diff * diff;
        }
        total += dist;
    }
    return total / (double)n_vectors;
}

int blc_vq_split_codebook(BLCVectorQuantizer* vq, double epsilon) {
    /** Split each codeword into two slightly perturbed copies.
     *  This doubles the codebook size, used in the splitting
     *  phase of the LBG algorithm. */
    if (!vq) return -1;
    int old_size = vq->codebook_size;
    int new_size = old_size * 2;

    if (new_size > 256) return -2;  /** Practical limit */

    double** new_codebook = (double**)calloc((size_t)new_size, sizeof(double*));
    int* new_counts = (int*)calloc((size_t)new_size, sizeof(int));
    if (!new_codebook || !new_counts) { free(new_codebook); free(new_counts); return -2; }

    for (int i = 0; i < new_size; i++) {
        new_codebook[i] = (double*)calloc((size_t)vq->dimension, sizeof(double));
    }

    for (int k = 0; k < old_size; k++) {
        for (int d = 0; d < vq->dimension; d++) {
            new_codebook[2*k][d]   = vq->codebook[k][d] * (1.0 + epsilon);
            new_codebook[2*k+1][d] = vq->codebook[k][d] * (1.0 - epsilon);
        }
    }

    /** Free old codebook */
    for (int i = 0; i < old_size; i++) free(vq->codebook[i]);
    free(vq->codebook);
    free(vq->cell_counts);

    vq->codebook = new_codebook;
    vq->cell_counts = new_counts;
    vq->codebook_size = new_size;
    vq->bits_per_sample = log2((double)new_size) / (double)vq->dimension;

    return 0;
}

double blc_vq_rate(const BLCVectorQuantizer* vq) {
    return vq ? vq->bits_per_sample : 0.0;
}

/* ================================================================
 * Huffman Entropy Coding
 *
 * Huffman algorithm (1952) produces optimal prefix codes
 * for a given probability distribution. Used here to approach
 * the entropy limit H = -Σ p_i log₂(p_i) after quantization.
 *
 * Algorithm:
 *   1. Start with each symbol as a leaf node with probability p_i
 *   2. Repeatedly combine two lowest-probability nodes
 *   3. Assign 0/1 to left/right branches
 *   4. Code length = tree depth of leaf
 * ================================================================ */

/** Internal Huffman tree node */
typedef struct {
    double prob;
    int    symbol;    /** -1 for internal nodes */
    int    left;
    int    right;
    int    parent;
} HuffNode;

int blc_huffman_build(BLCHuffmanCoder* huff, const double* probabilities) {
    if (!huff || !probabilities) return -1;
    int M = huff->alphabet_size;
    if (M < 2 || M > 256) return -2;

    /** Copy probabilities */
    huff->probabilities = (double*)calloc((size_t)M, sizeof(double));
    if (!huff->probabilities) return -3;
    memcpy(huff->probabilities, probabilities, (size_t)M * sizeof(double));

    /** Compute entropy */
    huff->entropy = blc_entropy(probabilities, M);

    /** Allocate tree nodes: 2M-1 nodes for M symbols */
    int total_nodes = 2 * M - 1;
    HuffNode* nodes = (HuffNode*)calloc((size_t)total_nodes, sizeof(HuffNode));
    if (!nodes) { free(huff->probabilities); return -3; }

    /** Initialize leaf nodes */
    int n_leaves = M;
    for (int i = 0; i < M; i++) {
        nodes[i].prob   = probabilities[i];
        nodes[i].symbol = i;
        nodes[i].left   = -1;
        nodes[i].right  = -1;
        nodes[i].parent = -1;
    }

    /** Build Huffman tree */
    int next_internal = M;
    while (n_leaves > 1) {
        /** Find two smallest */
        int min1 = -1, min2 = -1;
        double p1 = 2.0, p2 = 2.0;
        for (int i = 0; i < next_internal; i++) {
            if (nodes[i].parent == -1 && nodes[i].prob < p1) {
                p2 = p1; min2 = min1;
                p1 = nodes[i].prob; min1 = i;
            } else if (nodes[i].parent == -1 && nodes[i].prob < p2) {
                p2 = nodes[i].prob; min2 = i;
            }
        }
        if (min1 < 0 || min2 < 0) break;

        /** Create internal node */
        nodes[next_internal].prob   = p1 + p2;
        nodes[next_internal].symbol = -1;
        nodes[next_internal].left   = min1;
        nodes[next_internal].right  = min2;
        nodes[next_internal].parent = -1;
        nodes[min1].parent = next_internal;
        nodes[min2].parent = next_internal;
        n_leaves--;
        next_internal++;
    }

    /** Extract code lengths from tree depths */
    huff->code_lengths = (int*)calloc((size_t)M, sizeof(int));
    if (!huff->code_lengths) { free(nodes); return -3; }

    double avg_len = 0.0;
    huff->max_code_length = 0;
    for (int i = 0; i < M; i++) {
        int depth = 0;
        int node = i;
        while (nodes[node].parent != -1) {
            depth++;
            node = nodes[node].parent;
        }
        huff->code_lengths[i] = depth;
        if (depth > huff->max_code_length) huff->max_code_length = depth;
        avg_len += probabilities[i] * (double)depth;
    }
    huff->avg_code_length = avg_len;
    huff->is_optimal      = (avg_len <= huff->entropy + 1.0);

    free(nodes);
    return 0;
}

void blc_huffman_free(BLCHuffmanCoder* huff) {
    if (!huff) return;
    free(huff->probabilities);
    free(huff->code_lengths);
    huff->probabilities = NULL;
    huff->code_lengths  = NULL;
}

int blc_huffman_get_length(const BLCHuffmanCoder* huff, int symbol) {
    if (!huff || !huff->code_lengths || symbol < 0
        || symbol >= huff->alphabet_size) return 0;
    return huff->code_lengths[symbol];
}

int blc_huffman_encode_symbol(const BLCHuffmanCoder* huff,
                               int symbol, uint8_t* bits, int* offset) {
    /** Returns the length — actual bit packing would be done
     *  by the caller using the canonical Huffman code. */
    if (!huff) return 0;
    (void)bits; (void)offset;
    return blc_huffman_get_length(huff, symbol);
}

int blc_huffman_decode_symbol(const BLCHuffmanCoder* huff,
                               const uint8_t* bits, int* offset,
                               int total_bits) {
    /** For canonical Huffman: walk the tree using current bit.
     *  Simplified: return first symbol whose code length matches
     *  the pattern. In practice, a lookup table is used. */
    if (!huff || !bits || !offset) return 0;
    /** Placeholder — real implementation would walk a decoding tree */
    (void)total_bits;
    (*offset) += 8;  /** Consume 8 bits */
    return 0;
}

double blc_entropy(const double* probabilities, int n) {
    double H = 0.0;
    for (int i = 0; i < n; i++) {
        if (probabilities[i] > 0.0) {
            H -= probabilities[i] * log2(probabilities[i]);
        }
    }
    return H;
}

double blc_huffman_efficiency(const BLCHuffmanCoder* huff) {
    if (!huff || huff->avg_code_length <= 0.0) return 0.0;
    return huff->entropy / huff->avg_code_length;
}

/* ================================================================
 * Arithmetic Coding
 *
 * Arithmetic coding represents the entire message as a single
 * number in [0, 1). It achieves near-entropy compression,
 * typically within 0.1-0.5% of the entropy bound.
 *
 * Encoding: subdivide [0, 1) based on cumulative probabilities,
 * emit bits as soon as the MSB is fixed.
 *
 * @ref Witten, Neal, Cleary (1987), "Arithmetic coding for data
 *      compression", Communications of the ACM
 * ================================================================ */

void blc_arithmetic_init(BLCArithmeticEncoder* ae) {
    if (!ae) return;
    ae->low  = 0;
    ae->high = 0xFFFFFFFFFFFFFFFFULL;
    ae->range = ae->high - ae->low;
    ae->pending_bits = 0;
    ae->total_bits   = 0;
    memset(ae->symbol_count, 0, sizeof(ae->symbol_count));
}

int blc_arithmetic_encode(BLCArithmeticEncoder* ae,
                           double cum_low, double cum_high, double total) {
    if (!ae || total <= 0.0) return -1;

    uint64_t range = ae->high - ae->low;
    ae->high = ae->low + (uint64_t)((double)range * cum_high / total);
    ae->low  = ae->low + (uint64_t)((double)range * cum_low  / total);

    /** Emit bits as MSB stabilizes */
    while (1) {
        if ((ae->high >> 63) == (ae->low >> 63)) {
            /** MSB is fixed */
            /** In a real implementation, output the bit here */
            ae->total_bits++;
            ae->low  <<= 1;
            ae->high  = (ae->high << 1) | 1;
            /** Output any pending opposite bits */
            ae->pending_bits = 0;
        } else if ((ae->low >> 62) == 2 && (ae->high >> 62) == 1) {
            /** Underflow: rescale */
            ae->pending_bits++;
            ae->low  = (ae->low << 1) & 0x7FFFFFFFFFFFFFFFULL;
            ae->high = ((ae->high << 1) & 0x7FFFFFFFFFFFFFFFULL) | 0x8000000000000001ULL;
        } else {
            break;
        }
    }

    return 0;
}

int blc_arithmetic_finalize(BLCArithmeticEncoder* ae) {
    if (!ae) return 0;
    ae->total_bits += ae->pending_bits + 2;
    return ae->total_bits;
}

int blc_arithmetic_get_bits(const BLCArithmeticEncoder* ae) {
    return ae ? ae->total_bits : 0;
}

/* ================================================================
 * Bandwidth-Optimal Encoding Utilities
 * ================================================================ */

void blc_allocate_bits_eigenvalue(const double* eigenvalues, int n,
                                   int total_bits, int* allocation) {
    /** Allocate bits proportional to |Re(λ)| for unstable eigenvalues.
     *  Stable modes get minimum 1 bit, unstable modes get proportional share. */
    if (!eigenvalues || !allocation || n < 1) return;

    /** Sum positive real parts */
    double sum_pos = 0.0;
    for (int i = 0; i < n; i++) {
        if (eigenvalues[i] > 1e-12) {
            sum_pos += eigenvalues[i];
        }
    }

    if (sum_pos < 1e-12) {
        /** All eigenvalues stable: equal allocation */
        int each = total_bits / n;
        for (int i = 0; i < n; i++) allocation[i] = (each > 0) ? each : 1;
        return;
    }

    int allocated = 0;
    for (int i = 0; i < n; i++) {
        if (eigenvalues[i] > 1e-12) {
            double share = eigenvalues[i] / sum_pos;
            allocation[i] = (int)(share * (double)total_bits);
            if (allocation[i] < 2) allocation[i] = 2;
        } else {
            allocation[i] = 1;  /** 1 bit for stable modes */
        }
        allocated += allocation[i];
    }

    /** Distribute remainder */
    int rem = total_bits - allocated;
    int idx = 0;
    while (rem > 0 && idx < n * 2) {
        if (eigenvalues[idx % n] > 1e-12) {
            allocation[idx % n]++;
            rem--;
        }
        idx++;
    }
}

double blc_encoding_distortion_bound(double variance, double bits_per_sample) {
    /** Shannon lower bound for Gaussian source:
     *  D ≥ σ² · 2^{-2R}
     *  This is the theoretically minimal distortion for rate R. */
    if (bits_per_sample <= 0.0) return variance;
    return variance * pow(2.0, -2.0 * bits_per_sample);
}

double blc_encoding_bitrate_bound(double variance, double target_distortion) {
    /** R ≥ 0.5 · log₂(σ² / D) for D ≤ σ² */
    if (target_distortion <= 0.0 || variance <= 0.0) return 0.0;
    if (target_distortion >= variance) return 0.0;
    return 0.5 * log2(variance / target_distortion);
}

int blc_runlength_encode(const double* signal, int n,
                          double* run_values, int* run_lengths,
                          int max_runs) {
    /** Run-length encode for sparse signals.
     *  Returns number of runs. Each run is (value, count). */
    if (!signal || !run_values || !run_lengths || n < 1) return 0;

    int runs = 0;
    int i = 0;
    while (i < n && runs < max_runs) {
        double current_val = signal[i];
        int count = 1;
        i++;
        while (i < n && signal[i] == current_val && count < 65535) {
            count++;
            i++;
        }
        run_values[runs] = current_val;
        run_lengths[runs] = count;
        runs++;
    }
    return runs;
}