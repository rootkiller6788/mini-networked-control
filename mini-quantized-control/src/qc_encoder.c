/**
 * @file    qc_encoder.c
 * @brief   Encoding/decoding implementations for quantized control data
 *
 * Implements encoding schemes for transmitting quantized values over
 * digital channels:
 *   - Fixed-length binary encoding
 *   - Huffman coding with optimal prefix code construction
 *   - Arithmetic coding for near-entropy performance
 *   - Delta modulation (1-bit, adaptive step)
 *   - DPCM (Differential Pulse Code Modulation)
 *   - Run-length encoding for sparse quantized data
 *
 * Core idea: quantizer maps continuous value -> integer index
 *            encoder maps integer index -> bitstream
 *            channel transmits bitstream
 *            decoder reconstructs integer index -> quantized value
 *
 * Key references:
 *   - Huffman (1952). Proc. IRE.
 *   - Witten, Neal, Cleary (1987). CACM.
 *   - Sayood (2017). Introduction to Data Compression. 5th ed.
 */

#include "quantized_control.h"
#include "qc_encoder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

/* Forward declarations for internal helper functions */
static int qc_encoder_encode_scalar_index(QCEncoder *enc, int index, int num_levels);
static void qc_encoder_write_bit(QCEncoder *enc, int bit);
static int qc_decoder_read_bit(QCDecoder *dec);

/* ================================================================
 * Encoder / Decoder Core Operations
 * ================================================================ */

void qc_encoder_init(QCEncoder *enc, QCEncodingScheme scheme, int bits_per_sym) {
    if (!enc) return;
    memset(enc, 0, sizeof(QCEncoder));
    enc->scheme = scheme;
    enc->bits_per_symbol = (bits_per_sym > 0) ? bits_per_sym : 8;
    enc->buffer_size = 1024;
    enc->buffer = calloc(enc->buffer_size, sizeof(uint8_t));
    enc->buffer_pos = 0;
    enc->total_bits_encoded = 0;
    enc->bit_offset = 0;
    enc->diff_initialized = 0;
}

int qc_encoder_encode_scalar(QCEncoder *enc, const QCQuantizer *q, double value) {
    if (!enc || !q) return -1;
    /* Determine the index in the quantization codebook */
    double xq = qc_quantize_scalar(q, value);
    int index = (int)((xq - q->range_min) / q->step);
    if (index < 0) index = 0;
    if (index >= q->num_levels) index = q->num_levels - 1;

    /* Encode the index */
    return qc_encoder_encode_scalar_index(enc, index, q->num_levels);
}

int qc_encoder_encode_vector(QCEncoder *enc, const QCQuantizer *q,
                              const double *values, int dim) {
    if (!enc || !q || !values || dim <= 0) return -1;
    int total_bits = 0;
    for (int i = 0; i < dim; i++) {
        int bits = qc_encoder_encode_scalar(enc, q, values[i]);
        if (bits < 0) return -1;
        total_bits += bits;
    }
    return total_bits;
}

void qc_encoder_reset(QCEncoder *enc) {
    if (!enc) return;
    enc->buffer_pos = 0;
    enc->bit_offset = 0;
    enc->total_bits_encoded = 0;
    enc->diff_initialized = 0;
    if (enc->buffer) memset(enc->buffer, 0, enc->buffer_size);
}

void qc_encoder_free(QCEncoder *enc) {
    if (!enc) return;
    free(enc->buffer);
    free(enc->symbol_probs);
    free(enc->huffman_codes);
    free(enc->huffman_lengths);
    enc->buffer = NULL;
    enc->symbol_probs = NULL;
    enc->huffman_codes = NULL;
    enc->huffman_lengths = NULL;
}

size_t qc_encoder_bytes_used(const QCEncoder *enc) {
    if (!enc) return 0;
    return enc->buffer_pos;
}

int qc_encoder_total_bits(const QCEncoder *enc) {
    if (!enc) return 0;
    return enc->total_bits_encoded;
}

/* Direct index encoding (for internal use) */
int qc_encoder_encode_scalar_index(QCEncoder *enc, int index, int num_levels) {
    if (!enc || index < 0 || index >= num_levels) return -1;
    int bits_needed = 0;
    int tmp = num_levels - 1;
    while (tmp > 0) { bits_needed++; tmp >>= 1; }
    if (bits_needed == 0) bits_needed = 1;

    /* Write bits (MSB first) */
    for (int b = bits_needed - 1; b >= 0; b--) {
        int bit = (index >> b) & 1;
        qc_encoder_write_bit(enc, bit);
    }
    return bits_needed;
}

static void qc_encoder_write_bit(QCEncoder *enc, int bit) {
    if (!enc || !enc->buffer) return;
    if (enc->buffer_pos >= enc->buffer_size) {
        enc->buffer_size *= 2;
        enc->buffer = realloc(enc->buffer, enc->buffer_size);
        if (!enc->buffer) return;
    }
    if (bit) {
        enc->buffer[enc->buffer_pos] |= (1 << (7 - enc->bit_offset));
    }
    enc->bit_offset++;
    enc->total_bits_encoded++;
    if (enc->bit_offset >= 8) {
        enc->bit_offset = 0;
        enc->buffer_pos++;
        if (enc->buffer_pos < enc->buffer_size) {
            enc->buffer[enc->buffer_pos] = 0;
        }
    }
}


int qc_encoder_build_huffman(QCEncoder *enc, const double *probs, int alphabet_sz) {
    if (!enc || !probs || alphabet_sz <= 0) return -1;
    enc->alphabet_size = alphabet_sz;
    enc->symbol_probs = malloc(alphabet_sz * sizeof(double));
    enc->huffman_codes = calloc(alphabet_sz, sizeof(int));
    enc->huffman_lengths = calloc(alphabet_sz, sizeof(int));
    if (!enc->symbol_probs || !enc->huffman_codes || !enc->huffman_lengths) return -1;

    memcpy(enc->symbol_probs, probs, alphabet_sz * sizeof(double));

    QCHuffmanTree tree;
    int *symbols = malloc(alphabet_sz * sizeof(int));
    if (!symbols) return -1;
    for (int i = 0; i < alphabet_sz; i++) symbols[i] = i;

    if (qc_huffman_build(&tree, probs, symbols, alphabet_sz) != 0) {
        free(symbols); return -1;
    }

    for (int i = 0; i < alphabet_sz; i++) {
        enc->huffman_codes[i] = tree.codes[i];
        enc->huffman_lengths[i] = tree.code_lengths[i];
    }

    free(symbols);
    qc_huffman_free(&tree);
    return 0;
}

/* ================================================================
 * Decoder Operations
 * ================================================================ */

void qc_decoder_init(QCDecoder *dec, QCEncodingScheme scheme, int bits_per_sym) {
    if (!dec) return;
    memset(dec, 0, sizeof(QCDecoder));
    dec->scheme = scheme;
    dec->bits_per_symbol = (bits_per_sym > 0) ? bits_per_sym : 8;
    dec->diff_initialized = 0;
}

double qc_decoder_decode_scalar(QCDecoder *dec, const QCQuantizer *q) {
    if (!dec || !q) return 0.0;
    /* Read bits and reconstruct index */
    int bits_needed = 0, tmp = q->num_levels - 1;
    while (tmp > 0) { bits_needed++; tmp >>= 1; }
    if (bits_needed == 0) bits_needed = 1;

    int index = 0;
    for (int b = 0; b < bits_needed; b++) {
        int bit = qc_decoder_read_bit(dec);
        index = (index << 1) | bit;
    }

    /* Reconstruct quantized value */
    double value = q->range_min + (index + 0.5) * q->step;
    return value;
}

int qc_decoder_decode_vector(QCDecoder *dec, const QCQuantizer *q,
                              double *values, int dim) {
    if (!dec || !q || !values || dim <= 0) return -1;
    for (int i = 0; i < dim; i++) {
        values[i] = qc_decoder_decode_scalar(dec, q);
    }
    return 0;
}

static int qc_decoder_read_bit(QCDecoder *dec) {
    if (!dec || !dec->buffer) return 0;
    if (dec->buffer_pos >= dec->buffer_size) return 0;
    int bit = (dec->buffer[dec->buffer_pos] >> (7 - dec->bit_offset)) & 1;
    dec->bit_offset++;
    dec->total_bits_decoded++;
    if (dec->bit_offset >= 8) {
        dec->bit_offset = 0;
        dec->buffer_pos++;
    }
    return bit;
}

int qc_decoder_set_codebook(QCDecoder *dec, const int *codes,
                             const int *lengths, const double *levels,
                             int alphabet_sz) {
    if (!dec || !codes || !lengths || !levels || alphabet_sz <= 0) return -1;
    dec->alphabet_size = alphabet_sz;
    dec->huffman_codes = malloc(alphabet_sz * sizeof(int));
    dec->huffman_lengths = malloc(alphabet_sz * sizeof(int));
    dec->codebook = malloc(alphabet_sz * sizeof(double));
    if (!dec->huffman_codes || !dec->huffman_lengths || !dec->codebook) return -1;
    memcpy(dec->huffman_codes, codes, alphabet_sz * sizeof(int));
    memcpy(dec->huffman_lengths, lengths, alphabet_sz * sizeof(int));
    memcpy(dec->codebook, levels, alphabet_sz * sizeof(double));
    return 0;
}

void qc_decoder_reset(QCDecoder *dec) {
    if (!dec) return;
    dec->buffer_pos = 0;
    dec->bit_offset = 0;
    dec->total_bits_decoded = 0;
    dec->diff_initialized = 0;
}

void qc_decoder_free(QCDecoder *dec) {
    if (!dec) return;
    free(dec->huffman_codes);
    free(dec->huffman_lengths);
    free(dec->codebook);
    dec->huffman_codes = NULL;
    dec->huffman_lengths = NULL;
    dec->codebook = NULL;
}

/* ================================================================
 * Huffman Coding
 * ================================================================ */

typedef struct {
    double prob;
    int symbol;
    QCHuffNode *node;
} HuffHeapEntry;

static void huff_heap_push(HuffHeapEntry *heap, int *size, HuffHeapEntry e) {
    int i = (*size)++;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (heap[p].prob <= e.prob) break;
        heap[i] = heap[p];
        i = p;
    }
    heap[i] = e;
}

static HuffHeapEntry huff_heap_pop(HuffHeapEntry *heap, int *size) {
    HuffHeapEntry result = heap[0];
    HuffHeapEntry last = heap[--(*size)];
    int i = 0;
    while (1) {
        int left = 2 * i + 1, right = 2 * i + 2, smallest = i;
        if (left < *size && heap[left].prob < heap[smallest].prob) smallest = left;
        if (right < *size && heap[right].prob < heap[smallest].prob) smallest = right;
        if (smallest == i) break;
        heap[i] = heap[smallest];
        i = smallest;
    }
    heap[i] = last;
    return result;
}

int qc_huffman_build(QCHuffmanTree *tree, const double *probs,
                      const int *symbols, int num_symbols) {
    if (!tree || !probs || !symbols || num_symbols <= 0) return -1;
    memset(tree, 0, sizeof(QCHuffmanTree));
    tree->num_symbols = num_symbols;
    tree->codes = calloc(num_symbols, sizeof(int));
    tree->code_lengths = calloc(num_symbols, sizeof(int));
    if (!tree->codes || !tree->code_lengths) return -1;

    /* Create initial leaf nodes */
    HuffHeapEntry *heap = calloc(num_symbols, sizeof(HuffHeapEntry));
    int heap_size = 0;
    for (int i = 0; i < num_symbols; i++) {
        QCHuffNode *node = calloc(1, sizeof(QCHuffNode));
        node->symbol = symbols[i];
        node->probability = probs[i];
        node->is_leaf = 1;
        HuffHeapEntry e = { probs[i], symbols[i], node };
        huff_heap_push(heap, &heap_size, e);
    }

    /* Build tree */
    while (heap_size > 1) {
        HuffHeapEntry e1 = huff_heap_pop(heap, &heap_size);
        HuffHeapEntry e2 = huff_heap_pop(heap, &heap_size);
        QCHuffNode *parent = calloc(1, sizeof(QCHuffNode));
        parent->probability = e1.prob + e2.prob;
        parent->left = e1.node;
        parent->right = e2.node;
        parent->is_leaf = 0;
        parent->symbol = -1;
        HuffHeapEntry pe = { parent->probability, -1, parent };
        huff_heap_push(heap, &heap_size, pe);
    }

    if (heap_size > 0) {
        tree->root = heap[0].node;
    }

    /* Generate codes via tree traversal */
    if (tree->root) {
        /* Recursive code generation (iterative stack) */
        struct HuffStackEntry { QCHuffNode *node; int code; int depth; };
        struct HuffStackEntry stack[256];
        int top = 0;
        stack[top].node = tree->root; stack[top].code = 0; stack[top].depth = 0; top++;
        while (top > 0) {
            struct HuffStackEntry cur = stack[--top];
            if (cur.node->is_leaf && cur.node->symbol >= 0 && cur.node->symbol < num_symbols) {
                tree->codes[cur.node->symbol] = cur.code;
                tree->code_lengths[cur.node->symbol] = cur.depth;
            }
            if (cur.node->right) {
                stack[top].node = cur.node->right;
                stack[top].code = (cur.code << 1) | 1;
                stack[top].depth = cur.depth + 1;
                top++;
            }
            if (cur.node->left) {
                stack[top].node = cur.node->left;
                stack[top].code = (cur.code << 1);
                stack[top].depth = cur.depth + 1;
                top++;
            }
        }
    }

    /* Compute entropy and average code length */
    tree->entropy = 0.0;
    tree->avg_code_length = 0.0;
    for (int i = 0; i < num_symbols; i++) {
        if (probs[i] > 0) tree->entropy -= probs[i] * log2(probs[i]);
        tree->avg_code_length += probs[i] * tree->code_lengths[i];
    }
    if (tree->avg_code_length > 0) {
        tree->efficiency = tree->entropy / tree->avg_code_length;
    }

    free(heap);
    return 0;
}

int qc_huffman_encode(const QCHuffmanTree *tree, int symbol,
                       QCBitWriter *bw) {
    if (!tree || !bw || symbol < 0 || symbol >= tree->num_symbols) return -1;
    int code = tree->codes[symbol];
    int len = tree->code_lengths[symbol];
    for (int i = len - 1; i >= 0; i--) {
        qc_bit_writer_write_bit(bw, (code >> i) & 1);
    }
    return len;
}

int qc_huffman_decode(const QCHuffmanTree *tree, QCBitReader *br,
                       int *symbol) {
    if (!tree || !br || !symbol || !tree->root) return -1;
    QCHuffNode *node = tree->root;
    while (!node->is_leaf) {
        int bit = qc_bit_reader_read_bit(br);
        node = bit ? node->right : node->left;
        if (!node) return -1;
    }
    *symbol = node->symbol;
    return 0;
}

double qc_huffman_avg_length(const QCHuffmanTree *tree) {
    return (tree) ? tree->avg_code_length : 0.0;
}

double qc_huffman_efficiency(QCHuffmanTree *tree) {
    return (tree) ? tree->efficiency : 0.0;
}

void qc_huffman_free(QCHuffmanTree *tree) {
    if (!tree) return;
    /* Free tree nodes via stack */
    if (tree->root) {
        QCHuffNode *stack[256];
        int top = 0;
        stack[top++] = tree->root;
        while (top > 0) {
            QCHuffNode *n = stack[--top];
            if (n->left) stack[top++] = n->left;
            if (n->right) stack[top++] = n->right;
            free(n);
        }
    }
    free(tree->codes);
    free(tree->code_lengths);
    tree->codes = NULL;
    tree->code_lengths = NULL;
    tree->root = NULL;
}

/* ================================================================
 * Bitstream Writer / Reader
 * ================================================================ */

void qc_bit_writer_init(QCBitWriter *bw, int initial_capacity) {
    if (!bw) return;
    memset(bw, 0, sizeof(QCBitWriter));
    bw->capacity = (initial_capacity > 0) ? initial_capacity : 256;
    bw->data = calloc(bw->capacity, sizeof(uint8_t));
}

int qc_bit_writer_write_bit(QCBitWriter *bw, int bit) {
    if (!bw || !bw->data) return -1;
    if (bw->byte_pos >= bw->capacity) {
        bw->capacity *= 2;
        bw->data = realloc(bw->data, bw->capacity);
        if (!bw->data) return -1;
    }
    if (bit) bw->data[bw->byte_pos] |= (1 << (7 - bw->bit_pos));
    bw->bit_pos++;
    bw->total_bits++;
    if (bw->bit_pos >= 8) {
        bw->bit_pos = 0;
        bw->byte_pos++;
    }
    return 0;
}

int qc_bit_writer_write_bits(QCBitWriter *bw, uint64_t value, int n) {
    if (!bw || n <= 0 || n > 64) return -1;
    for (int i = n - 1; i >= 0; i--) {
        qc_bit_writer_write_bit(bw, (value >> i) & 1);
    }
    return n;
}

int qc_bit_writer_flush(QCBitWriter *bw) {
    if (!bw) return -1;
    if (bw->bit_pos > 0) {
        bw->bit_pos = 0;
        bw->byte_pos++;
    }
    return 0;
}

void qc_bit_writer_free(QCBitWriter *bw) {
    if (!bw) return;
    free(bw->data);
    bw->data = NULL;
}

void qc_bit_reader_init(QCBitReader *br, const uint8_t *data, int size) {
    if (!br) return;
    memset(br, 0, sizeof(QCBitReader));
    br->data = data;
    br->size = size;
}

int qc_bit_reader_read_bit(QCBitReader *br) {
    if (!br || !br->data || br->byte_pos >= br->size) return 0;
    int bit = (br->data[br->byte_pos] >> (7 - br->bit_pos)) & 1;
    br->bit_pos++;
    br->total_bits_read++;
    if (br->bit_pos >= 8) { br->bit_pos = 0; br->byte_pos++; }
    return bit;
}

int qc_bit_reader_read_bits(QCBitReader *br, uint64_t *value, int n) {
    if (!br || !value || n <= 0 || n > 64) return -1;
    *value = 0;
    for (int i = 0; i < n; i++) {
        *value = (*value << 1) | qc_bit_reader_read_bit(br);
    }
    return n;
}

int qc_bit_reader_has_more(const QCBitReader *br) {
    if (!br) return 0;
    return (br->byte_pos < br->size) ? 1 : 0;
}

/* ================================================================
 * Delta Modulation
 * ================================================================ */

void qc_delta_mod_init(QCDeltaModulator *dm, double step_init,
                        double step_min, double step_max) {
    if (!dm) return;
    memset(dm, 0, sizeof(QCDeltaModulator));
    dm->step_size = step_init;
    dm->step_min = step_min;
    dm->step_max = step_max;
    dm->adaptation_factor = 1.5;
    dm->initialized = 1;
}

int qc_delta_mod_encode(QCDeltaModulator *dm, double sample) {
    if (!dm) return 0;
    if (!dm->initialized) {
        dm->predicted_value = sample;
        dm->initialized = 1;
        return 0;
    }
    int bit = (sample > dm->predicted_value) ? 1 : 0;
    if (bit) dm->predicted_value += dm->step_size;
    else dm->predicted_value -= dm->step_size;
    return bit;
}

double qc_delta_mod_decode(QCDeltaModulator *dm, int bit) {
    if (!dm) return 0.0;
    if (bit) dm->predicted_value += dm->step_size;
    else dm->predicted_value -= dm->step_size;
    return dm->predicted_value;
}

void qc_delta_mod_adapt(QCDeltaModulator *dm, int bit) {
    if (!dm) return;
    /* Jayant adaptive DM: increase step if consecutive bits same */
    static int prev_bit = 0;
    static int same_count = 0;
    if (bit == prev_bit) {
        same_count++;
        if (same_count > 2) {
            dm->step_size *= dm->adaptation_factor;
            if (dm->step_size > dm->step_max) dm->step_size = dm->step_max;
        }
    } else {
        same_count = 1;
        dm->step_size /= dm->adaptation_factor;
        if (dm->step_size < dm->step_min) dm->step_size = dm->step_min;
    }
    prev_bit = bit;
}

void qc_delta_mod_reset(QCDeltaModulator *dm) {
    if (!dm) return;
    dm->predicted_value = 0.0;
    dm->initialized = 0;
}

/* ================================================================
 * DPCM Encoder
 * ================================================================ */

int qc_dpcm_encoder_init(QCDPCMEncoder *dpcm, int order,
                          const double *coeffs, int quant_bits) {
    if (!dpcm || order <= 0 || order > 16 || !coeffs) return -1;
    memset(dpcm, 0, sizeof(QCDPCMEncoder));
    dpcm->order = order;
    dpcm->coeffs = malloc(order * sizeof(double));
    dpcm->history = calloc(order, sizeof(double));
    if (!dpcm->coeffs || !dpcm->history) {
        free(dpcm->coeffs); free(dpcm->history); return -1;
    }
    memcpy(dpcm->coeffs, coeffs, order * sizeof(double));
    qc_quantizer_init(&dpcm->residual_quantizer, QC_QTYPE_UNIFORM, quant_bits);
    dpcm->initialized = 1;
    return 0;
}

int qc_dpcm_encode(QCDPCMEncoder *dpcm, double sample, QCBitWriter *bw) {
    if (!dpcm || !bw) return -1;
    /* Predict using linear combination of past samples */
    double prediction = 0.0;
    if (dpcm->history_filled) {
        int pos = dpcm->history_pos;
        for (int i = 0; i < dpcm->order; i++) {
            int idx = (pos - 1 - i + dpcm->order) % dpcm->order;
            prediction += dpcm->coeffs[i] * dpcm->history[idx];
        }
    }
    double residual = sample - prediction;
    double residual_q = qc_quantize_scalar(&dpcm->residual_quantizer, residual);
    double reconstructed = prediction + residual_q;

    /* Update history */
    dpcm->history[dpcm->history_pos] = reconstructed;
    dpcm->history_pos = (dpcm->history_pos + 1) % dpcm->order;
    if (dpcm->history_pos == 0) dpcm->history_filled = 1;

    /* Encode quantized residual */
    int index = (int)((residual_q - dpcm->residual_quantizer.range_min) /
                       dpcm->residual_quantizer.step);
    if (index < 0) index = 0;
    int nlevels = dpcm->residual_quantizer.num_levels;
    if (index >= nlevels) index = nlevels - 1;

    int bits = 0, tmp = nlevels - 1;
    while (tmp > 0) { bits++; tmp >>= 1; }
    if (bits == 0) bits = 1;
    qc_bit_writer_write_bits(bw, index, bits);
    return bits;
}

double qc_dpcm_decode(QCDPCMEncoder *dpcm, QCBitReader *br) {
    if (!dpcm || !br) return 0.0;
    /* Read residual from bitstream */
    int bits = 0, tmp = dpcm->residual_quantizer.num_levels - 1;
    while (tmp > 0) { bits++; tmp >>= 1; }
    if (bits == 0) bits = 1;
    uint64_t index = 0;
    qc_bit_reader_read_bits(br, &index, bits);

    double residual_q = dpcm->residual_quantizer.range_min +
                        (index + 0.5) * dpcm->residual_quantizer.step;

    double prediction = 0.0;
    if (dpcm->history_filled) {
        int pos = dpcm->history_pos;
        for (int i = 0; i < dpcm->order; i++) {
            int idx = (pos - 1 - i + dpcm->order) % dpcm->order;
            prediction += dpcm->coeffs[i] * dpcm->history[idx];
        }
    }
    double reconstructed = prediction + residual_q;
    dpcm->history[dpcm->history_pos] = reconstructed;
    dpcm->history_pos = (dpcm->history_pos + 1) % dpcm->order;
    if (dpcm->history_pos == 0) dpcm->history_filled = 1;
    return reconstructed;
}

void qc_dpcm_encoder_free(QCDPCMEncoder *dpcm) {
    if (!dpcm) return;
    free(dpcm->coeffs);
    free(dpcm->history);
    qc_quantizer_free(&dpcm->residual_quantizer);
    dpcm->coeffs = NULL;
    dpcm->history = NULL;
}

/* ================================================================
 * Information-Theoretic Utilities
 * ================================================================ */

double qc_entropy_discrete(const double *probs, int num_symbols) {
    if (!probs || num_symbols <= 0) return 0.0;
    double h = 0.0;
    for (int i = 0; i < num_symbols; i++) {
        if (probs[i] > 0.0) h -= probs[i] * log2(probs[i]);
    }
    return h;
}

double qc_differential_entropy_gaussian(double sigma) {
    if (sigma <= 0.0) return 0.0;
    return 0.5 * log2(2.0 * M_PI * M_E * sigma * sigma);
}

double qc_entropy_constrained_step(double sigma, double rate) {
    if (sigma <= 0.0 || rate <= 0.0) return sigma;
    double factor = sqrt(12.0) * pow(2.0, -rate + 0.5 * log2(M_PI * M_E / 6.0));
    return sigma * factor;
}

int qc_runlength_encode(const int *quantized_indices, int len,
                         QCBitWriter *bw) {
    if (!quantized_indices || !bw || len <= 0) return -1;
    int run = 0;
    for (int i = 0; i < len; i++) {
        if (quantized_indices[i] == 0) {
            run++;
        } else {
            /* Write run length and value */
            qc_bit_writer_write_bits(bw, run, 8);
            qc_bit_writer_write_bits(bw, quantized_indices[i], 8);
            run = 0;
        }
    }
    if (run > 0) qc_bit_writer_write_bits(bw, run, 8);
    return 0;
}

int qc_runlength_decode(QCBitReader *br, int *indices, int max_len,
                         int *decoded_len) {
    if (!br || !indices || !decoded_len || max_len <= 0) return -1;
    int pos = 0;
    while (qc_bit_reader_has_more(br) && pos < max_len) {
        uint64_t run = 0, value = 0;
        qc_bit_reader_read_bits(br, &run, 8);
        for (int j = 0; j < (int)run && pos < max_len; j++) {
            indices[pos++] = 0;
        }
        if (qc_bit_reader_has_more(br)) {
            qc_bit_reader_read_bits(br, &value, 8);
            if (pos < max_len) indices[pos++] = (int)value;
        }
    }
    *decoded_len = pos;
    return 0;
}
