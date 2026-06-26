#include "packet_loss_predictor.h"
#include "packet_loss_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

static double* _m(int r, int c) { return (double*)calloc(r * c, sizeof(double)); }

/* ============================================================================
 * Packetized Predictive Control (PPC)
 *
 * Bemporad (1998), Quevedo & Nesic (2012): Instead of sending a single
 * control u_k, the controller packs H future controls into each packet:
 *
 *   [u_{k|k}, u_{k+1|k}, ..., u_{k+H-1|k}]
 *
 * If current packet is lost, the actuator uses the next prediction
 * from the most recently received buffer (shifts by one step).
 *
 * Stability: guaranteed if consecutive losses ≤ H and the system
 * is open-loop stable over H steps after buffer exhaustion.
 *
 * Buffer generation: u_{t|k} = -L·(A_cl)^t·x̂_k
 * where A_cl = A - B·L is the closed-loop matrix.
 * ============================================================================ */

PacketizedPredictiveControl* pl_ppc_create(int horizon,
                                             const double* A, const double* B,
                                             const double* L,
                                             int n, int m,
                                             const double* x0) {
    PacketizedPredictiveControl* ppc = (PacketizedPredictiveControl*)calloc(1,
        sizeof(PacketizedPredictiveControl));
    if (!ppc) return NULL;

    ppc->horizon = (horizon > 0) ? horizon : 1;
    ppc->n = n; ppc->m = m;
    ppc->A = A; ppc->B = B; ppc->L = L;

    ppc->buffer = (double**)malloc(ppc->horizon * sizeof(double*));
    double* buf_data = _m(ppc->horizon, m);
    for (int i = 0; i < ppc->horizon; i++)
        ppc->buffer[i] = buf_data + i * m;

    ppc->current_input = (double*)calloc(m, sizeof(double));
    ppc->emergency_input = (double*)calloc(m, sizeof(double));
    ppc->x_current = (double*)malloc(n * sizeof(double));
    if (x0) memcpy(ppc->x_current, x0, n * sizeof(double));

    ppc->buffer_index = 0;
    ppc->consecutive_losses = 0;
    ppc->max_packet_loss_tolerated = ppc->horizon;
    ppc->buffer_exhausted = false;

    /* Precompute closed-loop matrix for prediction */
    ppc->buffer_exhausted = false;

    /* Generate initial buffer */
    pl_ppc_generate_buffer(ppc);
    ppc->buffer_index = 0;

    return ppc;
}

void pl_ppc_free(PacketizedPredictiveControl* ppc) {
    if (!ppc) return;
    free(ppc->buffer[0]);
    free(ppc->buffer);
    free(ppc->current_input);
    free(ppc->emergency_input);
    free(ppc->x_current);
    free(ppc);
}

void pl_ppc_generate_buffer(PacketizedPredictiveControl* ppc) {
    if (!ppc) return;

    /* Compute A_cl = A - B·L */
    double* Acl = _m(ppc->n, ppc->n);
    double* BL  = _m(ppc->n, ppc->n);
    for (int i = 0; i < ppc->n; i++)
        for (int j = 0; j < ppc->n; j++) {
            double s = 0.0;
            for (int l = 0; l < ppc->m; l++)
                s += ppc->B[i * ppc->m + l] * ppc->L[l * ppc->n + j];
            BL[i * ppc->n + j] = s;
            Acl[i * ppc->n + j] = ppc->A[i * ppc->n + j] - s;
        }

    /* Compute powers: A_cl^0 = I, A_cl^1, ..., A_cl^{H-1} */
    double* Apow = _m(ppc->n, ppc->n);
    for (int i = 0; i < ppc->n * ppc->n; i++) Apow[i] = (i % (ppc->n + 1) == 0) ? 1.0 : 0.0;

    for (int t = 0; t < ppc->horizon; t++) {
        /* x_pred = A_cl^t · x_current */
        double* x_pred = (double*)malloc(ppc->n * sizeof(double));
        for (int i = 0; i < ppc->n; i++) {
            x_pred[i] = 0.0;
            for (int j = 0; j < ppc->n; j++)
                x_pred[i] += Apow[i * ppc->n + j] * ppc->x_current[j];
        }

        /* u = -L·x_pred */
        for (int i = 0; i < ppc->m; i++) {
            double sum = 0.0;
            for (int j = 0; j < ppc->n; j++)
                sum += ppc->L[i * ppc->n + j] * x_pred[j];
            ppc->buffer[t][i] = -sum;
        }
        free(x_pred);

        /* Update Apow = Apow · A_cl for next iteration */
        if (t < ppc->horizon - 1) {
            double* new_Apow = _m(ppc->n, ppc->n);
            for (int i = 0; i < ppc->n; i++)
                for (int j = 0; j < ppc->n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < ppc->n; l++)
                        s += Apow[i * ppc->n + l] * Acl[l * ppc->n + j];
                    new_Apow[i * ppc->n + j] = s;
                }
            memcpy(Apow, new_Apow, ppc->n * ppc->n * sizeof(double));
            free(new_Apow);
        }
    }

    free(Acl); free(BL); free(Apow);
}

const double* pl_ppc_consume(PacketizedPredictiveControl* ppc,
                              bool new_packet_arrived) {
    if (!ppc) return NULL;

    ppc->total_packets_sent++;

    if (new_packet_arrived) {
        /* Regenerate buffer with updated state, reset index */
        pl_ppc_generate_buffer(ppc);
        ppc->buffer_index = 0;
        ppc->consecutive_losses = 0;
        ppc->buffer_exhausted = false;
        memcpy(ppc->current_input, ppc->buffer[0], ppc->m * sizeof(double));
    } else {
        /* Lost packet: advance buffer index */
        ppc->consecutive_losses++;
        ppc->total_packets_used_from_buffer++;

        if (ppc->consecutive_losses > ppc->max_consecutive_losses_seen)
            ppc->max_consecutive_losses_seen = ppc->consecutive_losses;

        if (ppc->buffer_index < ppc->horizon) {
            memcpy(ppc->current_input, ppc->buffer[ppc->buffer_index],
                   ppc->m * sizeof(double));
            ppc->buffer_index++;
        } else {
            /* Buffer exhausted: use emergency input (zero) */
            ppc->buffer_exhausted = true;
            memcpy(ppc->current_input, ppc->emergency_input,
                   ppc->m * sizeof(double));
        }
    }

    /* Update avg buffer utilization */
    double alpha = 0.01;
    ppc->avg_buffer_utilization = (1.0 - alpha) * ppc->avg_buffer_utilization
                                + alpha * ((double)ppc->buffer_index / ppc->horizon);

    return ppc->current_input;
}

void pl_ppc_update_state(PacketizedPredictiveControl* ppc, const double* x_new) {
    if (ppc && x_new) memcpy(ppc->x_current, x_new, ppc->n * sizeof(double));
}

void pl_ppc_reset(PacketizedPredictiveControl* ppc, const double* x0) {
    if (!ppc) return;
    if (x0) memcpy(ppc->x_current, x0, ppc->n * sizeof(double));
    ppc->buffer_index = 0;
    ppc->consecutive_losses = 0;
    ppc->buffer_exhausted = false;
    pl_ppc_generate_buffer(ppc);
}

void pl_ppc_print(const PacketizedPredictiveControl* ppc) {
    if (!ppc) return;
    printf("=== Packetized Predictive Control ===\n");
    printf("Horizon: %d  Buffer idx: %d  Exhausted: %s\n",
           ppc->horizon, ppc->buffer_index,
           ppc->buffer_exhausted ? "YES" : "NO");
    printf("Consecutive losses: %d (max: %d)\n",
           ppc->consecutive_losses, ppc->max_consecutive_losses_seen);
    printf("Avg buffer util: %.3f\n", ppc->avg_buffer_utilization);
}

/* ============================================================================
 * Model-Based Sensor Predictor
 *
 * Predicts missing sensor measurements using the system model:
 *   ŷ_k = C·A·x̂_{k-1}  (open-loop prediction)
 *
 * When a measurement arrives, corrects the predictor using
 * innovation: x̂ = x̂_pred + K_corr·(y - C·x̂_pred)
 *
 * Prediction quality degrades over consecutive losses:
 *   ||error_N|| ≈ ρ(A)^N·||error_0||
 * ============================================================================ */

ModelSensorPredictor* pl_msp_create(const double* A, const double* C,
                                      int n, int p,
                                      const double* x0, const double* P0) {
    ModelSensorPredictor* msp = (ModelSensorPredictor*)calloc(1,
        sizeof(ModelSensorPredictor));
    if (!msp) return NULL;

    msp->A = A; msp->C = C; msp->n = n; msp->p = p;

    msp->x_pred = (double*)malloc(n * sizeof(double));
    msp->y_pred = (double*)malloc(p * sizeof(double));
    msp->P_pred = _m(n, n);

    if (x0) memcpy(msp->x_pred, x0, n * sizeof(double));
    if (P0) memcpy(msp->P_pred, P0, n * n * sizeof(double));
    else for (int i = 0; i < n; i++) msp->P_pred[i * n + i] = 1.0;

    msp->consecutive_losses = 0;
    msp->correction_gain = 0.5;
    msp->drift_threshold = 1.0;
    msp->model_valid = true;

    return msp;
}

void pl_msp_free(ModelSensorPredictor* msp) {
    if (!msp) return;
    free(msp->x_pred); free(msp->y_pred); free(msp->P_pred);
    free(msp);
}

const double* pl_msp_predict_missing(ModelSensorPredictor* msp) {
    if (!msp) return NULL;

    /* x_pred = A·x_pred (open-loop propagation) */
    double* x_new = (double*)malloc(msp->n * sizeof(double));
    for (int i = 0; i < msp->n; i++) {
        x_new[i] = 0.0;
        for (int j = 0; j < msp->n; j++)
            x_new[i] += msp->A[i * msp->n + j] * msp->x_pred[j];
    }
    memcpy(msp->x_pred, x_new, msp->n * sizeof(double));
    free(x_new);

    /* y_pred = C·x_pred */
    for (int i = 0; i < msp->p; i++) {
        msp->y_pred[i] = 0.0;
        for (int j = 0; j < msp->n; j++)
            msp->y_pred[i] += msp->C[i * msp->n + j] * msp->x_pred[j];
    }

    /* Covariance inflation */
    for (int i = 0; i < msp->n; i++)
        msp->P_pred[i * msp->n + i] *= 1.2;

    msp->consecutive_losses++;

    return msp->y_pred;
}

void pl_msp_correct(ModelSensorPredictor* msp, const double* y_true) {
    if (!msp || !y_true) return;

    /* Compute prediction error */
    double* y_pred = (double*)malloc(msp->p * sizeof(double));
    for (int i = 0; i < msp->p; i++) {
        y_pred[i] = 0.0;
        for (int j = 0; j < msp->n; j++)
            y_pred[i] += msp->C[i * msp->n + j] * msp->x_pred[j];
    }

    double error = 0.0;
    for (int i = 0; i < msp->p; i++) {
        double e = y_true[i] - y_pred[i];
        error += e * e;
    }
    msp->prediction_error = sqrt(error);
    msp->accumulated_drift += msp->prediction_error;
    msp->model_valid = (msp->prediction_error < msp->drift_threshold);

    /* Correct state: x̂ += K_corr·C'(y - C·x̂) */
    for (int i = 0; i < msp->n; i++) {
        msp->x_pred[i] += 0.5 * msp->prediction_error;
    }

    msp->consecutive_losses = 0;
    free(y_pred);
}

const double* pl_msp_step(ModelSensorPredictor* msp,
                           const double* y, bool arrived) {
    if (arrived) pl_msp_correct(msp, y);
    return pl_msp_predict_missing(msp);
}

double pl_msp_get_prediction_error(const ModelSensorPredictor* msp) {
    return msp ? msp->prediction_error : -1.0;
}

bool pl_msp_is_model_valid(const ModelSensorPredictor* msp) {
    return msp ? msp->model_valid : false;
}

void pl_msp_print(const ModelSensorPredictor* msp) {
    if (!msp) return;
    printf("=== Model Sensor Predictor ===\n");
    printf("Consecutive losses: %d\n", msp->consecutive_losses);
    printf("Prediction error: %.6f  Drift: %.6f  Valid: %s\n",
           msp->prediction_error, msp->accumulated_drift,
           msp->model_valid ? "YES" : "NO");
}

/* ============================================================================
 * Hold Strategy Selector
 *
 * Adaptive strategy selection based on loss pattern and system dynamics.
 *
 * Selection logic:
 *   - Open-loop stable + isolated loss → ZERO_ORDER (simple, effective)
 *   - Open-loop stable + burst loss → PREDICTIVE (better compensation)
 *   - Open-loop unstable + any loss → PREDICTIVE or LQG_OPTIMAL
 *   - Very high loss (>50%) → ZERO_INPUT (safety first)
 * ============================================================================ */

HoldStrategySelector* pl_hss_create(bool open_loop_stable, double spectral_radius) {
    HoldStrategySelector* hss = (HoldStrategySelector*)calloc(1,
        sizeof(HoldStrategySelector));
    if (!hss) return NULL;

    hss->open_loop_stable = open_loop_stable;
    hss->spectral_radius = spectral_radius;
    hss->current_strategy = HOLD_ZERO_ORDER;
    hss->history_idx = 0;

    return hss;
}

void pl_hss_free(HoldStrategySelector* hss) { free(hss); }

HoldStrategy pl_hss_select(HoldStrategySelector* hss,
                            double loss_rate, double burst_index,
                            int consecutive_losses,
                            double* estimated_cost_out) {
    if (!hss) return HOLD_ZERO_INPUT;

    /* Safety override: if loss rate is very high, use zero input */
    if (loss_rate > 0.8) {
        hss->current_strategy = HOLD_ZERO_INPUT;
        if (estimated_cost_out) *estimated_cost_out = 1e6;
        return hss->current_strategy;
    }

    /* Estimate costs for each strategy (simplified model) */
    double cost_zoh, cost_pred, cost_lqg, cost_safety;

    cost_zoh = loss_rate * 10.0;
    cost_pred = loss_rate * 5.0 + (hss->open_loop_stable ? 0.0 : 2.0);
    cost_lqg  = loss_rate * 3.0 + 1.0;  /* Best but computational overhead */
    cost_safety = 100.0;                 /* Zero input is always worst for perf */

    /* Burst penalty */
    if (burst_index > 2.0) {
        cost_zoh *= 1.5;  /* ZOH degrades under bursts */
        cost_pred *= 1.1;  /* Predictive handles bursts better */
    }

    /* Consecutive loss penalty */
    if (consecutive_losses > 3) {
        cost_zoh *= (1.0 + 0.2 * consecutive_losses);
        cost_pred *= (1.0 + 0.05 * consecutive_losses);
    }

    /* Open-loop stability adjustment */
    if (!hss->open_loop_stable) {
        cost_zoh *= 2.0;   /* ZOH unreliable for unstable systems */
        cost_pred *= 1.2;
    }

    /* Select minimum cost strategy */
    HoldStrategy best = HOLD_ZERO_ORDER;
    double min_cost = cost_zoh;

    if (cost_pred < min_cost) { min_cost = cost_pred; best = HOLD_PREDICTIVE; }
    if (cost_lqg < min_cost)  { min_cost = cost_lqg;  best = HOLD_LQG_OPTIMAL; }
    if (cost_safety < min_cost) { min_cost = cost_safety; best = HOLD_ZERO_INPUT; }

    hss->current_strategy = best;
    hss->recent_loss_rate = loss_rate;
    hss->burst_index = burst_index;
    hss->consecutive_losses = consecutive_losses;
    hss->zero_order_cost = cost_zoh;
    hss->predictive_cost = cost_pred;
    hss->lqg_optimal_cost = cost_lqg;
    hss->safety_cost = cost_safety;

    if (estimated_cost_out) *estimated_cost_out = min_cost;

    return best;
}

void pl_hss_record_performance(HoldStrategySelector* hss,
                                HoldStrategy strategy, double cost) {
    if (!hss) return;
    int idx = hss->history_idx % 100;
    switch (strategy) {
        case HOLD_ZERO_ORDER:  hss->cost_zoh_history[idx] = cost; break;
        case HOLD_PREDICTIVE:  hss->cost_pred_history[idx] = cost; break;
        case HOLD_LQG_OPTIMAL: hss->cost_lqg_history[idx] = cost; break;
        default: break;
    }
    hss->history_idx++;
}

void pl_hss_print(const HoldStrategySelector* hss) {
    if (!hss) return;
    printf("=== Hold Strategy Selector ===\n");
    printf("Open-loop stable: %s  ρ(A)=%.4f\n",
           hss->open_loop_stable ? "YES" : "NO", hss->spectral_radius);
    printf("Current strategy: %d\n", hss->current_strategy);
    printf("Costs: ZOH=%.3f  Pred=%.3f  LQG=%.3f  Safety=%.3f\n",
           hss->zero_order_cost, hss->predictive_cost,
           hss->lqg_optimal_cost, hss->safety_cost);
}

/* ============================================================================
 * Reference Governor for NCS with Packet Loss
 *
 * Modifies reference to ensure constraint satisfaction under worst-case
 * packet loss scenario. If the predicted trajectory under H consecutive
 * losses violates constraints, the reference is scaled back.
 *
 * Reference: Bemporad (1998), "Predictive Control of Teleoperated
 * Constrained Systems with Bounded Time-Varying Delays"
 * ============================================================================ */

ReferenceGovernor* pl_rg_create(const double* A, const double* B, const double* C,
                                  int n, int m, int p,
                                  const double* reference, int ref_size,
                                  const double* y_min, const double* y_max,
                                  int max_loss) {
    ReferenceGovernor* rg = (ReferenceGovernor*)calloc(1, sizeof(ReferenceGovernor));
    if (!rg) return NULL;

    rg->A = A; rg->B = B; rg->C = C;
    rg->n = n; rg->m = m; rg->p = p;
    rg->ref_size = ref_size;
    rg->n_constraints = p;
    rg->max_consecutive_losses = max_loss;
    rg->safety_factor = 0.8;

    rg->reference = (double*)malloc(ref_size * sizeof(double));
    rg->modified_reference = (double*)malloc(ref_size * sizeof(double));
    rg->y_min = (double*)malloc(p * sizeof(double));
    rg->y_max = (double*)malloc(p * sizeof(double));
    rg->x_pred = (double*)malloc(n * sizeof(double));

    if (reference) memcpy(rg->reference, reference, ref_size * sizeof(double));
    if (y_min) memcpy(rg->y_min, y_min, p * sizeof(double));
    else for (int i = 0; i < p; i++) rg->y_min[i] = -1e6;
    if (y_max) memcpy(rg->y_max, y_max, p * sizeof(double));
    else for (int i = 0; i < p; i++) rg->y_max[i] = 1e6;

    return rg;
}

void pl_rg_free(ReferenceGovernor* rg) {
    if (!rg) return;
    free(rg->reference); free(rg->modified_reference);
    free(rg->y_min); free(rg->y_max); free(rg->x_pred);
    free(rg);
}

const double* pl_rg_adjust(ReferenceGovernor* rg, const double* x_current) {
    if (!rg) return NULL;

    /* Predict worst-case trajectory: no control for H steps */
    double* x_pred = (double*)malloc(rg->n * sizeof(double));
    memcpy(x_pred, x_current, rg->n * sizeof(double));

    bool constraint_violated = false;

    for (int step = 0; step < rg->max_consecutive_losses; step++) {
        /* Open-loop: x_next = A·x (no control input assumed) */
        double* x_next = (double*)malloc(rg->n * sizeof(double));
        for (int i = 0; i < rg->n; i++) {
            x_next[i] = 0.0;
            for (int j = 0; j < rg->n; j++)
                x_next[i] += rg->A[i * rg->n + j] * x_pred[j];
        }
        memcpy(x_pred, x_next, rg->n * sizeof(double));
        free(x_next);

        /* Check output constraints: y = C·x */
        for (int i = 0; i < rg->p; i++) {
            double y = 0.0;
            for (int j = 0; j < rg->n; j++)
                y += rg->C[i * rg->n + j] * x_pred[j];
            if (y < rg->y_min[i] || y > rg->y_max[i]) {
                constraint_violated = true;
                break;
            }
        }
        if (constraint_violated) break;
    }

    /* Adjust reference: scale by safety factor if violation predicted */
    for (int i = 0; i < rg->ref_size; i++) {
        if (constraint_violated)
            rg->modified_reference[i] = rg->safety_factor * rg->reference[i];
        else
            rg->modified_reference[i] = rg->reference[i];
    }

    free(x_pred);
    return rg->modified_reference;
}

void pl_rg_print(const ReferenceGovernor* rg) {
    if (!rg) return;
    printf("=== Reference Governor ===\n");
    printf("Max consecutive losses: %d  Safety factor: %.2f\n",
           rg->max_consecutive_losses, rg->safety_factor);
    printf("Constraints: y ∈ [%.2f, %.2f]\n", rg->y_min[0], rg->y_max[0]);
}

/* ============================================================================
 * Loss-Constrained Optimal Controller
 *
 * Minimizes expected LQG cost under Bernoulli packet loss.
 * For TCP-like: optimal policy is certainty-equivalent LQR.
 * For UDP-like (Gupta et al., 2007): gains depend on loss probability.
 *
 * The cost-to-go recursion:
 *   P_N = Q (terminal cost)
 *   P_t = Q + A'·[(1-p)·P_{t+1} + p·A'·P_{t+1}·A]·A  (UDP-like, simplified)
 *   K_t = (R + B'·P_{t+1}·B)^{-1}·B'·P_{t+1}·A
 * ============================================================================ */

LossConstrainedController* pl_lcc_create(const double* A, const double* B,
                                           const double* Q, const double* R,
                                           int n, int m, int horizon,
                                           double loss_prob,
                                           TransportProtocol protocol) {
    LossConstrainedController* lcc = (LossConstrainedController*)calloc(1,
        sizeof(LossConstrainedController));
    if (!lcc) return NULL;

    lcc->A = A; lcc->B = B; lcc->Q = Q; lcc->R = R;
    lcc->n = n; lcc->m = m;
    lcc->horizon = horizon;
    lcc->loss_probability = (loss_prob < 0.0) ? 0.0 :
                            (loss_prob > 1.0) ? 1.0 : loss_prob;
    lcc->protocol = protocol;

    lcc->gain_sequence = (double**)malloc(horizon * sizeof(double*));
    double* gs_data = _m(horizon * m, n);
    for (int i = 0; i < horizon; i++)
        lcc->gain_sequence[i] = gs_data + i * m * n;

    lcc->cost_to_go = (double**)malloc(horizon * sizeof(double*));
    double* ctg_data = _m(horizon * n, n);
    for (int i = 0; i < horizon; i++)
        lcc->cost_to_go[i] = ctg_data + i * n * n;

    lcc->u_optimal = (double*)calloc(m, sizeof(double));

    /* Compute cost-to-go and gains backward in time */
    double* P = _m(n, n);
    for (int i = 0; i < n * n; i++) P[i] = Q[i];  /* Terminal cost = Q */

    double* BT = _m(m, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) BT[j * n + i] = B[i * m + j];

    for (int t = horizon - 1; t >= 0; t--) {
        memcpy(lcc->cost_to_go[t], P, n * n * sizeof(double));

        /* K_t = (R + B'·P·B)^{-1}·B'·P·A */
        double* BPB = _m(m, m);
        double* PB  = _m(n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += P[i * n + l] * B[l * m + j];
                PB[i * m + j] = s;
            }
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += BT[i * n + l] * PB[l * m + j];
                BPB[i * m + j] = s + R[i * m + j];
            }

        double* PA = _m(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += P[i * n + l] * A[l * n + j];
                PA[i * n + j] = s;
            }

        double* BPA = _m(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < n; l++) s += BT[i * n + l] * PA[l * n + j];
                BPA[i * n + j] = s;
            }

        double* BPB_inv = (double*)malloc(m * m * sizeof(double));
        memcpy(BPB_inv, BPB, m * m * sizeof(double));
        /* Simple inversion (2×2 or general) for small m */
        if (m == 1) {
            if (fabs(BPB_inv[0]) > 1e-15) BPB_inv[0] = 1.0 / BPB_inv[0];
        }
        /* For general m, use inline inverse when small or skip */

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int l = 0; l < m; l++) s += BPB_inv[i * m + l] * BPA[l * n + j];
                lcc->gain_sequence[t][i * n + j] = s;
            }

        /* Update P for previous time step (UDP-like modified) */
        if (protocol == PROTO_UDP_LIKE) {
            /* P_new = Q + A'·[(1-p)·P + p·A'·P·A]·A */
            double* APA = _m(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < n; l++) s += A[l * n + i] * P[l * n + j];
                    PA[i * n + j] = s;
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < n; l++) s += PA[i*n+l] * A[l*n+j];
                    APA[i * n + j] = s;
                }
            double* P_inner = _m(n, n);
            for (int i = 0; i < n * n; i++)
                P_inner[i] = (1.0 - lcc->loss_probability) * P[i]
                           + lcc->loss_probability * APA[i];
            double* AT = _m(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) AT[i * n + j] = A[j * n + i];
            double* temp = _m(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < n; l++) s += AT[i*n+l] * P_inner[l*n+j];
                    temp[i * n + j] = s;
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < n; l++) s += temp[i*n+l] * A[l*n+j];
                    P[i * n + j] = s + Q[i * n + j];
                }
            free(APA); free(P_inner); free(AT); free(temp);
        } else {
            /* TCP-like: standard Riccati */
            /* P = Q + A'·P·A - A'·P·B·(R+B'·P·B)^{-1}·B'·P·A */
            /* Simplified: just use Q + A'PA for now */
            double* AT2 = _m(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) AT2[i * n + j] = A[j * n + i];
            double* AP2 = _m(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < n; l++) s += AT2[i*n+l] * P[l*n+j];
                    AP2[i * n + j] = s;
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    double s = 0.0;
                    for (int l = 0; l < n; l++) s += AP2[i*n+l] * A[l*n+j];
                    P[i * n + j] = s + Q[i * n + j];
                }
            free(AT2); free(AP2);
        }

        free(PB); free(BPB); free(PA); free(BPA); free(BPB_inv);
    }

    free(P); free(BT);

    return lcc;
}

void pl_lcc_free(LossConstrainedController* lcc) {
    if (!lcc) return;
    free(lcc->gain_sequence[0]);
    free(lcc->gain_sequence);
    free(lcc->cost_to_go[0]);
    free(lcc->cost_to_go);
    free(lcc->u_optimal);
    free(lcc);
}

const double* pl_lcc_compute(LossConstrainedController* lcc,
                              const double* x_estimate) {
    if (!lcc || !x_estimate) return NULL;

    /* u_optimal = -K_0·x_estimate (first gain in sequence) */
    for (int i = 0; i < lcc->m; i++) {
        lcc->u_optimal[i] = 0.0;
        for (int j = 0; j < lcc->n; j++)
            lcc->u_optimal[i] -= lcc->gain_sequence[0][i * lcc->n + j]
                                * x_estimate[j];
    }

    return lcc->u_optimal;
}

void pl_lcc_print(const LossConstrainedController* lcc) {
    if (!lcc) return;
    printf("=== Loss-Constrained Controller ===\n");
    printf("Horizon: %d  Loss prob: %.4f  Protocol: %s\n",
           lcc->horizon, lcc->loss_probability,
           lcc->protocol == PROTO_TCP_LIKE ? "TCP-like" : "UDP-like");
}
