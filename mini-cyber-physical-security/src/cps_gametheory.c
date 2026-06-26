#include "cps_gametheory.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CPS_EPS 1e-12
#define PI 3.14159265358979323846

/* ============================================================================
 * Game-Theoretic CPS Security (L8: Advanced Topics)
 *
 * Models the attacker-defender interaction using game theory.
 *
 * In a zero-sum game, the defender minimizes a cost function J(u,a)
 * while the attacker maximizes it. The saddle-point equilibrium (u*,a*)
 * satisfies: J(u*,a) <= J(u*,a*) <= J(u,a*) for all u, a.
 *
 * This provides a principled framework for:
 *   - Optimal defense resource allocation
 *   - Worst-case attack analysis
 *   - Security investment planning
 *
 * Reference: Zhu & Basar (2015), "Game-Theoretic Methods for Robustness,
 *            Security, and Resilience of Cyberphysical Control Systems"
 * ============================================================================ */

/* ============================================================================
 * 2x2 Zero-Sum Game Solver (L5: Algorithms)
 *
 * Solves the classic 2x2 zero-sum matrix game analytically.
 * The game value V and optimal strategies are computed by solving
 * a 2x2 linear system derived from the equilibrium conditions.
 *
 * Payoff matrix (attacker's payoff, defender minimizes):
 *     Attacker
 *        L     R
 * D  U [a11, a12]
 *    D [a21, a22]
 *
 * Defender's strategy: (p, 1-p) for Up, Down
 * Attacker's strategy: (q, 1-q) for Left, Right
 * ============================================================================ */

void cps_solve_2x2_zerosum(double a11, double a12, double a21, double a22,
                            double* saddle_value,
                            double* defender_p1,
                            double* attacker_q1) {
    /* Compute the saddle point analytically.
     * Defender's expected payoff: V_def = p*q*a11 + p*(1-q)*a12
     *                                    + (1-p)*q*a21 + (1-p)*(1-q)*a22 */
    double denom = a11 - a12 - a21 + a22;

    if (fabs(denom) > CPS_EPS) {
        /* Mixed strategy equilibrium exists */
        double p_opt = (a22 - a12) / denom;
        double q_opt = (a22 - a21) / denom;

        /* Clamp to [0,1] (pure strategy if outside) */
        if (p_opt < 0.0) p_opt = 0.0;
        if (p_opt > 1.0) p_opt = 1.0;
        if (q_opt < 0.0) q_opt = 0.0;
        if (q_opt > 1.0) q_opt = 1.0;

        *defender_p1 = p_opt;
        *attacker_q1 = q_opt;
        *saddle_value = p_opt * q_opt * a11 + p_opt * (1.0 - q_opt) * a12
                      + (1.0 - p_opt) * q_opt * a21
                      + (1.0 - p_opt) * (1.0 - q_opt) * a22;
    } else {
        /* No unique mixed strategy — use pure minimax */
        double min_row_max = (a11 < a12) ? a11 : a12;
        double min_row_max2 = (a21 < a22) ? a21 : a22;
        if (min_row_max > min_row_max2) {
            *defender_p1 = 0.0;  /* Choose row 1 */
            *saddle_value = min_row_max2;
            *attacker_q1 = (a21 < a22) ? 0.0 : 1.0;
        } else {
            *defender_p1 = 1.0;
            *saddle_value = min_row_max;
            *attacker_q1 = (a11 < a12) ? 0.0 : 1.0;
        }
    }
}

/* ============================================================================
 * General Zero-Sum Game via Linear Programming (L8: Advanced)
 *
 * The minimax theorem (von Neumann, 1928) states that every finite
 * zero-sum game has a saddle point in mixed strategies.
 *
 * We solve the LP: max v subject to:
 *   sum_i x_i * A_{ij} >= v  for all j
 *   sum_i x_i = 1
 *   x_i >= 0
 *
 * Using a simplified simplex-like iterative method.
 * Alternatively, the fictitious play algorithm (Brown, 1951) is used
 * here for its simplicity and guaranteed convergence.
 * ============================================================================ */

int cps_solve_zerosum_lp(CPSZeroSumGame* game, int max_iter) {
    if (!game || game->n_defender_actions <= 0
        || game->n_attacker_actions <= 0) return -1;

    int n_d = game->n_defender_actions;
    int n_a = game->n_attacker_actions;

    /* Fictitious play: iteratively update empirical strategy */
    double* def_cumulative = (double*)calloc((size_t)n_d, sizeof(double));
    double* att_cumulative = (double*)calloc((size_t)n_a, sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        /* Attacker best-responds to defender's empirical strategy */
        double* def_empirical = (double*)malloc(
            (size_t)n_d * sizeof(double));
        for (int i = 0; i < n_d; i++)
            def_empirical[i] = def_cumulative[i]
                / (iter + 1.0 + CPS_EPS);

        int att_best = 0;
        double att_best_payoff = -1e300;
        for (int j = 0; j < n_a; j++) {
            double expected = 0.0;
            for (int i = 0; i < n_d; i++)
                expected += def_empirical[i]
                    * game->payoff_matrix[i * n_a + j];
            if (expected > att_best_payoff) {
                att_best_payoff = expected;
                att_best = j;
            }
        }
        att_cumulative[att_best] += 1.0;

        /* Defender best-responds to attacker's empirical strategy */
        double* att_empirical = (double*)malloc(
            (size_t)n_a * sizeof(double));
        for (int j = 0; j < n_a; j++)
            att_empirical[j] = att_cumulative[j]
                / (iter + 1.0 + CPS_EPS);

        int def_best = 0;
        double def_best_payoff = 1e300;
        for (int i = 0; i < n_d; i++) {
            double expected = 0.0;
            for (int j = 0; j < n_a; j++)
                expected += att_empirical[j]
                    * game->payoff_matrix[i * n_a + j];
            if (expected < def_best_payoff) {
                def_best_payoff = expected;
                def_best = i;
            }
        }
        def_cumulative[def_best] += 1.0;

        free(def_empirical);
        free(att_empirical);
    }

    /* Final strategies */
    double total_iter = (double)max_iter;
    if (game->defender_mixed) free(game->defender_mixed);
    if (game->attacker_mixed) free(game->attacker_mixed);
    game->defender_mixed = (double*)malloc(
        (size_t)n_d * sizeof(double));
    game->attacker_mixed = (double*)malloc(
        (size_t)n_a * sizeof(double));

    for (int i = 0; i < n_d; i++)
        game->defender_mixed[i] = def_cumulative[i] / total_iter;
    for (int j = 0; j < n_a; j++)
        game->attacker_mixed[j] = att_cumulative[j] / total_iter;

    /* Compute saddle value */
    game->saddle_value = 0.0;
    for (int i = 0; i < n_d; i++)
        for (int j = 0; j < n_a; j++)
            game->saddle_value += game->defender_mixed[i]
                * game->attacker_mixed[j]
                * game->payoff_matrix[i * n_a + j];
    game->solved = 1;

    free(def_cumulative);
    free(att_cumulative);
    return 0;
}

int cps_zerosum_nash_equilibrium(CPSZeroSumGame* game) {
    return cps_solve_zerosum_lp(game, 1000);
}

/* ============================================================================
 * Stackelberg Security Game (L8: Advanced)
 *
 * In a Stackelberg game, the defender (leader) commits to a strategy
 * first, and the attacker (follower) best-responds. The Strong
 * Stackelberg Equilibrium (SSE) assumes the attacker breaks ties
 * in favor of the defender.
 *
 * This models security scenarios where the defender deploys a
 * fixed security posture (e.g., sensor placement, patrolling schedule)
 * that the attacker can observe before acting.
 * ============================================================================ */

int cps_stackelberg_sse(CPSZeroSumGame* game,
                         double* defender_strategy,
                         double* attacker_best_response,
                         double* sse_value) {
    if (!game || !defender_strategy || !sse_value) return -1;
    int n_d = game->n_defender_actions;
    int n_a = game->n_attacker_actions;

    /* Try each attacker pure strategy as best response.
     * For each attacker action j, solve:
     *   maximize sum_i x_i * A_{ij}
     *   subject to: sum_i x_i * A_{ij} >= sum_i x_i * A_{ik} for all k
     *               sum_i x_i = 1, x_i >= 0 */
    double best_sse = -1e300;
    int best_att_action = 0;
    double* best_def_strat = (double*)malloc(
        (size_t)n_d * sizeof(double));

    for (int j = 0; j < n_a; j++) {
        /* Simplified: defender puts all probability on best pure
         * strategy for this attacker response */
        int best_def = 0;
        double best_pay = game->payoff_matrix[0 * n_a + j];
        for (int i = 1; i < n_d; i++) {
            double pay = game->payoff_matrix[i * n_a + j];
            if (pay < best_pay) { best_pay = pay; best_def = i; }
        }

        /* Check if attacker j is indeed a best response */
        int is_best_response = 1;
        for (int k = 0; k < n_a && is_best_response; k++) {
            if (k == j) continue;
            if (game->payoff_matrix[best_def * n_a + k]
                > game->payoff_matrix[best_def * n_a + j] + CPS_EPS)
                is_best_response = 0;
        }

        if (is_best_response && best_pay < best_sse + CPS_EPS) {
            /* Defender minimizes attacker payoff */
            if (best_pay > best_sse) {
                best_sse = best_pay;
                best_att_action = j;
                for (int i = 0; i < n_d; i++)
                    best_def_strat[i] = (i == best_def) ? 1.0 : 0.0;
            }
        }
    }

    if (best_sse > -1e299) {
        for (int i = 0; i < n_d; i++)
            defender_strategy[i] = best_def_strat[i];
        if (attacker_best_response) {
            for (int j = 0; j < n_a; j++)
                attacker_best_response[j] = (j == best_att_action)
                    ? 1.0 : 0.0;
        }
        *sse_value = best_sse;
        free(best_def_strat);
        return 0;
    }

    free(best_def_strat);
    return -1;
}

double cps_defender_utility_breach(CPSZeroSumGame* game,
                                    const double* defender_strategy) {
    if (!game || !defender_strategy) return 0.0;
    int n_d = game->n_defender_actions;
    int n_a = game->n_attacker_actions;

    /* Attacker best-responds */
    double max_payoff = -1e300;
    for (int j = 0; j < n_a; j++) {
        double expected = 0.0;
        for (int i = 0; i < n_d; i++)
            expected += defender_strategy[i]
                * game->payoff_matrix[i * n_a + j];
        if (expected > max_payoff) max_payoff = expected;
    }
    return max_payoff;
}

/* ============================================================================
 * Dynamic (Multi-stage) Security Game (L8: Advanced)
 * ============================================================================ */

void cps_dynamic_game_init(CPSDynamicGame* game, int n_states, int n_stages) {
    if (!game) return;
    game->n_states = n_states;
    game->n_stages = n_stages;
    game->discount_factor = 0.95;
    game->state_values = (double*)calloc(
        (size_t)(n_states * (n_stages + 1)), sizeof(double));
    game->defender_policy = (int*)malloc(
        (size_t)(n_states * n_stages) * sizeof(int));
    game->attacker_policy = (int*)malloc(
        (size_t)(n_states * n_stages) * sizeof(int));
    for (int i = 0; i < n_states * n_stages; i++) {
        game->defender_policy[i] = 0;
        game->attacker_policy[i] = 0;
    }
}

void cps_dynamic_game_solve(CPSDynamicGame* game) {
    if (!game || game->n_states <= 0 || game->n_stages <= 0) return;

    /* Backward induction: V_T(s) = 0, then
     * V_k(s) = min_{a_d} max_{a_a} [
     *   r(s,a_d,a_a) + gamma * sum_{s'} P(s'|s,a_d,a_a) * V_{k+1}(s') ]
     *
     * With simplified transition model: s' = s (state doesn't change
     * for single-stage analysis) */
    int S = game->n_states;
    int T = game->n_stages;
    double gamma = game->discount_factor;

    for (int s = 0; s < S; s++)
        game->state_values[s * (T + 1) + T] = 0.0;

    for (int k = T - 1; k >= 0; k--) {
        for (int s = 0; s < S; s++) {
            /* For each state, solve a zero-sum game at this stage */
            int n_d = 2; /* Assume binary choices for simplicity */
            int n_a = 2;
            double payoff[4]; /* 2x2 matrix */

            /* Payoff depends on state: higher state = more vulnerable */
            double state_factor = 1.0 + 0.5 * s;
            payoff[0] = 1.0 * state_factor;  /* def=0, att=0 */
            payoff[1] = 5.0 * state_factor;  /* def=0, att=1 */
            payoff[2] = 2.0 * state_factor;  /* def=1, att=0 */
            payoff[3] = 3.0 * state_factor;  /* def=1, att=1 */

            double p_def, q_att, val;
            cps_solve_2x2_zerosum(payoff[0], payoff[1],
                                   payoff[2], payoff[3],
                                   &val, &p_def, &q_att);

            game->state_values[s * (T + 1) + k] = val
                + gamma * game->state_values[s * (T + 1) + k + 1];
            game->defender_policy[s * T + k] = (p_def > 0.5) ? 1 : 0;
            game->attacker_policy[s * T + k] = (q_att > 0.5) ? 1 : 0;
        }
    }
}

void cps_dynamic_game_policy(CPSDynamicGame* game, int state, int stage,
                              int* def_action, int* att_action) {
    if (!game || state < 0 || state >= game->n_states
        || stage < 0 || stage >= game->n_stages) {
        if (def_action) *def_action = 0;
        if (att_action) *att_action = 0;
        return;
    }
    int T = game->n_stages;
    if (def_action) *def_action = game->defender_policy[state * T + stage];
    if (att_action) *att_action = game->attacker_policy[state * T + stage];
}

/* ============================================================================
 * Attack-Defense Cost Models (L8: Advanced)
 * ============================================================================ */

double cps_lqg_cost_under_attack(const double* A, const double* B,
                                  const double* Q, const double* R,
                                  const double* Sigma_a,
                                  int n, int m) {
    if (!A || !B || n <= 0 || m <= 0) return 0.0;

    /* Trace approximation of LQG cost with attack covariance Sigma_a
     * J = trace(P * (B*Sigma_a*B')) + trace(Q * Sigma_x)
     * where P solves: P = A'*P*A + Q - A'*P*B*(R+B'*P*B)^{-1}*B'*P*A */
    double cost = 0.0;
    for (int i = 0; i < n; i++) {
        double sigma_ii = Sigma_a ? fabs(Sigma_a[i * n + i]) : 1.0;
        cost += sigma_ii * (Q ? fabs(Q[0]) : 1.0);
    }
    return cost;
}

double cps_optimal_detection_threshold(double cost_fp, double cost_fn,
                                        double cost_invest_coeff,
                                        double attack_prior,
                                        double detection_snr) {
    /* Minimize: C_total = C_fp * P_fp(th) + C_fn * P_fn(th) + C_inv * th
     * where P_fp(th) = P(N(0,1) > th), P_fn(th) = P(N(SNR,1) < th)
     * Approximate via Gaussian CDF */
    double best_threshold = 3.0;
    double best_cost = 1e300;

    for (int i = 1; i <= 100; i++) {
        double th = 0.1 * i;
        double p_fp = 0.5 * (1.0 - erf(th / sqrt(2.0)));
        double p_fn = 0.5 * (1.0 + erf((th - detection_snr) / sqrt(2.0)));
        double total = cost_fp * p_fp * (1.0 - attack_prior)
                     + cost_fn * p_fn * attack_prior
                     + cost_invest_coeff * th;

        if (total < best_cost) {
            best_cost = total;
            best_threshold = th;
        }
    }
    return best_threshold;
}

/* ============================================================================
 * Bayesian Attacker Belief (L8: Advanced)
 * ============================================================================ */

void cps_bayesian_attacker_belief(double* belief, int n_types,
                                   const double* observation_likelihood,
                                   int observed_action) {
    if (!belief || n_types <= 0) return;

    /* Bayesian update: P(type | obs) proportional to P(obs | type) * P(type) */
    double total = 0.0;
    for (int t = 0; t < n_types; t++) {
        double likelihood = observation_likelihood
            ? observation_likelihood[t * n_types + observed_action]
            : 1.0 / n_types;
        belief[t] *= likelihood + CPS_EPS;
        total += belief[t];
    }
    if (total > CPS_EPS) {
        for (int t = 0; t < n_types; t++)
            belief[t] /= total;
    }
}

void cps_defender_robust_strategy(const CPSZeroSumGame** games,
                                   int n_types,
                                   const double* type_belief,
                                   double* robust_strategy) {
    if (!games || !type_belief || !robust_strategy || n_types <= 0)
        return;

    /* Bayesian-Nash equilibrium: maximize expected utility */
    int n_d = games[0]->n_defender_actions;

    /* Take weighted average of strategies */
    for (int i = 0; i < n_d; i++) {
        robust_strategy[i] = 0.0;
        for (int t = 0; t < n_types; t++) {
            if (games[t] && games[t]->defender_mixed)
                robust_strategy[i] += type_belief[t]
                    * games[t]->defender_mixed[i];
        }
    }

    /* Normalize */
    double sum = 0.0;
    for (int i = 0; i < n_d; i++) sum += robust_strategy[i];
    if (sum > CPS_EPS)
        for (int i = 0; i < n_d; i++) robust_strategy[i] /= sum;
}

/* ============================================================================
 * Lifecycle Management
 * ============================================================================ */

void cps_gameplayer_init(CPSGamePlayer* player, int n_actions,
                          int is_attacker) {
    if (!player) return;
    player->n_actions = n_actions;
    player->is_attacker = is_attacker;
    player->payoff = 0.0;
    player->resource_budget = 1.0;
    player->strategy = (double*)calloc((size_t)n_actions, sizeof(double));
    player->action_costs = (double*)calloc((size_t)n_actions,
                                            sizeof(double));
    if (n_actions > 0) player->strategy[0] = 1.0;
}

void cps_gameplayer_free(CPSGamePlayer* player) {
    if (!player) return;
    free(player->strategy);
    free(player->action_costs);
}

void cps_zerosum_game_init(CPSZeroSumGame* game, int n_def, int n_att) {
    if (!game) return;
    game->n_defender_actions = n_def;
    game->n_attacker_actions = n_att;
    game->saddle_value = 0.0;
    game->solved = 0;
    game->payoff_matrix = (double*)calloc(
        (size_t)(n_def * n_att), sizeof(double));
    game->defender_mixed = (double*)calloc((size_t)n_def, sizeof(double));
    game->attacker_mixed = (double*)calloc((size_t)n_att, sizeof(double));
    if (n_def > 0) game->defender_mixed[0] = 1.0;
    if (n_att > 0) game->attacker_mixed[0] = 1.0;
}

void cps_zerosum_game_free(CPSZeroSumGame* game) {
    if (!game) return;
    free(game->payoff_matrix);
    free(game->defender_mixed);
    free(game->attacker_mixed);
}

void cps_dynamic_game_free(CPSDynamicGame* game) {
    if (!game) return;
    free(game->stage_games);
    free(game->state_values);
    free(game->defender_policy);
    free(game->attacker_policy);
    if (game->transition) {
        for (int i = 0; i < game->n_states; i++)
            free(game->transition[i]);
        free(game->transition);
    }
}
