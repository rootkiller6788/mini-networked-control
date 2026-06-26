#ifndef CPS_GAMETHEORY_H
#define CPS_GAMETHEORY_H

#include "cps_security_core.h"

/* ============================================================================
 * Game-Theoretic CPS Security (L8: Advanced Topics)
 *
 * Reference: Zhu & Basar (2015) "Game-Theoretic Methods for Robustness,
 *            Security, and Resilience of Cyberphysical Control Systems"
 *
 * Models the attacker-defender interaction as a dynamic game:
 *   - Zero-sum game: attacker minimizes cost, defender maximizes
 *   - Stackelberg game: defender commits first (leader), attacker follows
 *   - Stochastic game: state transitions depend on both players
 * ============================================================================ */

/* --- Player Types --- */

typedef struct {
    double* strategy;           /* Mixed/probability strategy vector */
    int n_actions;              /* Number of pure actions */
    double payoff;              /* Expected payoff */
    int is_attacker;            /* 1 = attacker, 0 = defender */
    double resource_budget;     /* Total attack/defense resource */
    double* action_costs;       /* Cost per action */
} CPSGamePlayer;

/* Zero-sum game matrix */
typedef struct {
    double* payoff_matrix;      /* Attacker's payoff matrix (flattened) */
    int n_defender_actions;     /* Defender's pure strategies */
    int n_attacker_actions;     /* Attacker's pure strategies */
    double saddle_value;        /* Value of the game */
    double* defender_mixed;     /* Optimal mixed strategy for defender */
    double* attacker_mixed;     /* Optimal mixed strategy for attacker */
    int solved;                 /* Whether saddle point has been found */
} CPSZeroSumGame;

/* Dynamic (multi-stage) security game */
typedef struct {
    CPSZeroSumGame* stage_games;/* Game at each stage */
    int n_stages;
    double discount_factor;     /* Future payoff discount */
    double* state_values;       /* Value function */
    int* defender_policy;       /* Optimal policy per state */
    int* attacker_policy;
    int n_states;               /* Number of system states */
    double** transition;        /* State transition probs */
} CPSDynamicGame;

/* --- Zero-Sum Game Solving (L5: Algorithms) --- */

/* Solve a 2x2 zero-sum game analytically
 * Player 1 (defender) chooses row, Player 2 (attacker) chooses column
 * payoff = attacker's payoff (defender minimizes this) */
void cps_solve_2x2_zerosum(double a11, double a12, double a21, double a22,
                            double* saddle_value,
                            double* defender_p1,
                            double* attacker_q1);

/* Solve general zero-sum game via linear programming
 * Uses the simplex method for the minimax LP:
 * max v s.t. sum_i x_i * A_{ij} >= v for all j, sum_i x_i = 1, x_i >= 0 */
int cps_solve_zerosum_lp(CPSZeroSumGame* game, int max_iter);

/* Compute Nash equilibrium of a zero-sum game
 * Returns 0 on success, -1 if no saddle point found */
int cps_zerosum_nash_equilibrium(CPSZeroSumGame* game);

/* --- Stackelberg Security Game (L8: Advanced) --- */

/* Strong Stackelberg Equilibrium (SSE)
 * Defender (leader) commits to a mixed strategy;
 * Attacker (follower) best-responds.
 * Computed via multiple LPs over attacker's best-response regions. */
int cps_stackelberg_sse(CPSZeroSumGame* game,
                         double* defender_strategy,
                         double* attacker_best_response,
                         double* sse_value);

/* Defender utility under attacker's best response */
double cps_defender_utility_breach(CPSZeroSumGame* game,
                                    const double* defender_strategy);

/* --- Multistage Dynamic Game (L8: Advanced) --- */

/* Initialize a dynamic security game */
void cps_dynamic_game_init(CPSDynamicGame* game, int n_states, int n_stages);

/* Solve via backward induction (finite horizon)
 * V_k(s) = max_{a_d} min_{a_a} [ r(s, a_d, a_a) + gamma * V_{k+1}(s') ] */
void cps_dynamic_game_solve(CPSDynamicGame* game);

/* Get optimal action for state s at stage t */
void cps_dynamic_game_policy(CPSDynamicGame* game, int state, int stage,
                              int* def_action, int* att_action);

/* --- Attack-Defense Cost Models --- */

/* LQG cost under attack:
 * J = E[ sum (x_k'*Q*x_k + u_k'*R*u_k) ] under attack distribution */
double cps_lqg_cost_under_attack(const double* A, const double* B,
                                  const double* Q, const double* R,
                                  const double* Sigma_a, int n, int m);

/* Defender's security investment optimization:
 * Choose detection threshold to minimize:
 *   C_total = C_fp * P_fp + C_fn * P_fn + C_invest(thresh)
 * where C_fp = cost of false positive (unnecessary response)
 *       C_fn = cost of false negative (missed attack) */
double cps_optimal_detection_threshold(double cost_fp, double cost_fn,
                                        double cost_invest_coeff,
                                        double attack_prior,
                                        double detection_snr);

/* --- Stochastic Game with Imperfect Information --- */

/* Bayesian update of attacker type belief
 * Given observed actions, update belief about attacker capability */
void cps_bayesian_attacker_belief(double* belief, int n_types,
                                   const double* observation_likelihood,
                                   int observed_action);

/* Defender's optimal strategy under type uncertainty */
void cps_defender_robust_strategy(const CPSZeroSumGame** games,
                                   int n_types, const double* type_belief,
                                   double* robust_strategy);

/* --- Utilities --- */

void cps_gameplayer_init(CPSGamePlayer* player, int n_actions, int is_attacker);
void cps_gameplayer_free(CPSGamePlayer* player);
void cps_zerosum_game_init(CPSZeroSumGame* game, int n_def, int n_att);
void cps_zerosum_game_free(CPSZeroSumGame* game);
void cps_dynamic_game_free(CPSDynamicGame* game);

#endif /* CPS_GAMETHEORY_H */
