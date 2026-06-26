#include "cps_security_core.h"
#include "cps_gametheory.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 3: Game-Theoretic Defense Strategy for CPS ===\n\n");

    /* Defender has 3 strategies: normal operation, heightened alert, lockdown
     * Attacker has 3 strategies: no attack, stealthy attack, surge attack
     * Payoff matrix = attacker's expected damage (defender minimizes) */
    CPSZeroSumGame game;
    cps_zerosum_game_init(&game, 3, 3);
    double* P = game.payoff_matrix;
    /*       NoAtk  Stealth  Surge */
    /* Norm*/ P[0]=0; P[1]=8; P[2]=15;
    /* Alrt*/ P[3]=1; P[4]=5; P[5]=8;
    /* Lock*/ P[6]=3; P[7]=4; P[8]=6;

    printf("Payoff matrix (attacker damage / defender cost):\n");
    printf("           No Attack  Stealth   Surge\n");
    printf(" Normal    %9.1f %9.1f %9.1f\n", P[0], P[1], P[2]);
    printf(" Alert     %9.1f %9.1f %9.1f\n", P[3], P[4], P[5]);
    printf(" Lockdown  %9.1f %9.1f %9.1f\n", P[6], P[7], P[8]);

    int ret = cps_solve_zerosum_lp(&game, 2000);
    printf("\nSolution via fictitious play (%d iterations):\n", 2000);
    printf("Saddle value (expected damage): %.3f\n", game.saddle_value);
    printf("Defender mixed strategy: [%.3f, %.3f, %.3f]\n",
           game.defender_mixed[0], game.defender_mixed[1],
           game.defender_mixed[2]);
    printf("Attacker mixed strategy:  [%.3f, %.3f, %.3f]\n",
           game.attacker_mixed[0], game.attacker_mixed[1],
           game.attacker_mixed[2]);
    printf("Game solved: %s\n", game.solved ? "YES" : "NO");

    /* Stackelberg equilibrium */
    double def_sse[3], att_sse[3], sse_val;
    ret = cps_stackelberg_sse(&game, def_sse, att_sse, &sse_val);
    printf("\nStackelberg Equilibrium (defender leads):\n");
    printf("SSE value: %.3f\n", sse_val);
    printf("Defender SSE strategy: [%.1f, %.1f, %.1f]\n",
           def_sse[0], def_sse[1], def_sse[2]);
    printf("Attacker best response: [%.1f, %.1f, %.1f]\n",
           att_sse[0], att_sse[1], att_sse[2]);

    /* Optimal detection threshold */
    double opt_thresh = cps_optimal_detection_threshold(
        1000.0, 50000.0, 100.0, 0.01, 3.0);
    printf("\nOptimal detection threshold: %.2f\n", opt_thresh);
    printf("(Balances false positive cost, missed detection cost, and investment)\n");

    cps_zerosum_game_free(&game);
    return 0;
}
