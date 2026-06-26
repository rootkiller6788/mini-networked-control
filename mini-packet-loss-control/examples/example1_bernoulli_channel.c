#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 1: Bernoulli Packet Loss Channel ===\n\n");

    /* Create a Bernoulli channel with 20% loss rate */
    BernoulliChannel* ch = pl_bernoulli_create(0.20, 12345);
    printf("Bernoulli channel: p=%.2f\n", ch->loss_probability);

    /* Simulate 1000 transmissions */
    PacketHistory* hist = pl_history_create(2000);
    int received = 0, lost = 0;
    for (int i = 0; i < 1000; i++) {
        PacketStatus st = pl_bernoulli_transmit(ch);
        if (st == PACKET_RECEIVED) received++;
        else lost++;
        pl_history_record(hist, (unsigned long)i, i * 0.001, st);
    }

    printf("Results: %d received, %d lost (empirical rate=%.3f, target=0.20)\n",
           received, lost, (double)lost / 1000.0);

    /* Stability check for different system spectral radii */
    printf("\nStability analysis (Sinopoli 2004, Theorem 2):\n");
    double rhos[] = {0.5, 1.0, 1.2, 1.5, 2.0};
    for (int i = 0; i < 5; i++) {
        bool stable = pl_bernoulli_is_stable(0.20, rhos[i]);
        double p_c = (rhos[i] > 1.0) ? (1.0 - 1.0/(rhos[i]*rhos[i])) : 1.0;
        printf("  rho(A)=%.1f: critical_p=%.3f, p=0.20 -> %s\n",
               rhos[i], p_c, stable ? "STABLE" : "UNSTABLE");
    }

    /* Show loss statistics */
    printf("\nLoss statistics:\n");
    printf("  Recent loss rate: %.4f\n", pl_history_loss_rate(hist, 100));
    printf("  Burst index: %.4f (expect ~0.20 for Bernoulli)\n",
           pl_history_burst_index(hist));

    pl_history_free(hist);
    pl_bernoulli_free(ch);
    printf("\nExample 1 complete.\n");
    return 0;
}
