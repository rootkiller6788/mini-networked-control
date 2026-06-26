#include "ebc_core.h"
#include "ebc_trigger.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
 * example2_trigger_comparison.c -- Trigger Type Comparison
 *
 * Compares different event-triggering conditions on the same system:
 *   1. Absolute error trigger
 *   2. Relative error trigger
 *   3. Mixed threshold trigger
 *
 * Shows how different trigger types affect the number of events
 * and control performance.
 *
 * L6: Canonical problem -- comparison of triggering strategies
 */

int main(void) {
    printf("=== Example 2: Trigger Type Comparison ===\n\n");

    int n = 2, m = 1;
    double x0[] = {2.0, 1.0};

    /* Test state evolution (simulated as decaying oscillation) */
    double x_vals[] = {2.0, 1.0, 1.5, 0.5, 1.0, 0.2, 0.5, -0.1, 0.1, -0.2, 0.0, 0.0};

    double x_last[] = {2.0, 1.0};

    /* Create system for each trigger test */
    EBC_System* sys = ebc_system_create(n, m, EBC_CONTINUOUS_ETC);
    ebc_system_set_state(sys, x0);

    /* Absolute trigger */
    EBC_TriggerParams tp_abs = ebc_trigger_default();
    tp_abs.type = EBC_ABSOLUTE_ERROR;
    tp_abs.epsilon = 0.5;

    /* Relative trigger */
    EBC_TriggerParams tp_rel = ebc_trigger_make_relative(0.3);

    /* Mixed trigger */
    EBC_TriggerParams tp_mix = ebc_trigger_mixed(0.2, 0.1);

    printf("State evolution with 3 trigger types:\n");
    printf("  t      x1      x2      Abs?  Rel?  Mix?\n");
    printf("  -----  ------  ------  ----  ----  ----\n");

    int event_abs = 0, event_rel = 0, event_mix = 0;

    /* Simulate through state sequence */
    for (int k = 0; k < 6; k++) {
        sys->x[0] = x_vals[k * 2];
        sys->x[1] = x_vals[k * 2 + 1];
        sys->x_last[0] = x_last[0];
        sys->x_last[1] = x_last[1];

        bool abs_fire = ebc_trigger_absolute(sys, &tp_abs);
        bool rel_fire = ebc_trigger_relative(sys, &tp_rel);
        bool mix_fire = ebc_trigger_mixed_threshold(sys, &tp_mix);

        if (abs_fire) { event_abs++; memcpy(x_last, sys->x, n * sizeof(double)); }
        /* For relative and mixed, track independently */
        /* Simplified: just record whether they would fire */

        printf("  %.1f   %.4f  %.4f   %s    %s    %s\n",
               k * 0.5, sys->x[0], sys->x[1],
               abs_fire ? "YES" : "no ",
               rel_fire ? "YES" : "no ",
               mix_fire ? "YES" : "no ");
    }

    printf("\nAbsolute threshold triggered %d times\n", event_abs);

    ebc_system_free(sys);
    printf("\nExample 2 complete.\n");
    return 0;
}
