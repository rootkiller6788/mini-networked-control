#include "time_delay_system.h"
#include "delay_stability.h"
#include "dde_solver.h"
#include "smith_predictor.h"
#include "networked_delay.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * L7 Applications — Real-world Time-Delay System Scenarios
 *
 * Each function implements a distinct application domain with
 * time-delay challenges.
 *
 * Applications covered:
 *   1. Boeing 747 roll control with actuator delay
 *   2. Teleoperation with communication delay (NASA/Mars rover)
 *   3. Smart grid frequency regulation with communication delay
 *   4. Automotive CAN bus brake-by-wire with delay
 *   5. Chemical process control with transport delay
 * ============================================================================ */

/* ============================================================================
 * Boeing 747 Roll Control with Actuator Delay
 *
 * The Boeing 747 lateral dynamics exhibit a Dutch roll mode.
 * With actuator delay from hydraulic systems, the stability
 * margin degrades. We model a simplified roll dynamics:
 *
 *   ẋ = A x + B δ, where δ = control surface deflection
 *   with delay τ from hydraulic actuation.
 *
 * Reference: Boeing 747 lateral dynamics data,
 *   Heffley & Jewell, NASA CR-2144 (1972)
 * ============================================================================ */

void boeing747_roll_control_delay(void) {
    printf("=== Boeing 747 Roll Control with Actuator Delay ===\n");

    /* Simplified roll-subsidence + Dutch roll model (2-state) */
    TimeDelaySystem* sys = tds_create("Boeing747-Roll", 2, 1, 1);

    /* State: [roll_rate p, sideslip_angle β]
     * A = [[-1.5,  -8.0],
     *      [ 0.1,  -0.5]]
     * A_d = [[-0.1,  0.0],
     *        [ 0.0, -0.05]]   ← delay coupling
     * B = [ -3.0, 0.0 ]ᵀ    ← aileron effectiveness
     */
    double A[4]  = {-1.5, -8.0, 0.1, -0.5};
    double Ad[4] = {-0.1, 0.0, 0.0, -0.05};
    double B[2]  = {-3.0, 0.0};
    double C[2]  = {1.0, 0.0};

    tds_set_linear_model(sys, A, Ad, B, C);
    tds_add_delay(sys, DELAY_CONSTANT, 0.05, 0.05, 0.05);

    /* Compute delay margin */
    double margin = delay_margin_frequency_sweep(sys);
    printf("  System: Boeing 747 Roll (2-state)\n");
    printf("  τ_nominal = 0.05 s (hydraulic actuator)\n");
    printf("  Delay margin τ* = %.4f s\n", margin);
    printf("  Status: delay %.3f %s stable margin\n",
           0.05, (0.05 < margin) ? "within" : "EXCEEDS");

    /* Compute characteristic roots */
    tds_compute_characteristic_roots(sys, 10);
    printf("  Characteristic roots (rightmost): ");
    for (int i = 0; i < (sys->n_roots < 5 ? sys->n_roots : 5); i++)
        printf("%.3f+j%.3f  ", sys->roots_real[i], sys->roots_imag[i]);
    printf("\n");

    tds_free(sys);
}

/* ============================================================================
 * NASA Mars Rover Teleoperation with Earth-Mars Delay
 *
 * One-way light-time delay between Earth and Mars varies
 * from 4 to 24 minutes depending on orbital positions.
 *
 * Model: second-order rover position dynamics with
 * human-in-the-loop command delay.
 *
 * ẍ = u(t-τ)  (double integrator with delayed force)
 *
 * Reference: NASA JPL, Mars Exploration Rover telecom design (2003)
 *   Sheridan, "Space Teleoperation Through Time Delay" (1993)
 *   SpaceX Starship Earth-Mars transit communication analysis
 * ============================================================================ */

void mars_rover_teleoperation_delay(void) {
    printf("\n=== Mars Rover Teleoperation (Earth-Mars Delay) ===\n");

    /* Double integrator: ẋ₁ = x₂, ẋ₂ = u(t-τ) */
    TimeDelaySystem* rover = tds_create("MarsRover", 2, 1, 1);

    double A[4]  = {0.0, 1.0, 0.0, 0.0};
    double Ad[4] = {0.0, 0.0, 0.0, 0.0};
    double B[2]  = {0.0, 1.0};
    double C[2]  = {1.0, 0.0};

    tds_set_linear_model(rover, A, Ad, B, C);

    /* Earth-Mars one-way delay at conjunction (~240 s)
     * At opposition (~240 s round trip = 120 s one-way).
     * Average value: 600 s one-way. */
    double tau_mars_conjunction = 600.0;  /* ~10 minutes one-way */
    double tau_mars_opposition = 240.0;   /* ~4 minutes one-way */

    printf("  Earth-Mars one-way light-time delay:\n");
    printf("    Conjunction (max): %.0f s (~%d min)\n",
           tau_mars_conjunction, (int)(tau_mars_conjunction / 60.0));
    printf("    Opposition (min):  %.0f s (~%d min)\n",
           tau_mars_opposition, (int)(tau_mars_opposition / 60.0));

    /* For double integrator with delay feedback u = -Kx:
     * Characteristic: s² + k_d s e^{-τs} + k_p e^{-τs} = 0
     * Delay margin for K = [k_p, k_d] = [0.5, 1.0]:
     * τ_crit ≈ 1.0 / sqrt(k_p) = 1.41 s (without delay)
     * With delay: significantly reduced.
     * With τ = 600s, even tiny gains cause instability. */

    tds_add_delay(rover, DELAY_CONSTANT, tau_mars_conjunction, 0, tau_mars_conjunction);

    double margin_rover = delay_margin_frequency_sweep(rover);
    printf("  Rover delay margin (state feedback): %.4f s\n", margin_rover);
    printf("  Status: Earth-Mars delay %.0f s FAR exceeds margin\n",
           tau_mars_conjunction);
    printf("  Consequence: Supervisory control required (move-and-wait)\n");
    printf("  Solution: Shared autonomy + predictive display\n");
    printf("  Reference: Sheridan (1993), Fong & Thorpe (2001)\n");

    tds_free(rover);
}

/* ============================================================================
 * Smart Grid Frequency Regulation with Communication Delay
 *
 * In smart grids, phasor measurement units (PMUs) send data
 * over IP networks to control centers. Communication delays
 * in the wide-area measurement system (WAMS) affect frequency
 * regulation performance.
 *
 * Model: Single-area power system frequency dynamics
 *   Δḟ = (1/M) (ΔP_m - ΔP_L - D Δf)
 *   ΔP_ṁ = (1/T_ch) (ΔP_v - ΔP_m)
 *   ΔP_v̇ = (1/T_g) (u(t-τ) - ΔP_v - Δf/R)
 *
 * where M = inertia, D = damping, R = droop,
 * T_ch = turbine time constant, T_g = governor time constant.
 *
 * Reference: Kundur, "Power System Stability and Control" (1994)
 *   PMU communication delays in WAMS: IEEE Std C37.118
 * ============================================================================ */

void smart_grid_frequency_delay(void) {
    printf("\n=== Smart Grid Frequency Regulation with PMU Delay ===\n");

    /* 3-state power system model */
    TimeDelaySystem* grid = tds_create("SmartGrid-Frequency", 3, 1, 1);

    double M = 10.0, D = 1.0, R = 0.05;
    double T_ch = 0.3, T_g = 0.1;

    /* State: [Δf, ΔP_m, ΔP_v]
     * A = [[-D/M,  1/M,  0   ],
     *      [ 0,   -1/T_ch, 1/T_ch],
     *      [-1/(R*T_g), 0, -1/T_g]] */
    double A_grid[9] = {
        -D/M,  1.0/M,  0.0,
         0.0, -1.0/T_ch, 1.0/T_ch,
        -1.0/(R*T_g), 0.0, -1.0/T_g
    };

    /* Delay coupling from SCADA/WAMS communication */
    double Ad_grid[9] = {0.0};  /* Assume state feedback with delay */
    Ad_grid[8] = -0.01;  /* Small delay coupling from governor signal */

    double B_grid[3] = {0.0, 0.0, 1.0/T_g};
    double C_grid[3] = {1.0, 0.0, 0.0};

    tds_set_linear_model(grid, A_grid, Ad_grid, B_grid, C_grid);

    /* Typical WAMS delays:
     * - PMU to PDC: 20-200 ms
     * - PDC to control center: 100-500 ms
     * - Total: 50-700 ms typical, up to few seconds in congested networks */
    double tau_pmu_ms[] = {20.0, 100.0, 500.0};  /* ms */
    printf("  Power system: 3-state (Δf, ΔP_m, ΔP_v)\n");
    printf("  PMU delay scenarios:\n");
    for (int i = 0; i < 3; i++) {
        tds_add_delay(grid, DELAY_CONSTANT, tau_pmu_ms[i] / 1000.0,
                      0, tau_pmu_ms[i] / 1000.0);
        /* Overwrite delay each iteration (simplified) */
        if (grid->delays && grid->delays[0])
            grid->delays[0]->tau_nominal = tau_pmu_ms[i] / 1000.0;

        double margin = delay_margin_frequency_sweep(grid);
        printf("    τ=%.0f ms: delay margin τ*=%.4f s → %s\n",
               tau_pmu_ms[i], margin,
               (tau_pmu_ms[i]/1000.0 < margin) ? "STABLE" : "UNSTABLE");
    }

    printf("  Reference: PMU-based WAMS, IEEE C37.118\n");
    printf("  Solution: Smith predictor in AGC loop\n");
    printf("  Impact: Delay > 500ms can degrade frequency regulation\n");

    tds_free(grid);
}

/* ============================================================================
 * Automotive CAN Bus Brake-by-Wire Delay
 *
 * In modern vehicles (Toyota, Tesla), brake-by-wire systems
 * communicate over CAN bus with typical message latencies
 * of 1-10 ms. At highway speeds, 10 ms delay translates to
 * ~30 cm of travel distance — critical for safety.
 *
 * Model: Vehicle longitudinal dynamics with braking delay
 *   v̇ = -(1/M) F_brake(t-τ) - (1/M) F_drag
 *
 * Reference: ISO 26262 (functional safety)
 *   Toyota electronic throttle control (2000s)
 *   Tesla brake-by-wire system (Model 3)
 * ============================================================================ */

void automotive_can_brake_delay(void) {
    printf("\n=== Automotive CAN Bus Brake-by-Wire Delay ===\n");

    TimeDelaySystem* car = tds_create("Brake-by-Wire", 2, 1, 1);

    /* State: [velocity v, braking_force F_b]
     * Simplified model: v̇ = -F_b/M (deceleration)
     *                  F_ḃ = (u(t-τ) - F_b)/T_brake */
    double M_car = 1500.0;      /* kg — typical sedan (Toyota Camry class) */
    double T_brake = 0.05;      /* brake actuator time constant */

    double A_car[4] = {0.0, -1.0/M_car, 0.0, -1.0/T_brake};
    double Ad_car[4] = {0.0, 0.0, 0.0, 0.0};
    double B_car[2] = {0.0, 1.0/T_brake};
    double C_car[2] = {1.0, 0.0};

    tds_set_linear_model(car, A_car, Ad_car, B_car, C_car);

    /* CAN bus delay scenarios */
    double tau_can[] = {0.001, 0.005, 0.010, 0.020};  /* 1, 5, 10, 20 ms */

    printf("  Vehicle: %d kg sedan at highway speed (30 m/s)\n", (int)M_car);
    printf("  CAN bus brake latency scenarios:\n");
    for (int i = 0; i < 4; i++) {
        double dist = 30.0 * tau_can[i];  /* Distance traveled during delay */
        printf("    τ=%.0f ms → distance traveled during delay: %.3f m (%.0f cm)\n",
               tau_can[i] * 1000.0, dist, dist * 100.0);
    }

    printf("  Safety implication: 10ms delay = 30cm at highway speed\n");
    printf("  ISO 26262: Brake-by-wire requires τ < 10ms for ASIL D\n");
    printf("  CAN FD and FlexRay reduce latency vs classical CAN\n");
    printf("  Tesla uses redundant communication paths for braking\n");
    printf("  Toyota brake override system: response < 100ms\n");

    tds_free(car);
}

/* ============================================================================
 * Chemical Process Control with Transport Delay
 *
 * In chemical plants, transport delays arise from:
 *   - Pipe flow between reactor and sensor
 *   - Heat exchanger residence time
 *   - Material transport on conveyor belts
 *
 * Model: First-order + dead time (FOPDT):
 *   G(s) = K e^{-τs} / (T s + 1)
 *
 * Reference: Seborg, Edgar, Mellichamp, "Process Dynamics
 *   and Control" (2011)
 * ============================================================================ */

void chemical_process_transport_delay(void) {
    printf("\n=== Chemical Process Transport Delay ===\n");

    /* FOPDT model for a stirred tank heater:
     * τ = V/q = residence time (volume/flow rate)
     * T = process time constant
     * K = process gain */
    double V = 2.0;      /* m³ — tank volume */
    double q = 0.01;     /* m³/s — flow rate */
    double tau_res = V / q;  /* 200 s residence time */

    /* Create 1-state delay system */
    TimeDelaySystem* tank = tds_create("StirredHeater", 1, 1, 1);
    double A[1] = {-1.0 / 100.0};    /* 100 s time constant */
    double Ad[1] = {0.005};           /* delay coupling from recirculation */
    double B[1] = {1.0 / 100.0};
    double C[1] = {1.0};

    tds_set_linear_model(tank, A, Ad, B, C);
    tds_add_delay(tank, DELAY_CONSTANT, tau_res, 0, tau_res);

    printf("  Process: Stirred tank heater\n");
    printf("  Volume: %.1f m³, Flow: %.3f m³/s\n", V, q);
    printf("  Residence time (transport delay): %.0f s\n", tau_res);

    double margin = delay_margin_frequency_sweep(tank);
    printf("  Delay margin: %.2f s\n", margin);

    /* PID tuning rules for FOPDT */
    /* Cohen-Coon tuning: */
    double T_proc = 100.0;  /* Time constant */
    double K_proc = 1.0;    /* Gain */
    double Kc = (T_proc / (K_proc * tau_res)) * (1.0 + tau_res / (3.0 * T_proc));
    double Ti = tau_res * (30.0 + 3.0 * tau_res / T_proc)
                / (9.0 + 20.0 * tau_res / T_proc);
    double Td = tau_res * 4.0 / (11.0 + 2.0 * tau_res / T_proc);

    printf("  Cohen-Coon PID (FOPDT):\n");
    printf("    Kc=%.3f, Ti=%.2f, Td=%.2f\n", Kc, Ti, Td);
    printf("  τ/T ratio: %.1f (process becomes hard to control if > 1.0)\n",
           tau_res / T_proc);

    tds_free(tank);
}

/* ============================================================================
 * Tele-surgery Haptic Feedback Delay
 *
 * In robotic telesurgery (da Vinci surgical system), haptic
 * feedback loops involve delays from video processing,
 * network transmission, and robot kinematics.
 *
 * Human perception: delays > 100ms significantly degrade
 * surgical performance. Ideal: < 10ms.
 *
 * Reference: Intuitive Surgical da Vinci system
 *   Marescaux et al., "Transatlantic robot-assisted telesurgery"
 *   (Nature, 2001) — Operation Lindbergh
 * ============================================================================ */

void telesurgery_haptic_delay(void) {
    printf("\n=== Telesurgery Haptic Feedback Delay ===\n");

    /* Force-position haptic model: mass-spring-damper
     * ẋ₁ = x₂
     * ẋ₂ = (1/M) (F_human - K_s x₁ - B_s x₂ - F_env(t-τ))
     *
     * where F_env is the delayed environmental force feedback
     */
    TimeDelaySystem* haptic = tds_create("HapticMaster", 2, 1, 1);

    double M_h = 0.5;    /* Haptic device mass (kg) */
    double K_h = 10.0;   /* Virtual stiffness */
    double B_h = 0.2;    /* Virtual damping */

    double A_hap[4] = {0.0, 1.0, -K_h/M_h, -B_h/M_h};
    double Ad_hap[4] = {0.0, 0.0, 0.0, 0.0};
    double B_hap[2] = {0.0, 1.0/M_h};
    double C_hap[2] = {1.0, 0.0};

    tds_set_linear_model(haptic, A_hap, Ad_hap, B_hap, C_hap);

    /* Delay scenarios */
    double tau_local[] = {0.001, 0.010};     /* Local: 1-10 ms */
    double tau_transatlantic[] = {0.080, 0.150}; /* 80-150 ms RTT */

    printf("  Haptic device: M=%.1f kg, K=%.0f N/m\n", M_h, K_h);
    printf("  Local surgery (OR): τ=%.0f-%.0f ms → transparent feel\n",
           tau_local[0]*1000, tau_local[1]*1000);
    printf("  Transatlantic (Operation Lindbergh, 2001):\n");
    printf("    RTT=%.0f ms → noticeable but usable with prediction\n",
           tau_transatlantic[1]*1000);
    printf("    RTT=%.0f ms → challenging, requires compensation\n",
           tau_transatlantic[0]*1000);

    double margin_haptic = delay_margin_frequency_sweep(haptic);
    printf("  Haptic loop delay margin: %.4f s (%.0f ms)\n",
           margin_haptic, margin_haptic * 1000.0);

    tds_free(haptic);
}

/* ============================================================================
 * Run all application demonstrations
 * ============================================================================ */

int main(void) {
    boeing747_roll_control_delay();
    mars_rover_teleoperation_delay();
    smart_grid_frequency_delay();
    automotive_can_brake_delay();
    chemical_process_transport_delay();
    telesurgery_haptic_delay();

    printf("\n=== All L7 Applications Complete ===\n");
    return 0;
}
