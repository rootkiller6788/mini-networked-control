#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include "packet_loss_estimator.h"
#include "packet_loss_analysis.h"
#include "packet_loss_predictor.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=========================================\n");
    printf("  Packet Loss in Networked Control\n");
    printf("  Demo — All Channel Models\n");
    printf("=========================================\n\n");

    /* Demo all channel models */
    printf("--- Bernoulli Channel ---\n");
    BernoulliChannel* bc = pl_bernoulli_create(0.25, 42);
    for (int i = 0; i < 1000; i++) pl_bernoulli_transmit(bc);
    printf("  Target p=0.25, empirical=%.4f\n", pl_bernoulli_get_loss_rate(bc));
    pl_bernoulli_free(bc);

    printf("\n--- Gilbert-Elliott Channel ---\n");
    GilbertElliottChannel* ge = pl_gilbert_elliott_create(0.05, 0.30, 0.02, 0.50, 42);
    pl_gilbert_elliott_compute_steady_state(ge);
    printf("  Steady-state loss rate: %.4f\n", ge->steady_loss_rate);
    printf("  Burstiness: %.4f\n", pl_gilbert_elliott_burstiness(ge));
    for (int i = 0; i < 500; i++) pl_gilbert_elliott_transmit(ge);
    printf("  Empirical loss: %.4f\n", ge->observed_loss_rate);
    pl_gilbert_elliott_free(ge);

    printf("\n--- LQR Controller Demo ---\n");
    LTISystem* sys = pl_lti_create(2, 1, 1, 2);
    double A[4] = {1.0, 0.1, 0.0, 0.9};
    double B[2] = {0.0, 0.1};
    pl_lti_set_A(sys, A, 2); pl_lti_set_B(sys, B);
    for (int i = 0; i < 4; i++) sys->Q[i] = (i % 3 == 0) ? 1.0 : 0.0;
    double Rc[1] = {0.1};
    LQRSolution* sol = pl_lqr_solve(sys, Rc, 500, 1e-6);
    printf("  LQR: %d iterations, rho_cl=%.6f\n",
           sol->iterations, pl_lqr_spectral_radius(sol, sys));
    printf("  L = [%.4f, %.4f]\n", sol->L[0], sol->L[1]);
    pl_lqr_free(sol);
    pl_lti_free(sys);

    printf("\n--- Intermittent Kalman Filter ---\n");
    double Ak[4] = {1.0, 0.1, 0.0, 0.9};
    double Ck[2] = {1.0, 0.0};
    double Qk[4] = {0.01, 0.0, 0.0, 0.01};
    double Rk[1] = {0.01};
    double x0[2] = {0.0, 0.0}, P0[4] = {1.0, 0.0, 0.0, 1.0};
    IntermittentKalmanFilter* ikf = pl_ikf_create(Ak, Ck, Qk, Rk, 2, 1, x0, P0, 0.75);
    printf("  gamma_c ∈ [%.4f, %.4f]\n", ikf->critical_gamma_lower, ikf->critical_gamma_upper);
    printf("  Mean-stable: %s\n", pl_ikf_is_stable(ikf) ? "YES" : "NO");
    pl_ikf_free(ikf);

    printf("\n--- Critical Probability ---\n");
    double Al[4] = {1.2, 0.0, 0.0, 1.1};
    double Bl[2] = {1.0, 0.5};
    double Cl[2] = {1.0, 0.0};
    double Ll[2] = {0.5, 0.3};
    CriticalProbabilityAnalysis* cpa = pl_critical_prob_analyze(Al, Bl, Cl, Ll, 2, 1, 1);
    printf("  p_c(sensor)=%.4f, p_c(actuator)=%.4f, p_c(joint)=%.4f\n",
           cpa->p_c_sensor, cpa->p_c_actuator, cpa->p_c_joint);
    pl_critical_prob_free(cpa);

    printf("\n--- PPC Buffer ---\n");
    double Ap[4] = {1.0, 0.1, 0.0, 0.9};
    double Bp[2] = {0.0, 0.1};
    double Lp[2] = {2.0, 0.5};
    double x0p[2] = {1.0, 0.0};
    PacketizedPredictiveControl* ppc = pl_ppc_create(5, Ap, Bp, Lp, 2, 1, x0p);
    printf("  Horizon=5, buffer[0]=[%.4f]\n", ppc->buffer[0][0]);
    pl_ppc_free(ppc);

    printf("\n=========================================\n");
    printf("  Demo complete.\n");
    printf("=========================================\n");
    return 0;
}
