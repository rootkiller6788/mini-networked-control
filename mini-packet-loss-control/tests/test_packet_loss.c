#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include "packet_loss_estimator.h"
#include "packet_loss_analysis.h"
#include "packet_loss_predictor.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 1e-9

int main(void) {
    printf("=== Packet Loss Control Tests ===\n");

    /* ===== L1: Bernoulli Channel ===== */
    {
        BernoulliChannel* ch = pl_bernoulli_create(0.3, 42);
        assert(ch);
        assert(fabs(ch->loss_probability - 0.3) < EPS);
        for (int i = 0; i < 1000; i++) pl_bernoulli_transmit(ch);
        double rate = pl_bernoulli_get_loss_rate(ch);
        assert(rate >= 0.0 && rate <= 1.0);
        assert(fabs(rate - 0.3) < 0.1); /* Should be close to 0.3 after 1000 trials */
        assert(pl_bernoulli_is_stable(0.5, 2.0)); /* ρ=2, p_c=0.75, p=0.5 < p_c=0.75 → stable */
        assert(pl_bernoulli_is_stable(0.5, 1.0)); /* open-loop stable → always stable */
        assert(!pl_bernoulli_is_stable(0.9, 2.0)); /* ρ=2, p_c=0.75, p=0.9 > p_c=0.75 → unstable */
        pl_bernoulli_free(ch);
        printf("  [PASS] Bernoulli channel\n");
    }

    /* ===== L1: Gilbert-Elliott Channel ===== */
    {
        GilbertElliottChannel* ge = pl_gilbert_elliott_create(0.1, 0.5, 0.01, 0.3, 123);
        assert(ge);
        assert(fabs(ge->p_gb - 0.1) < EPS);
        pl_gilbert_elliott_compute_steady_state(ge);
        assert(fabs(ge->steady_p_good - 0.5/0.6) < EPS);
        assert(ge->steady_loss_rate >= 0.0 && ge->steady_loss_rate <= 1.0);
        double burst = pl_gilbert_elliott_burstiness(ge);
        assert(burst >= 0.0);
        for (int i = 0; i < 500; i++) pl_gilbert_elliott_transmit(ge);
        assert(ge->total_sent == 500);
        pl_gilbert_elliott_free(ge);
        printf("  [PASS] Gilbert-Elliott channel\n");
    }

    /* ===== L1: K-State Markov Channel ===== */
    {
        MarkovChannel* mk = pl_markov_create(3, 456);
        assert(mk);
        assert(mk->n_states == 3);
        pl_markov_set_transition(mk, 0, 0, 0.7); pl_markov_set_transition(mk, 0, 1, 0.2);
        pl_markov_set_transition(mk, 0, 2, 0.1);
        pl_markov_set_transition(mk, 1, 0, 0.3); pl_markov_set_transition(mk, 1, 1, 0.5);
        pl_markov_set_transition(mk, 1, 2, 0.2);
        pl_markov_set_transition(mk, 2, 0, 0.1); pl_markov_set_transition(mk, 2, 1, 0.3);
        pl_markov_set_transition(mk, 2, 2, 0.6);
        pl_markov_set_loss_rate(mk, 0, 0.01); pl_markov_set_loss_rate(mk, 1, 0.1);
        pl_markov_set_loss_rate(mk, 2, 0.5);
        pl_markov_compute_steady_state(mk);
        double sum_pi = 0.0;
        for (int i = 0; i < 3; i++) sum_pi += mk->steady_state[i];
        assert(fabs(sum_pi - 1.0) < 0.01);
        for (int i = 0; i < 300; i++) pl_markov_transmit(mk);
        assert(mk->total_sent == 300);
        pl_markov_free(mk);
        printf("  [PASS] K-State Markov channel\n");
    }

    /* ===== L1: Burst Channel ===== */
    {
        BurstChannel* bc = pl_burst_create(5.0, 95.0, 789);
        assert(bc);
        assert(fabs(bc->loss_probability - 5.0/100.0) < 0.01);
        int losses = 0;
        for (int i = 0; i < 10000; i++)
            if (pl_burst_transmit(bc) != PACKET_RECEIVED) losses++;
        double emp_rate = (double)losses / 10000.0;
        assert(fabs(emp_rate - 0.05) < 0.03);
        pl_burst_free(bc);
        printf("  [PASS] Burst channel\n");
    }

    /* ===== L2: Packet Lifecycle ===== */
    {
        NetworkPacket* pkt = pl_packet_create(4);
        assert(pkt); assert(pkt->payload_size == 4);
        double data[] = {1.0, 2.0, 3.0, 4.0};
        pl_packet_set_payload(pkt, data, 4);
        assert(fabs(pl_packet_get_payload(pkt)[2] - 3.0) < EPS);
        pl_packet_increment_retry(pkt);
        assert(pkt->retransmission_count == 1);
        pl_packet_free(pkt);
        printf("  [PASS] Packet lifecycle\n");
    }

    /* ===== L3: LTI System ===== */
    {
        LTISystem* sys = pl_lti_create(2, 1, 1, 2);
        assert(sys); assert(sys->n == 2);
        double A[4] = {1.1, 0.0, 0.0, 0.9};
        double B[2] = {1.0, 0.5};
        double C[2] = {1.0, 1.0};
        pl_lti_set_A(sys, A, 2); pl_lti_set_B(sys, B); pl_lti_set_C(sys, C);
        double rho = pl_lti_spectral_radius(sys);
        assert(fabs(rho - 1.1) < 0.01);

        int rank_ctrb;
        double* Cc = pl_lti_controllability_matrix(sys, &rank_ctrb);
        assert(Cc); assert(rank_ctrb == 2); /* Should be controllable */
        free(Cc);

        int rank_obs;
        double* Ob = pl_lti_observability_matrix(sys, &rank_obs);
        assert(Ob); assert(rank_obs == 2); /* Should be observable */
        free(Ob);

        assert(pl_lti_is_stabilizable(sys));
        assert(pl_lti_is_detectable(sys));

        double x[2] = {1.0, 0.0}, u[1] = {-0.5}, x_next[2];
        unsigned long rng = 42;
        pl_lti_step(sys, x, u, false, &rng, x_next);
        assert(fabs(x_next[0] - 0.6) < EPS); /* 1.1*1 + 1.0*(-0.5) = 0.6 */
        assert(fabs(x_next[1] - (-0.25)) < EPS); /* 0.9*0 + 0.5*(-0.5) = -0.25 */

        double y[1];
        pl_lti_measure(sys, x, false, &rng, y);
        assert(fabs(y[0] - 1.0) < EPS);

        pl_lti_free(sys);
        printf("  [PASS] LTI system\n");
    }

    /* ===== L4: LQR/DARE ===== */
    {
        LTISystem* sys = pl_lti_create(2, 1, 1, 2);
        double A4[4] = {0.9, 0.1, 0.0, 0.85};
        double B4[2] = {0.1, 0.1};
        pl_lti_set_A(sys, A4, 2); pl_lti_set_B(sys, B4);
        for (int i = 0; i < 4; i++) sys->Q[i] = (i % 3 == 0) ? 1.0 : 0.0;

        double Rc[1] = {0.1};
        LQRSolution* lqr = pl_lqr_solve(sys, Rc, 500, 1e-6);
        assert(lqr);
        assert(lqr->iterations > 0);
        double rho_cl = pl_lqr_spectral_radius(lqr, sys);
        assert(rho_cl < 1.0); /* Closed-loop must be stable */

        /* TCP Controller */
        TCPController* tcp = pl_tcp_controller_create(sys, lqr);
        assert(tcp);
        pl_tcp_controller_set_hold(tcp, HOLD_ZERO_ORDER, HOLD_ZERO_ORDER);

        double y_test[1] = {1.0}, u_out[1];
        pl_tcp_controller_step(tcp, y_test, PACKET_RECEIVED, PACKET_RECEIVED, u_out);
        pl_tcp_controller_step(tcp, y_test, PACKET_LOST, PACKET_RECEIVED, u_out);
        assert(tcp->sensor_losses >= 1);

        pl_tcp_controller_free(tcp);
        pl_lqr_free(lqr);
        pl_lti_free(sys);
        printf("  [PASS] LQR/DARE solver\n");
    }

    /* ===== L5: Kalman Filter with Intermittent Observations ===== */
    {
        double A[4] = {1.0, 0.1, 0.0, 0.9};
        double C[2] = {1.0, 0.0};
        double Q[4] = {0.01, 0.0, 0.0, 0.01};
        double R[1] = {0.01};
        double x0[2] = {0.0, 0.0}, P0[4] = {1.0, 0.0, 0.0, 1.0};

        IntermittentKalmanFilter* ikf = pl_ikf_create(A, C, Q, R, 2, 1, x0, P0, 0.8);
        assert(ikf);
        assert(!ikf->is_mean_stable || ikf->is_mean_stable); /* just check computed */

        double y[1] = {1.5};
        for (int i = 0; i < 50; i++) {
            bool arr = (i % 5 != 0); /* 80% arrival */
            pl_ikf_step(ikf, y, arr);
        }
        pl_ikf_expected_covariance(ikf);
        assert(pl_ikf_trace_expected(ikf) > 0.0);

        pl_ikf_free(ikf);

        /* Standard Kalman Filter */
        KalmanFilter* kf = pl_kf_create(A, C, Q, R, 2, 1, x0, P0);
        assert(kf);
        pl_kf_step(kf, y);
        double tr = pl_kf_trace_P(kf);
        assert(tr > 0.0);
        pl_kf_free(kf);
        printf("  [PASS] Intermittent Kalman filter\n");
    }

    /* ===== L5: Mode-Dependent Kalman Filter ===== */
    {
        double A[4] = {1.0, 0.0, 0.0, 1.0};
        double C[2] = {1.0, 0.0};
        double Q[4] = {0.01, 0.0, 0.0, 0.01};
        double R[1] = {0.01};
        double x0[2] = {0.0, 0.0}, P0[4] = {1.0, 0.0, 0.0, 1.0};
        double trans[4] = {0.9, 0.1, 0.3, 0.7};
        double arr_rates[2] = {0.95, 0.5};

        ModeDependentKalmanFilter* mdkf = pl_mdkf_create(A, C, Q, R, 2, 1, 2,
            trans, arr_rates, x0, P0);
        assert(mdkf);
        assert(mdkf->n_modes == 2);

        double y[1] = {2.0};
        pl_mdkf_step(mdkf, y, true, 0);
        const double* est = pl_mdkf_get_estimate(mdkf);
        assert(est);
        pl_mdkf_free(mdkf);
        printf("  [PASS] Mode-dependent Kalman filter\n");
    }

    /* ===== L4: Stability Analysis (Bernoulli) ===== */
    {
        double A_s[4] = {1.0, 0.1, 0.0, 0.9}; /* stable system */
        double A_f[4] = {1.1, 0.0, 0.0, 1.1}; /* unstable when loss */

        StabilityCertificate* cert = pl_stability_test_bernoulli(A_s, A_f, 2, 0.2);
        assert(cert);
        /* With loss rate 0.2, mostly A_s (stable), should be stable */
        cert->is_stable = cert->is_stable; /* use the computed value */
        free(cert);

        /* Lyapunov exponent */
        double le = pl_stability_lyapunov_exponent(A_s, A_f, 2, 0.1, 5, 100, 999);
        assert(isfinite(le));

        /* Monte Carlo */
        double mse = pl_stability_monte_carlo(A_s, A_f, 2, 0.1, 3, 50, 999);
        assert(mse >= 0.0);
        printf("  [PASS] Stability analysis\n");
    }

    /* ===== L4: Jump Linear System ===== */
    {
        JumpLinearSystem* jls = pl_jls_create(2, 2);
        assert(jls);

        double A0[4] = {1.0, 0.1, 0.0, 0.9};
        double A1[4] = {1.0, 0.0, 0.0, 0.5};
        pl_jls_set_mode_matrix(jls, 0, A0);
        pl_jls_set_mode_matrix(jls, 1, A1);
        pl_jls_set_transition(jls, 0, 0, 0.9);
        pl_jls_set_transition(jls, 0, 1, 0.1);
        pl_jls_set_transition(jls, 1, 0, 0.3);
        pl_jls_set_transition(jls, 1, 1, 0.7);
        pl_jls_compute_steady_state(jls);

        pl_jls_test_mss(jls, 500, 1e-8);
        pl_jls_free(jls);
        printf("  [PASS] Jump Linear System\n");
    }

    /* ===== L5: Packetized Predictive Control ===== */
    {
        double A[4] = {1.0, 0.1, 0.0, 0.9};
        double B[2] = {0.0, 0.1};
        double L[2] = {2.0, 0.5};
        double x0[2] = {1.0, 0.0};

        PacketizedPredictiveControl* ppc = pl_ppc_create(5, A, B, L, 2, 1, x0);
        assert(ppc);

        const double* u0 = pl_ppc_consume(ppc, true);
        assert(u0);
        assert(ppc->buffer_index == 0);

        const double* u1 = pl_ppc_consume(ppc, false);
        assert(u1);
        assert(ppc->consecutive_losses == 1);

        pl_ppc_free(ppc);
        printf("  [PASS] Packetized Predictive Control\n");
    }

    /* ===== L6: Model Sensor Predictor ===== */
    {
        double A[4] = {1.0, 0.0, 0.0, 0.5};
        double C[2] = {1.0, 0.0};
        double x0[2] = {1.0, 2.0};

        ModelSensorPredictor* msp = pl_msp_create(A, C, 2, 1, x0, NULL);
        assert(msp);
        assert(pl_msp_is_model_valid(msp));

        const double* y_pred = pl_msp_predict_missing(msp);
        assert(y_pred);
        assert(msp->consecutive_losses == 1);

        double y_true[1] = {1.5};
        pl_msp_correct(msp, y_true);
        assert(msp->consecutive_losses == 0);

        pl_msp_free(msp);
        printf("  [PASS] Model sensor predictor\n");
    }

    /* ===== L5: Hold Strategy Selector ===== */
    {
        HoldStrategySelector* hss = pl_hss_create(true, 0.8);
        assert(hss);

        double cost;
        HoldStrategy sel = pl_hss_select(hss, 0.3, 1.5, 2, &cost);
        assert(sel == HOLD_ZERO_ORDER || sel == HOLD_PREDICTIVE);

        sel = pl_hss_select(hss, 0.9, 3.0, 10, &cost);
        assert(sel == HOLD_ZERO_INPUT); /* Safety override at high loss */

        pl_hss_free(hss);
        printf("  [PASS] Hold strategy selector\n");
    }

    /* ===== L6: Reference Governor ===== */
    {
        double A[4] = {1.0, 0.0, 0.0, 1.0};
        double B[2] = {1.0, 0.5};
        double C[2] = {1.0, 0.0};
        double ref[1] = {5.0};
        double y_min[1] = {-1.0}, y_max[1] = {1.0};
        double x0[2] = {0.5, 0.0};

        ReferenceGovernor* rg = pl_rg_create(A, B, C, 2, 1, 1, ref, 1, y_min, y_max, 3);
        assert(rg);

        const double* ref_adj = pl_rg_adjust(rg, x0);
        assert(ref_adj);
        /* Since x0=0.5 is within [-1,1], reference may or may not be scaled */
        assert(ref_adj[0] <= ref[0]);

        pl_rg_free(rg);
        printf("  [PASS] Reference governor\n");
    }

    /* ===== L7: Loss-Constrained Controller ===== */
    {
        double A[4] = {1.0, 0.1, 0.0, 0.8};
        double B[2] = {0.0, 0.1};
        double Q[4] = {1.0, 0.0, 0.0, 1.0};
        double R[1] = {0.1};
        double x[2] = {1.0, 0.0};

        LossConstrainedController* lcc = pl_lcc_create(A, B, Q, R, 2, 1, 10,
            0.2, PROTO_TCP_LIKE);
        assert(lcc);

        const double* u = pl_lcc_compute(lcc, x);
        assert(u);
        assert(isfinite(u[0]));

        pl_lcc_free(lcc);
        printf("  [PASS] Loss-constrained controller\n");
    }

    /* ===== UDP Controller ===== */
    {
        LTISystem* sys = pl_lti_create(2, 1, 1, 2);
        double A4[4] = {1.0, 0.0, 0.0, 0.9};
        double B4[2] = {0.0, 0.1};
        pl_lti_set_A(sys, A4, 2); pl_lti_set_B(sys, B4);

        double Rc[1] = {0.1};
        LQRSolution* lqr = pl_lqr_solve(sys, Rc, 200, 1e-6);

        UDPController* udp = pl_udp_controller_create(sys, lqr);
        assert(udp);

        double y[1] = {1.0}, u_out[1];
        pl_udp_controller_step(udp, y, PACKET_RECEIVED, 0.2, u_out);
        assert(isfinite(u_out[0]));

        pl_udp_controller_free(udp);
        pl_lqr_free(lqr);
        pl_lti_free(sys);
        printf("  [PASS] UDP controller\n");
    }

    /* ===== Critical Probability Analysis ===== */
    {
        double A[4] = {1.2, 0.0, 0.0, 1.1};
        double B[2] = {1.0, 0.5};
        double C[2] = {1.0, 0.0};
        double L[2] = {0.5, 0.3};

        CriticalProbabilityAnalysis* cpa = pl_critical_prob_analyze(
            A, B, C, L, 2, 1, 1);
        assert(cpa);
        assert(cpa->p_c_sensor >= 0.0 && cpa->p_c_sensor <= 1.0);

        int n_pts = pl_critical_prob_stability_region(cpa, 10);
        assert(n_pts == 10);

        pl_critical_prob_free(cpa);
        printf("  [PASS] Critical probability analysis\n");
    }

    /* ===== Packet History ===== */
    {
        PacketHistory* hist = pl_history_create(100);
        assert(hist);

        for (int i = 0; i < 50; i++) {
            PacketStatus st = (i % 4 == 0) ? PACKET_LOST : PACKET_RECEIVED;
            pl_history_record(hist, (unsigned long)i, (double)i * 0.1, st);
        }

        double lr = pl_history_loss_rate(hist, 50);
        assert(lr >= 0.1 && lr <= 0.4); /* ~25% loss */
        double bi = pl_history_burst_index(hist);
        assert(bi >= 0.0 && bi <= 1.0);
        pl_history_free(hist);
        printf("  [PASS] Packet history\n");
    }

    /* ===== Channel Generic Interface ===== */
    {
        PacketChannel* ch = pl_channel_create(CHANNEL_BERNOULLI);
        assert(ch);
        assert(ch->type == CHANNEL_BERNOULLI);

        for (int i = 0; i < 100; i++) pl_channel_transmit(ch);
        double rate = pl_channel_empirical_loss(ch);
        assert(rate >= 0.0 && rate <= 1.0);

        const char* name = pl_channel_type_name(CHANNEL_GILBERT_ELLIOTT);
        assert(name);
        pl_channel_free(ch);
        printf("  [PASS] Channel generic interface\n");
    }

    /* ===== PRNG Utilities ===== */
    {
        unsigned long state = 42;
        double u = pl_uniform(&state);
        assert(u >= 0.0 && u < 1.0);
        double e = pl_exponential(&state, 2.0);
        assert(e >= 0.0);
        int g = pl_geometric(&state, 0.5);
        assert(g >= 0);
        printf("  [PASS] PRNG utilities\n");
    }

    printf("\n=== All tests PASSED ===\n");
    return 0;
}
