#include "packet_loss_core.h"
#include "packet_loss_controller.h"
#include "packet_loss_estimator.h"
#include <stdio.h>
#include <time.h>

#define BENCH_ITERS 100000
#define BENCH_WARMUP 1000

double time_now_ms(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC * 1000.0;
}

int main(void) {
    printf("=== Packet Loss Control Benchmarks ===\n");

    /* Benchmark: Bernoulli channel throughput */
    {
        BernoulliChannel* ch = pl_bernoulli_create(0.2, 42);
        for (int i = 0; i < BENCH_WARMUP; i++) pl_bernoulli_transmit(ch);
        double start = time_now_ms();
        for (int i = 0; i < BENCH_ITERS; i++) pl_bernoulli_transmit(ch);
        double elapsed = time_now_ms() - start;
        printf("Bernoulli transmit: %.2f Mpkts/s  (%d in %.1f ms)\n",
               BENCH_ITERS / elapsed / 1000.0, BENCH_ITERS, elapsed);
        pl_bernoulli_free(ch);
    }

    /* Benchmark: Gilbert-Elliott channel throughput */
    {
        GilbertElliottChannel* ge = pl_gilbert_elliott_create(0.1, 0.3, 0.01, 0.3, 42);
        for (int i = 0; i < BENCH_WARMUP; i++) pl_gilbert_elliott_transmit(ge);
        double start = time_now_ms();
        for (int i = 0; i < BENCH_ITERS; i++) pl_gilbert_elliott_transmit(ge);
        double elapsed = time_now_ms() - start;
        printf("Gilbert-Elliott transmit: %.2f Mpkts/s\n",
               BENCH_ITERS / elapsed / 1000.0);
        pl_gilbert_elliott_free(ge);
    }

    /* Benchmark: Kalman filter step */
    {
        double A[4] = {1.0, 0.1, 0.0, 0.9};
        double C[2] = {1.0, 0.0};
        double Q[4] = {0.01, 0.0, 0.0, 0.01};
        double R[1] = {0.01};
        double x0[2] = {0.0, 0.0}, P0[4] = {1.0, 0.0, 0.0, 1.0};
        KalmanFilter* kf = pl_kf_create(A, C, Q, R, 2, 1, x0, P0);

        double y[1] = {2.0};
        for (int i = 0; i < BENCH_WARMUP; i++) pl_kf_step(kf, y);
        double start = time_now_ms();
        for (int i = 0; i < BENCH_ITERS; i++) pl_kf_step(kf, y);
        double elapsed = time_now_ms() - start;
        printf("Kalman filter step: %.2f ksteps/s\n",
               BENCH_ITERS / elapsed);
        pl_kf_free(kf);
    }

    /* Benchmark: LTI step */
    {
        LTISystem* sys = pl_lti_create(4, 2, 2, 4);
        double A[16] = {0}; for (int i = 0; i < 4; i++) A[i*5] = 0.9;
        pl_lti_set_A(sys, A, 4);
        double x[4] = {1.0, 0.0, 0.0, 0.0}, u[2] = {0.5, -0.3};
        unsigned long rng = 42;

        double x_next[4];
        for (int i = 0; i < BENCH_WARMUP; i++) pl_lti_step(sys, x, u, false, &rng, x_next);
        double start = time_now_ms();
        for (int i = 0; i < BENCH_ITERS; i++) pl_lti_step(sys, x, u, false, &rng, x_next);
        double elapsed = time_now_ms() - start;
        printf("LTI step (4-state): %.2f ksteps/s\n", BENCH_ITERS / elapsed);
        pl_lti_free(sys);
    }

    /* Benchmark: LQR solve (small system) */
    {
        LTISystem* sys = pl_lti_create(2, 1, 1, 2);
        double A[4] = {1.0, 0.1, 0.0, 0.9}, B[2] = {0.0, 0.1};
        pl_lti_set_A(sys, A, 2); pl_lti_set_B(sys, B);
        double Rc[1] = {0.1};
        LQRSolution* sol = pl_lqr_solve(sys, Rc, 100, 1e-6);
        printf("LQR solve (2-state): %d iterations\n", sol->iterations);
        pl_lqr_free(sol);
        pl_lti_free(sys);
    }

    printf("=== Benchmarks complete ===\n");
    return 0;
}
