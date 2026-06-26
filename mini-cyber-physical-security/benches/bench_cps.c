#include "cps_security_core.h"
#include "cps_detection.h"
#include <stdio.h>
#include <time.h>

int main(void) {
    clock_t start, end;
    int n_iter = 10000;

    printf("=== CPS Security Benchmarks ===\n\n");

    /* Matrix multiply */
    double A[100] = {0}; for(int i=0;i<10;i++) A[i*10+i]=1.0;
    double B[100] = {0}; for(int i=0;i<10;i++) B[i*10+i]=1.0;
    double C[100];
    start = clock();
    for (int i = 0; i < n_iter; i++)
        cps_matrix_multiply(C, A, B, 10, 10, 10);
    end = clock();
    printf("Matrix 10x10 multiply x%d: %.3f ms\n", n_iter,
           1000.0 * (end - start) / CLOCKS_PER_SEC);

    /* Chi-squared test */
    CPSDetector det;
    cps_detector_init(&det, CPS_DETECT_CHI2, 3.841);
    double res[2] = {0.5, -0.3};
    double cov[2] = {1.0, 1.0};
    start = clock();
    for (int i = 0; i < n_iter * 10; i++)
        cps_chi2_test(&det, res, cov, 2);
    end = clock();
    printf("Chi2 test x%d: %.3f ms\n", n_iter * 10,
           1000.0 * (end - start) / CLOCKS_PER_SEC);

    /* System step */
    CPSSecuritySystem* sys = cps_security_create(2, 1, 1);
    double A2[4] = {1,0.1,0,1}, B2[2]={0,0.1}, C2[2]={1,0}, x0[2]={1,0};
    cps_set_state_matrix(sys, A2, 2, 2);
    cps_set_input_matrix(sys, B2, 2, 1);
    cps_set_output_matrix(sys, C2, 1, 2);
    cps_set_initial_state(sys, x0);
    double u[1] = {0.0};
    start = clock();
    for (int i = 0; i < n_iter; i++)
        cps_step(sys, u);
    end = clock();
    printf("System step x%d: %.3f ms\n", n_iter,
           1000.0 * (end - start) / CLOCKS_PER_SEC);
    cps_security_free(sys);
    cps_detector_free_internals(&det);

    printf("\n=== Benchmarks complete ===\n");
    return 0;
}
