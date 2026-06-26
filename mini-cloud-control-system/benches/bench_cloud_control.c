#include "cloud_control_core.h"
#include "cloud_control_delay.h"
#include "cloud_control_resource.h"
#include "cloud_control_network.h"
#include <stdio.h>
#include <time.h>
#define ITER 10000
int main(void) {
    printf("=== Cloud Control Benchmarks ===
");
    CloudControlSystem *ccs = ccs_create("bench", 4, 2, 2, CC_MODE_CLOUD_ONLY);
    ccs_random_init(ccs, 4, 2, 2, 42);
    /* Benchmark control loop */
    clock_t start = clock();
    for (int i = 0; i < ITER; i++) {
        double m[2] = {1.0, 2.0}, u[2];
        ccs_compute_control(ccs, u);
        ccs_apply_control(ccs, u, 0.001);
    }
    clock_t end = clock();
    double dt = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Control loop: %d iters in %.3fs = %.0f Hz
", ITER, dt, ITER/dt);
    /* Benchmark delay statistics */
    double delays[1000];
    delay_trace_generate(DELAY_MODEL_EXPONENTIAL, 5000, 0, 1000, delays);
    start = clock();
    for (int i = 0; i < 1000; i++) {
        DelayStatistics s; delay_stats_compute(delays, 1000, &s);
    }
    end = clock();
    dt = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Delay stats: 1000x1000 samples in %.3fs
", dt);
    /* Benchmark topology */
    EdgeCloudTopology *topo = topology_create("bench");
    for (int i = 0; i < 100; i++) {
        char nid[16]; sprintf(nid, "N%d", i); topology_add_node(topo, nid);
    }
    for (int i = 0; i < 99; i++) {
        char src[16], dst[16]; sprintf(src,"N%d",i); sprintf(dst,"N%d",i+1);
        topology_add_link(topo, src, dst, 100, 5000);
    }
    start = clock();
    for (int i = 0; i < 1000; i++) {
        MultiPathRoute r; topology_find_shortest_path(topo, "N0", "N99", &r);
    }
    end = clock();
    dt = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Shortest path (100 nodes): 1000 runs in %.3fs
", dt);
    topology_free(topo); ccs_free(ccs);
    return 0;
}
