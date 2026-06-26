#include "cloud_control_core.h"
#include "cloud_control_resource.h"
#include <stdio.h>
#include <math.h>
int main(void){
printf("=== Cloud Resource Scaling Demo ===\n");
printf("Elastic scaling and workload prediction.\n\n");
ResourcePool*pool=rpool_create("cloud-pool");
rpool_add_capacity(pool,RESOURCE_CPU,16,2);rpool_add_capacity(pool,RESOURCE_MEMORY,32768,4096);rpool_add_capacity(pool,RESOURCE_NETWORK,10000,1000);
printf("Pool: 16CPU 32GB 10Gbps\n\n");
ElasticController*ec=elastic_create(pool,SCALE_REACTIVE);elastic_set_thresholds(ec,0.80,0.30);
WorkloadPredictor*wp=wp_create();wp_set_parameters(wp,0.3,0.1,0.05,12);
printf("Simulating varying workload:\n");
double bl[]={0.3,0.5,0.8,0.9,0.7,0.4,0.3,0.6,0.85,0.95};int n=10;
for(int i=0;i<n;i++){ResourceType t[2]={RESOURCE_CPU,RESOURCE_MEMORY};double a[2]={bl[i]*4,bl[i]*4096};rpool_allocate(pool,"task",t,a,2,QOS_MEDIUM);rpool_rebalance(pool);wp_add_sample(wp,pool->total_utilization,i*60);elastic_evaluate(ec,i*60);printf("  t=%dmin: util=%.2f instances=%d forecast=%.2f\n",i,pool->total_utilization,elastic_get_instance_count(ec),wp_forecast_holtwinters(wp,1));rpool_release(pool,"task");}
printf("\nScale-up: %llu Scale-down: %llu\n",(unsigned long long)ec->scale_up_count,(unsigned long long)ec->scale_down_count);
double utils[]={0.25,0.30,0.15,0.10};
printf("\nScheduling: 4 tasks RM=%s EDF=%s\n",sched_rate_monotonic_bound(utils,4)?"PASS":"FAIL",sched_edf_bound(utils,4)?"PASS":"FAIL");
double c,m,b;printf("\nQoS mapping:\n");
qos_resource_map(QOS_CRITICAL,&c,&m,&b);printf("  CRITICAL: %.1f CPU %.0f MB %.0f Mbps\n",c,m,b);
qos_resource_map(QOS_HIGH,&c,&m,&b);printf("  HIGH: %.1f CPU %.0f MB %.0f Mbps\n",c,m,b);
elastic_free(ec);wp_free(wp);rpool_free(pool);return 0;}
