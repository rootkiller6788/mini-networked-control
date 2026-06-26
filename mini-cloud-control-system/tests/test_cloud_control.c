#include "cloud_control_core.h"
#include "cloud_control_delay.h"
#include "cloud_control_resource.h"
#include "cloud_control_network.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

static int passed = 0, failed = 0;
#define T(n) printf("  TEST: %s... ", n)
#define P() do { printf("PASS\n"); passed++; } while(0)
#define F(m) do { printf("FAIL: %s\n", m); failed++; } while(0)
#define C(c) do { if (!(c)) { F(#c); return; } } while(0)


static void t_create_free(void) {
    T("ccs_create/free");
    CloudControlSystem *ccs = ccs_create("t1", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    C(ccs != NULL); C(ccs->n==2 && ccs->m==1 && ccs->p==1);
    ccs_free(ccs); P();
}
static void t_create_null(void) {
    T("ccs_create NULL id");
    C(ccs_create(NULL, 2, 1, 1, CC_MODE_CLOUD_ONLY) == NULL); P();
}
static void t_plant_model(void) {
    T("ccs_set_plant_model");
    CloudControlSystem *ccs = ccs_create("pm", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    C(ccs != NULL);
    double A[4]={-2,0,1,-3}, B[2]={1,0.5}, C[2]={1,0}, D[1]={0};
    C(ccs_set_plant_model(ccs,A,B,C,D)==0);
    C(fabs(ccs->A[0][0]+2.0)<1e-9);
    ccs_free(ccs); P();
}
static void t_controller(void) {
    T("ccs_set_controller");
    CloudControlSystem *ccs = ccs_create("c", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    double K[2]={0.5,0.3}, L[2]={1,0.8};
    C(ccs_set_controller(ccs,K,L)==0);
    C(fabs(ccs->K[0][0]-0.5)<1e-9); C(fabs(ccs->L[0][0]-1.0)<1e-9);
    ccs_free(ccs); P();
}
static void t_compute_control(void) {
    T("ccs_compute_control");
    CloudControlSystem *ccs = ccs_create("ctl", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    double A[4]={-2,0,1,-3}, B[2]={1,0.5}, C[2]={1,0}, K[2]={1,0.5}, L[2]={0.5,0.3};
    ccs_set_plant_model(ccs,A,B,C,NULL); ccs_set_controller(ccs,K,L);
    ccs->x_hat[0]=1.0; ccs->x_hat[1]=0.5;
    double u[1]; C(ccs_compute_control(ccs,u)==0);
    C(fabs(u[0]-(-1.25))<1e-9); ccs_free(ccs); P();
}
static void t_apply_control(void) {
    T("ccs_apply_control");
    CloudControlSystem *ccs = ccs_create("ac", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    double A[4]={-1,0,0,-2}, B[2]={1,0}, C[2]={1,0};
    ccs_set_plant_model(ccs,A,B,C,NULL);
    ccs->plant_state.x[0]=2.0; ccs->plant_state.x[1]=3.0;
    double u[1]={1.0}; C(ccs_apply_control(ccs,u,0.1)==0);
    C(fabs(ccs->plant_state.x[0]-(2.0+0.1*(-1.0*2.0+1.0)))<1e-6);
    ccs_free(ccs); P();
}
static void t_step(void) {
    T("ccs_step");
    CloudControlSystem *ccs = ccs_create("st", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    double A[4]={-1,0,0,-2}, B[2]={1,0}, C[2]={1,0}, K[2]={0.5,0.2}, L[2]={0.3,0.1};
    ccs_set_plant_model(ccs,A,B,C,NULL); ccs_set_controller(ccs,K,L);
    double m[1]={1.0}; C(ccs_step(ccs,m,0.0,0.01,5000.0)==0);
    C(ccs->total_cycles==1); ccs_free(ccs); P();
}
static void t_stability(void) {
    T("ccs_is_stable");
    CloudControlSystem *ccs = ccs_create("sb", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    double A[4]={-3,0,0,-4}, B[2]={1,0.5}, C[2]={1,0}, K[2]={0,0}, L[2]={0,0};
    ccs_set_plant_model(ccs,A,B,C,NULL); ccs_set_controller(ccs,K,L);
    C(ccs_is_stable(ccs,1000.0)==1); ccs_free(ccs); P();
}
static void t_mati(void) {
    T("ccs_max_allowable_delay");
    CloudControlSystem *ccs = ccs_create("mt", 2, 1, 1, CC_MODE_CLOUD_ONLY);
    double A[4]={-3,0,0,-4}, B[2]={1,0.5}, C[2]={1,0}, K[2]={0,0}, L[2]={0,0};
    ccs_set_plant_model(ccs,A,B,C,NULL); ccs_set_controller(ccs,K,L);
    C(ccs_max_allowable_delay(ccs)>0); ccs_free(ccs); P();
}
static void t_delay_stats(void) {
    T("delay_stats_compute");
    double d[]={1000,2000,1500,3000,1200,1800,2500,900,1100,2100};
    DelayStatistics s; C(delay_stats_compute(d,10,&s)==0);
    C(s.sample_count==10); C(fabs(s.mean_us-1710.0)<200.0);
    C(s.min_us==900.0); C(s.max_us==3000.0); P();
}

static void t_delay_gen(void) {
    T("delay_generate");
    C(fabs(delay_generate(DELAY_MODEL_CONSTANT,5000,0)-5000.0)<1e-9);
    double du = delay_generate(DELAY_MODEL_UNIFORM,1000,3000);
    C(du>=1000.0 && du<=3000.0);
    C(delay_generate(DELAY_MODEL_EXPONENTIAL,5000,0)>=0); P();
}
static void t_delay_ema(void) {
    T("delay_estimate_online");
    double ema=1000.0; double nv=delay_estimate_online(&ema,2000.0,0.3);
    C(fabs(nv-(0.3*2000+0.7*1000))<1e-9); P();
}
static void t_delay_jitter(void) {
    T("delay_jitter_compute");
    double d[]={1000,1100,1050,1200,1150};
    C(delay_jitter_compute(d,5)>=0); P();
}
static void t_delay_stable(void) {
    T("delay_is_stable");
    CloudControlSystem *ccs = ccs_create("ds",2,1,1,CC_MODE_CLOUD_ONLY);
    double A[4]={-3,0,0,-4},B[2]={1,0.5},C[2]={1,0},K[2]={0,0},L[2]={0,0};
    ccs_set_plant_model(ccs,A,B,C,NULL); ccs_set_controller(ccs,K,L);
    C(delay_is_stable(ccs,1000.0)==1); ccs_free(ccs); P();
}
static void t_smith(void) {
    T("smith_predict");
    CloudControlSystem *ccs = ccs_create("sp",2,1,1,CC_MODE_CLOUD_ONLY);
    double A[4]={-1,0,0,-2},B[2]={1,0},C[2]={1,0},K[2]={0.1,0},L[2]={0.1,0};
    ccs_set_plant_model(ccs,A,B,C,NULL); ccs_set_controller(ccs,K,L);
    SmithPredictor *sp = smith_create(ccs); C(sp!=NULL);
    double ym[1]={1.5},u[1]={0.5},yc[1];
    C(smith_predict(sp,ym,u,0.01,yc)==0);
    smith_free(sp); ccs_free(ccs); P();
}
static void t_rpool(void) {
    T("rpool_allocate");
    ResourcePool *p = rpool_create("p1"); C(p!=NULL);
    C(rpool_add_capacity(p,RESOURCE_CPU,4.0,0.5)==0);
    C(rpool_add_capacity(p,RESOURCE_MEMORY,8192.0,1024.0)==0);
    ResourceType t[2]={RESOURCE_CPU,RESOURCE_MEMORY};
    double a[2]={1.0,2048.0};
    C(rpool_allocate(p,"task1",t,a,2,QOS_HIGH)==0);
    rpool_free(p); P();
}
static void t_sched(void) {
    T("RM/EDF bounds");
    double U1[]={0.2,0.3,0.1}; C(sched_rate_monotonic_bound(U1,3)==1);
    double U2[]={0.5,0.4,0.3}; C(sched_edf_bound(U2,3)==0);
    double U3[]={0.3,0.3,0.3}; C(sched_edf_bound(U3,3)==1); P();
}
static void t_resp_time(void) {
    T("sched_response_time");
    double C[]={1,1,2}, Tp[]={4,6,10}; int pri[]={0,1,2};
    C(sched_response_time(C,Tp,pri,3,0)>0); P();
}
static void t_elastic(void) {
    T("elastic_scale_up");
    ResourcePool *p = rpool_create("e1"); rpool_add_capacity(p,RESOURCE_CPU,8.0,0);
    ElasticController *ec = elastic_create(p,SCALE_REACTIVE); C(ec!=NULL);
    C(elastic_get_instance_count(ec)==1);
    C(elastic_scale_up(ec)==0); C(elastic_get_instance_count(ec)==2);
    elastic_free(ec); rpool_free(p); P();
}
static void t_workload(void) {
    T("wp_forecast");
    WorkloadPredictor *wp = wp_create(); C(wp!=NULL);
    for(int i=0;i<20;i++) wp_add_sample(wp,50.0+i*0.5,i*1.0);
    C(wp_forecast_holtwinters(wp,1.0)>0); C(wp_get_mape(wp)>=0);
    wp_free(wp); P();
}
static void t_topology(void) {
    T("topology shortest path");
    EdgeCloudTopology *topo = topology_create("tt"); C(topo!=NULL);
    topology_add_node(topo,"c"); topology_add_node(topo,"e1"); topology_add_node(topo,"e2");
    topology_add_link(topo,"e1","c",100,5000); topology_add_link(topo,"e2","c",50,8000);
    topology_add_link(topo,"c","e1",100,5000); topology_add_link(topo,"c","e2",50,8000);
    C(topology_is_connected(topo)==1);
    MultiPathRoute r; C(topology_find_shortest_path(topo,"e1","c",&r)==0);
    C(r.total_latency_us>0); topology_free(topo); P();
}
static void t_pktbuf(void) {
    T("pktbuf enqueue/dequeue");
    PacketBuffer *b = pktbuf_create("b1",100); C(b!=NULL);
    ControlPacket p; memset(&p,0,sizeof(p)); p.header.seq_num=1;
    C(pktbuf_enqueue(b,&p)==0); C(pktbuf_count(b)==1);
    ControlPacket p2; C(pktbuf_dequeue(b,&p2)==0);
    C(p2.header.seq_num==1); pktbuf_free(b); P();
}
static void t_bwalloc(void) {
    T("bw_alloc reserve");
    BandwidthAllocator *ba = bw_alloc_create("l1",1000.0); C(ba!=NULL);
    C(bw_alloc_reserve(ba,"c1",200.0,QOS_HIGH)==0);
    C(fabs(bw_alloc_get_available(ba)-800.0)<1e-9); bw_alloc_free(ba); P();
}
static void t_holdover(void) {
    T("edge holdover");
    EdgeNode e; memset(&e,0,sizeof(e)); e.max_holdover_time_us=1e6;
    edge_holdover_enter(&e); C(edge_holdover_check(&e)==1);
    edge_holdover_exit(&e); C(edge_holdover_check(&e)==0); P();
}
static void t_mode(void) {
    T("ccs_switch_mode");
    CloudControlSystem *ccs = ccs_create("ms",2,1,1,CC_MODE_CLOUD_ONLY);
    ccs_switch_mode(ccs,CC_MODE_EDGE_ONLY); C(ccs->mode==CC_MODE_EDGE_ONLY);
    ccs_free(ccs); P();
}
static void t_strings(void) {
    T("string conversions");
    C(strcmp(ccs_mode_string(CC_MODE_REDUNDANT),"REDUNDANT")==0);
    C(strcmp(ccs_state_string(CC_STATE_DEGRADED),"DEGRADED")==0);
    C(strcmp(net_condition_string(NET_CONDITION_CONGESTED),"CONGESTED")==0); P();
}
static void t_compare(void) {
    T("ccs_compare");
    CloudControlSystem *a=ccs_create("a",2,1,1,CC_MODE_CLOUD_ONLY);
    CloudControlSystem *b=ccs_create("b",2,1,1,CC_MODE_CLOUD_ONLY);
    double A[4]={-1,0,0,-2},B[2]={1,0},Cv[2]={1,0};
    ccs_set_plant_model(a,A,B,Cv,NULL); ccs_set_plant_model(b,A,B,Cv,NULL);
    C(ccs_compare(a,b)==1); b->A[0][0]=-999; C(ccs_compare(a,b)==0);
    ccs_free(a); ccs_free(b); P();
}
static void t_random(void) {
    T("ccs_random_init");
    CloudControlSystem *ccs = ccs_create("r",3,1,2,CC_MODE_CLOUD_ONLY);
    C(ccs_random_init(ccs,3,1,2,12345)==0);
    C(ccs->n==3&&ccs->m==1&&ccs->p==2); ccs_free(ccs); P();
}
static void t_perf(void) {
    T("ccs_compute_performance");
    CloudControlSystem *ccs = ccs_create("pf",2,1,1,CC_MODE_CLOUD_ONLY);
    double A[4]={-5,0,0,-6},B[2]={1,0},Cv[2]={1,0},K[2]={0.5,0},L[2]={0.2,0};
    ccs_set_plant_model(ccs,A,B,Cv,NULL); ccs_set_controller(ccs,K,L);
    ccs->reference.x[0]=1.0;
    for(int i=0;i<100;i++){double m[1]={ccs->plant_state.y[0]};ccs_step(ccs,m,i*0.01,0.01,1000);}
    C(ccs_compute_performance(ccs)==0); C(ccs->ise>0); ccs_free(ccs); P();
}
static void t_netcond(void) {
    T("network conditions");
    C(net_detect_condition(0.02,500,10000)==NET_CONDITION_NORMAL);
    C(net_detect_condition(0.6,1000,5000)==NET_CONDITION_DDOS);
    C(net_condition_delay_factor(NET_CONDITION_CONGESTED)>1.0); P();
}
static void t_cloud_lat(void) {
    T("cloud_latency_estimate");
    C(cloud_latency_estimate("us-east","us-east")<5000.0);
    C(cloud_latency_estimate("us-east","eu-west")>cloud_latency_estimate("us-east","us-east")); P();
}
static void t_pkt_ser(void) {
    T("pkt_serialize/deserialize");
    ControlPacket p1,p2; double y[2]={1.5,-0.5};
    pkt_create_measurement(&p1,"s1",y,2,42);
    uint8_t buf[512]; int len=pkt_serialize(&p1,buf,sizeof(buf));
    C(len>0); C(pkt_deserialize(buf,len,&p2)>0);
    C(p2.header.seq_num==42); C(fabs(p2.value[0]-1.5)<1e-9); P();
}
static void t_qos_map(void) {
    T("qos_resource_map");
    double c,me,bw; qos_resource_map(QOS_CRITICAL,&c,&me,&bw);
    C(c==2.0&&me==2048.0); qos_resource_map(QOS_BACKGROUND,&c,&me,&bw);
    C(c==0.1); P();
}

static void t_delay_fit(void) {
    T("delay_fit_exponential/normal");
    double d[]={1000,1200,900,1100,1300,1050,950,1150};
    double rate,mean,std;
    C(delay_fit_exponential(d,8,&rate)==0 && rate>0);
    C(delay_fit_normal(d,8,&mean,&std)==0 && fabs(mean-1081.25)<10.0); P();
}
static void t_delay_score(void) {
    T("delay_network_quality_score");
    double d[]={1000,1200,1100,900,1300}; DelayStatistics s;
    delay_stats_compute(d,5,&s);
    double sc=delay_network_quality_score(&s,2000.0);
    C(sc>0.0 && sc<=1.0); P();
}
static void t_delay_packet_corr(void) {
    T("delay_packet_loss_correlate");
    double d[]={1000,5000,1100,8000,1200};
    int lost[]={0,1,0,1,0};
    double corr=delay_packet_loss_correlate(d,lost,5);
    C(corr>=-1.0 && corr<=1.0); P();
}
static void t_network_link(void) {
    T("netlink_update_metrics");
    NetworkLink link; memset(&link,0,sizeof(link));
    link.bandwidth_mbps=100.0;
    netlink_update_metrics(&link,5000,500,0.02,0.6);
    C(link.latency_us>0); C(netlink_is_congested(&link)==0);
    C(netlink_effective_bandwidth(&link)<100.0); P();
}
static void t_mproute(void) {
    T("mproute functions");
    MultiPathRoute route; memset(&route,0,sizeof(route));
    NetworkLink link; memset(&link,0,sizeof(link));
    link.bandwidth_mbps=50.0; link.latency_us=3000; link.packet_loss_rate=0.01;
    C(mproute_add_link(&route,&link)==0);
    C(mproute_total_latency(&route)>0); C(mproute_effective_bw(&route)>0);
    MultiPathRoute r2; memset(&r2,0,sizeof(r2));
    NetworkLink l2; memset(&l2,0,sizeof(l2));
    l2.bandwidth_mbps=100.0; l2.latency_us=1000; l2.packet_loss_rate=0.001;
    mproute_add_link(&r2,&l2);
    C(mproute_compare(&r2,&route)==1); P();
}
static void t_rpool_ext(void) {
    T("rpool extended functions");
    ResourcePool *p = rpool_create("ext"); rpool_add_capacity(p,RESOURCE_CPU,8.0,2.0);
    ResourceType t[1]={RESOURCE_CPU}; double a[1]={2.0};
    rpool_allocate(p,"t1",t,a,1,QOS_HIGH);
    C(rpool_get_allocation_count(p)==1); C(rpool_get_overload_count(p)==0);
    rpool_rebalance(p); C(rpool_get_slack(p)>=0);
    rpool_release(p,"t1"); C(rpool_get_allocation_count(p)==0);
    rpool_free(p); P();
}
static void t_wp_ext(void) {
    T("wp extended functions");
    WorkloadPredictor *wp = wp_create();
    for(int i=0;i<30;i++) wp_add_sample(wp,50.0+(i%10)*10.0,i*1.0);
    double m=wp_mean_workload(wp); C(m>0);
    double v=wp_variance_workload(wp); C(v>=0);
    double pmr=wp_peak_to_mean_ratio(wp); C(pmr>=1.0);
    wp_set_parameters(wp,0.4,0.2,0.1,10);
    wp_free(wp); P();
}
static void t_elastic_ext(void) {
    T("elastic extended functions");
    ResourcePool *p = rpool_create("ee"); rpool_add_capacity(p,RESOURCE_CPU,16.0,0);
    ElasticController *ec = elastic_create(p,SCALE_REACTIVE);
    elastic_set_thresholds(ec,0.75,0.25);
    elastic_scale_up(ec); elastic_scale_up(ec);
    C(elastic_cost_estimate(ec,0.5)>0);
    elastic_free(ec); rpool_free(p); P();
}
static void t_delay_ext(void) {
    T("delay extended functions");
    double d[]={1000,2000,3000,1500,2500}; DelayStatistics s;
    delay_stats_compute(d,5,&s);
    double bnd=delay_predict_bound(&s,0.95); C(bnd>s.mean_us);
    double th=delay_adaptive_threshold(&s,2.0); C(th>s.mean_us);
    C(delay_is_anomalous(&s,50000.0,3.0)==1);
    double cost=delay_control_cost(NULL,d,5); C(cost<0); P();
}
int main(void) {
    printf("=== Cloud Control System Test Suite ===\n\n");
    t_create_free();
    t_create_null();
    t_plant_model();
    t_controller();
    t_compute_control();
    t_apply_control();
    t_step();
    t_stability();
    t_mati();
    t_delay_stats();
    t_delay_gen();
    t_delay_ema();
    t_delay_jitter();
    t_delay_stable();
    t_smith();
    t_rpool();
    t_sched();
    t_resp_time();
    t_elastic();
    t_workload();
    t_topology();
    t_pktbuf();
    t_bwalloc();
    t_holdover();
    t_mode();
    t_strings();
    t_compare();
    t_random();
    t_perf();
    t_netcond();
    t_cloud_lat();
    t_pkt_ser();
    t_qos_map();
    t_delay_fit();
    t_delay_score();
    t_delay_packet_corr();
    t_network_link();
    t_mproute();
    t_rpool_ext();
    t_wp_ext();
    t_elastic_ext();
    t_delay_ext();
    printf("\nResults: %d passed, %d failed\n", passed, failed);
    return failed > 0 ? 1 : 0;
}
