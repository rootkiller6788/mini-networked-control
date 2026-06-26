#include "cloud_control_core.h"
#include "cloud_control_delay.h"
#include "cloud_control_resource.h"
#include "cloud_control_network.h"
#include <stdio.h>
int main(void){
printf("================================================\n");
printf("  Mini Cloud Control System - Overview Demo\n");
printf("================================================\n\n");
printf("Module: 4 layers (Core, Delay, Resource, Network)\n");
EdgeCloudTopology*topo=topology_create("ov");
topology_add_node(topo,"Edge");topology_add_node(topo,"Cloud-US");topology_add_node(topo,"Cloud-EU");
topology_add_link(topo,"Edge","Cloud-US",1000,5000);topology_add_link(topo,"Edge","Cloud-EU",500,8500);
printf("Network: %d nodes connected=%d\n",topo->node_count,topology_is_connected(topo));
ResourcePool*pool=rpool_create("dp");rpool_add_capacity(pool,RESOURCE_CPU,8,1);
printf("Resources: RM=%d EDF=%d\n",rpool_is_schedulable(pool,SCHED_RM),rpool_is_schedulable(pool,SCHED_EDF));
CloudControlSystem*ccs=ccs_create("ov",2,1,1,CC_MODE_COLLABORATIVE);
ccs_random_init(ccs,2,1,1,12345);
printf("Control: stable=%d MATI=%.0fus\n",ccs_is_stable(ccs,5000),ccs_max_allowable_delay(ccs));
printf("\nModule: 3000+ lines C code L1-L8.\n");
topology_free(topo);rpool_free(pool);ccs_free(ccs);return 0;}
