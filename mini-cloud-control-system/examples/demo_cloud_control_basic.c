#include "cloud_control_core.h"
#include "cloud_control_delay.h"
#include <stdio.h>
#include <math.h>
int main(void){
printf("=== Cloud Control Basic Demo ===\n");
printf("Simulating 2nd-order plant via cloud with delay.\n\n");
CloudControlSystem*ccs=ccs_create("demo",2,1,1,CC_MODE_CLOUD_ONLY);
double A[4]={0,1,0,0},B[2]={0,1},C[2]={1,0},K[2]={1.0,1.5},L[2]={2.0,1.0};
ccs_set_plant_model(ccs,A,B,C,NULL);ccs_set_controller(ccs,K,L);
double ref[2]={5.0,0.0};ccs_set_reference(ccs,ref);
ccs->plant_state.x[0]=0;ccs->plant_state.x[1]=0;
printf("Plant: double integrator, target=5.0\nCloud delay=5ms RTT\n\n");
printf("%-8s %-12s %-12s %-12s\n","Step","Position","Velocity","Control");
for(int i=0;i<50;i++){double m[1]={ccs->plant_state.x[0]};ccs_step(ccs,m,i*0.02,0.02,5000);double u[1];ccs_compute_control(ccs,u);printf("%-8d %-12.4f %-12.4f %-12.4f\n",i,ccs->plant_state.x[0],ccs->plant_state.x[1],u[0]);}
ccs_compute_performance(ccs);
printf("\nPerformance: settling=%.3fs overshoot=%.1f%% ISE=%.4f\n",ccs->settling_time,ccs->overshoot_percent,ccs->ise);
ccs_free(ccs);return 0;}
