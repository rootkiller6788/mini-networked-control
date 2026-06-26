#include "cloud_control_core.h"
#include "cloud_control_delay.h"
#include <stdio.h>
#include <math.h>
int main(void){
printf("=== Delay Compensation Demo ===\n");
printf("Smith predictor vs uncompensated.\n\n");
CloudControlSystem*ccs=ccs_create("dd",2,1,1,CC_MODE_CLOUD_ONLY);
double A[4]={0.5,1,0,0.5},B[2]={0,1},C[2]={1,0},K[2]={2.5,2.0},L[2]={1.5,1.0};
ccs_set_plant_model(ccs,A,B,C,NULL);ccs_set_controller(ccs,K,L);
double ref[2]={10,0};ccs_set_reference(ccs,ref);ccs->plant_state.x[0]=0;
SmithPredictor*sp=smith_create(ccs);smith_set_delay(sp,10000);
printf("Delay:10ms Plant:unstable open-loop\n\n");
printf("%-8s %-14s %-14s %-14s\n","Step","NoComp","Smith","Control");
CloudControlSystem*ccs2=ccs_create("dd2",2,1,1,CC_MODE_CLOUD_ONLY);
ccs_set_plant_model(ccs2,A,B,C,NULL);ccs_set_controller(ccs2,K,L);ccs_set_reference(ccs2,ref);
for(int i=0;i<40;i++){double m1[1]={ccs->plant_state.x[0]},m2[1]={ccs2->plant_state.x[0]};delay_compensated_step(ccs,sp,m1,i*0.02,0.02,COMPENSATE_SMITH);delay_compensated_step(ccs2,NULL,m2,i*0.02,0.02,COMPENSATE_NONE);double u[1];ccs_compute_control(ccs,u);printf("%-8d %-14.4f %-14.4f %-14.4f\n",i,ccs2->plant_state.x[0],ccs->plant_state.x[0],u[0]);}
ccs_compute_performance(ccs);ccs_compute_performance(ccs2);
printf("\nWith Smith: ISE=%.4f settling=%.3fs\n",ccs->ise,ccs->settling_time);
printf("Without Smith: ISE=%.4f settling=%.3fs\n",ccs2->ise,ccs2->settling_time);
smith_free(sp);ccs_free(ccs);ccs_free(ccs2);return 0;}
