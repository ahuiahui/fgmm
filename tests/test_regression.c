#include "fgmm.h"
//#include "regression.h"
#define _USE_MATH_DEFINES  // force visual studio to define M_PI ... sigh .. 
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

/**
 * do some regression on a sine dataset 
 */

int main(int argc, char ** argv)
{
 
  int data_len = 1000;
  int i=0;
  float * data = (float *)malloc( sizeof(float)*data_len*2);
  struct gmm * gmm;
  float likelihood;
  struct fgmm_reg  * regression;
  float input,output,cvar;
  float error=0.;
  float x,y;
  float in[2];
  float out[2];
  float cov[3];
  float gt[2];
  float var = 0;
  int test_points = 1000;
  FILE * rfuke;
  
  srand ( time(NULL) );
  // generating data : 

  for(;i<data_len;i++)
    {
      data[2*i] = ((float)i)/data_len;
      data[2*i+1] = sin( 2 * data[2*i] * M_PI) + (0.1 * rand())/RAND_MAX;
      // printf("%f %f \n",data[2*i],data[2*i+1]); 
    }

  fgmm_alloc(&gmm,10,2);
  fgmm_init_random(gmm,data,data_len);

  fgmm_em(gmm,data,data_len,&likelihood,1e-3);
  
  // fgmm_dump(&gmm);

  
  fgmm_regression_alloc_simple(&regression, gmm, 1);
  fgmm_regression_init(regression);

  for(i=0;i<100;i++)
    {
      input = ((float) i )/100;
      fgmm_regression(regression,&input,&output,&cvar);
      error += pow(output-sin(2*input*M_PI),2);
      }
  error /= 100.;
  error = sqrtf(error);
  printf("RMS error on the toy - sine dataset :: %f\n",error);
  assert(error < .25);
  printf(".. pass \n");
  
  fgmm_free(&gmm);
  fgmm_regression_free(&regression);
  exit(0);
  free(data);
 
  /* multidim 4 d data set */ 
  // data_len = 100000;
  data = malloc( sizeof(float) * data_len * 4);
 
  for(i=0;i<data_len;i++)
    {
      x = (M_PI*2 / 500) *  (rand()%500) ;
      y =  (M_PI*2/500)  * (rand()%500 );
      data[4*i] = x ;
      data[4*i+1] = y ;
      data[4*i+2] = sin(x) + cos(y) +  (0.1 * rand())/RAND_MAX; 
      data[4*i+3] = cos(2*x) -  sin(y + 2) +  (0.1 * rand())/RAND_MAX; 
      //printf("%f %f %f %f \n",data[4*i], data[4*i+1], data[4*i+2], data[4*i+3]);
    }
  
  fgmm_alloc(&gmm , 16 , 4 );
  fgmm_init_random(gmm,data,data_len);
  fgmm_em(gmm,data,data_len,&likelihood,1e-2);
 
  fgmm_regression_alloc_simple(&regression, gmm, 2);
  fgmm_regression_init(regression);

  error = 0.;
  //rfuke = fopen("result.txt","w");
  for(i=0;i<test_points;i++)
    {
      x = (M_PI*2 / 500) *  (rand()%500) ;
      y =  (M_PI*2/500)  * (rand()%500 );

      in[0] = x;
      in[1] = y;

      fgmm_regression(regression,in,out,cov);
      
      gt[0] =   sin(x) + cos(y);
      gt[1] =   cos(2*x) - sin(y+2);
      //fprintf(rfuke,"%f %f %f %f %f %f \n",x,y,gt[0],gt[1],out[0],out[1]);

      error += pow(out[0] - gt[0] ,2);
      error += pow(out[1] - gt[1] ,2);
      var += gt[0]*gt[0] + gt[1]*gt[1];
      }

  error /= test_points;
  var /= test_points ;
  error = sqrtf(error);
  var = sqrtf(var);
  //fclose(rfuke);
  printf("nRMS error on the 4d sine dataset :: %f\n",error);
  assert(error < .75);
  printf(".. pass \n");
  return 0;
}
				     
