#include "gmm.h"
//#include "regression.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/**
 * do some regression on a sine dataset 
 */

int main(int argc, char ** argv)
{
  int data_len = 1000;
  int i=0;
  float * data = (float *)malloc( sizeof(float)*data_len*2);
  
  // generating data : 

  for(;i<data_len;i++)
    {
      data[2*i] = ((float)i)/data_len;
      data[2*i+1] = sin( 2 * data[2*i] * M_PI) + (0.1 * rand())/RAND_MAX;
      // printf("%f %f \n",data[2*i],data[2*i+1]); 
    }

  struct gmm gmm;
  gmm_alloc(&gmm,10,2);
  gmm_init_random(&gmm,data,data_len);
  float likelihood;
  em(&gmm,data,data_len,&likelihood);
  
  gmm_dump(&gmm);

  struct gmm_reg  * regression;
  
  gmm_regression_alloc_simple(&regression, &gmm, 1);
  gmm_regression_init(regression);
  /*
  float sample[2];
  FILE * sample_file = fopen("samples.txt","w");
  for(i=0;i<10000;i++)
    {
      gmm_draw_sample(&gmm,sample);
      fprintf(sample_file,"%f %f\n",sample[0],sample[1]);
    }
  fclose(sample_file);
  */


  FILE * reg_file = fopen("regression.txt","w");
  float input,output;
  for(i=0;i<100;i++)
    {
      input = ((float) i )/100;
      gmm_regression(regression,&input,&output);
      fprintf(reg_file,"%f %f %f\n",input,output,sin(2*input*M_PI));
      }
  fclose(reg_file);
  return 0;
}
				     
