#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> // gettimeofday .. 
#include "em.h"

int main(int argc,char ** argv)
{
  int num_state = 5;
  float *  data;
  data = (float *) malloc(300000*sizeof(float));
  //   data = (float ** ) malloc(3000*3*sizeof(float)]; 
  int state_i;
  struct gaussian * GMM;
  GMM = (struct gaussian *) malloc(sizeof(struct gaussian) * num_state);
    
  /* random initialization */ 
  for(state_i=0;state_i<num_state;state_i++)
    {
      init_gaussian(&GMM[state_i],3);
      GMM[state_i].prior = 1./num_state;
      // dump(&GMM[state_i]);
    }


  int i=0;
  FILE * infile;
  infile=fopen("test.txt","r");
  if(infile == NULL)
    printf("can't open file -- will gently segfault later :)) \n");
  float * pdata = data;
  for(i=0;i<100000;i++)
    {
      fscanf(infile,"%f  %f  %f\n",pdata,pdata +1,pdata +2);
      pdata +=3 ;
    }
  
  printf("end loading file\n");
  float lik;
  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  int iterations = em(GMM,data,3,100000,num_state,&lik);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/iterations);
  printf("%d iterations %f\n",iterations,lik);

  for(state_i=0;state_i<num_state;state_i++)
    {
      printf("Gaussian %d ::\n",state_i);
      dump(&GMM[state_i]);
    }
  
  return 0;
}

