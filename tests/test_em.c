#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> // gettimeofday .. 
#include <time.h>

#include "em.h"
//#include "gmm.h"

#define N_STATES 6
#define DIM 3
#define N_DATA 10000

int main(int argc,char ** argv)
{
  srand(time(NULL));
  float *  data;
  data = (float *) malloc(DIM*N_DATA*sizeof(float));
  
  struct gmm GMM;
  gmm_alloc(&GMM,N_STATES,DIM);

  int state_i=0;
  int i=0,j=0;
  /* random initialization */ 
  for(state_i=0;state_i<N_STATES;state_i++)
    {
      GMM.gauss[state_i].prior = 1./N_STATES;
      for(j=0;j<DIM;j++)
	GMM.gauss[state_i].mean[j] = ((float)rand()/RAND_MAX)*2. - 1.;
      //dump(&GMM[state_i]);
    }
  FILE * infile;
  infile=fopen("test.txt","r");
  if(infile == NULL)
    printf("can't open file -- will gently segfault later :)) \n");
  float * pdata = data;
  for(i=0;i<N_DATA;i++)
    {
      for(j=0;j<DIM;j++)
	{
	  fscanf(infile,"%f",pdata);
	  pdata++;
	}
    }
  
  printf("end loading file\n");
  float lik;
  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  int iterations = em(&GMM,data,N_DATA,&lik);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/iterations);
  printf("%d iterations %f\n",iterations,lik);

  for(state_i=0;state_i<N_STATES;state_i++)
    {
      printf("Gaussian %d ::\n",state_i);
      dump(&GMM.gauss[state_i]);
    }
  /*  
  float samp[3];
  for(i=0;i<N_DATA;i++)
    {
      gmm_draw_sample(&GMM,samp);
      printf("%f %f %f\n",samp[0],samp[1],samp[2]);
    } 
  */
  gmm_free(&GMM);
  free(data);
  return 0;
}

