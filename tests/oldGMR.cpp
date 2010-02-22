#include "gmr.h"
#include <sys/time.h>

#define N_STATES 6
#define DIM 3
#define N_DATA 10000


int main(int arc,char ** argv)
{
  Matrix data(N_DATA,DIM);
  int i,j;
  FILE * infile;
  infile=fopen("test.txt","r");
  if(infile == NULL)
    printf("can't open file -- will gently segfault later :)) \n");
  float pdata;
  for(i=0;i<N_DATA;i++)
    {
      for(j=0;j<DIM;j++)
	{
	  fscanf(infile,"%f",&pdata);
	  data(i,j) = pdata;
	}
    }

  GaussianMixture gmm;
  gmm.initEM_random(N_STATES,data);

  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  gmm.doEM(data,0,10);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/10);
  gmm.debug();
  return 0;
}
