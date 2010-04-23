#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> // gettimeofday .. 
#include <time.h>
#include <string.h>

#include "em.h"

static char usage[] = "test_em infile #states \n infile text with space separated floats \n # state number of gaussians .. ";

int main(int argc,char ** argv)
{
  if(argc < 3)
    {
      printf("%s",usage);
      exit(0);
    }

  char * filename = argv[1];
  int n_states = atoi(argv[2]);

  srand(time(NULL));

  //  data = (float *) malloc(DIM*N_DATA*sizeof(float));

  FILE * infile;

  infile=fopen(filename,"r");
  if(infile == NULL)
    {
      printf("can't open file\n");
      exit(0);
    }


  char line[1024];
  char *pline;

  int n_data=0;
  int dim=0;

  int i=0;
  while(!feof(infile))
    {
      fgets(line,1024,infile);
      if(dim==0)
	{
	  pline = strtok(line," ");
	  while(pline != NULL)
	    {
	      dim++;
	      pline=strtok(NULL," ");
	    }
	}
      n_data++;
    }
  n_data--;
  dim--;
  rewind(infile);
  printf("%d , %d\n",n_data,dim);
  float *  data;
  data = (float *) malloc(sizeof(float) * dim * n_data);
  float * pdata = data;
  
  int j;
  for(i=0;i<n_data;i++)
    {
      for(j=0;j<dim;j++)
	{
	  fscanf(infile,"%f",pdata);
	  pdata++;
	}
    }

  
  struct gmm GMM;
  gmm_alloc(&GMM,n_states,dim);

  int state_i=0;
  //int i=0,j=0;
  /* random initialization */ 
  for(state_i=0;state_i<n_states;state_i++)
    {
      GMM.gauss[state_i].prior = 1./n_states;
      for(j=0;j<dim;j++)
	GMM.gauss[state_i].mean[j] = ((float)rand()/RAND_MAX)*2. - 1.;
      //dump(&GMM[state_i]);
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

