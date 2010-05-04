#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> // gettimeofday .. 
#include <time.h>
#include <string.h>

#include "gmm.h"

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
  while(1)
    {
      if(fgets(line,1024,infile)==NULL) // eof
	{
	  break;
	}
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
	  if(fscanf(infile,"%f",pdata)<1)
	    {
	      printf("error reading file");
	      exit(1);
	    }
	  pdata++;
	}
    }

  
  struct gmm GMM;
  gmm_alloc(&GMM,n_states,dim);
  /* random initialization */ 
  gmm_init_random(&GMM,data,n_data);
  
  printf("end loading file\n");
  float lik;
  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  int iterations = em(&GMM,data,n_data,&lik,1e-4);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/iterations);
  printf("%d iterations %f\n",iterations,lik);
  gmm_dump(&GMM);
  
  FILE * sample_file = fopen("samples.txt","w");

  float samp[dim];
  for(i=0;i<100;i++)
    {
      gmm_draw_sample(&GMM,samp);
      for(j=0;j<dim;j++)
	fprintf(sample_file,"%f ",samp[j]);
      fprintf(sample_file,"\n");
    } 
  
  gmm_free(&GMM);
  free(data);
  return 0;
}

