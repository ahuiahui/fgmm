#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#ifndef _MSC_VER
  #include <sys/time.h> // gettimeofday .. 
#endif


#include "fgmm.h"

static char usage[] = "test_em infile #states \n infile text with space separated floats \n # state number of gaussians .. ";

int main(int argc,char ** argv)
{
  char * filename = argv[1];
  int n_states = atoi(argv[2]);
  char line[1024];
  char *pline;
  float *  data, *pdata;
  int n_data=0;
  int dim=0;

  int i=0,j=0,iterations=0;
  struct gmm * GMM;
  float lik;
  
  FILE * infile;
  FILE * sample_file;
  float * samp;
  if(argc < 3)
    {
      printf("%s",usage);
      exit(0);
    }

 
  srand(time(NULL));

  //  data = (float *) malloc(DIM*N_DATA*sizeof(float));


  infile=fopen(filename,"r");
  if(infile == NULL)
    {
      printf("can't open file\n");
      exit(0);
    }


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
  data = (float *) malloc(sizeof(float) * dim * n_data);
  pdata = data;

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

  printf("end loading file\n");
  fgmm_alloc(&GMM,n_states,dim);
  
#ifndef _MSC_VER
  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  fgmm_init_kmeans(GMM,data,n_data);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f s for kmenas\n",(t1.tv_sec*1000. + t1.tv_usec*.001));
#else  // no getimeofday .. 
  
  /* random initialization */ 
  fgmm_init_kmeans(GMM,data,n_data);

#endif

  fgmm_dump(GMM);

  // fgmm_em(GMM,data,n_data,&lik,1e-4,NULL);
  fgmm_dump(GMM);
  sample_file = fopen("samples.txt","w");
  samp = (float *) malloc(sizeof(float) * dim );
  for(i=0;i<100;i++)
    {
      fgmm_draw_sample(GMM,samp);
      for(j=0;j<dim;j++)
	fprintf(sample_file,"%f ",samp[j]);
      fprintf(sample_file,"\n");
    }   
  fgmm_free(&GMM);
  free(data);
  free(samp);
  return 0;
}

