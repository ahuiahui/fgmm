#include <stdlib.h>
#include <stdio.h> 
#include <time.h>
#include <string.h>
#include "fgmm.h"

static char usage[] = "test_em infile #states \n infile text with space separated floats \n # state number of gaussians .. ";

int main(int argc,char ** argv)
{
  char * filename;
  int n_states; 
  FILE * infile;
  char line[1024];
  char *pline;

  int n_data=0;
  int dim=0;

  int i=0,j=0;
  _fgmm_real *data, *pdata,*samp;
  struct gmm * GMM;
  
  if(argc < 3)
    {
      printf("%s",usage);
      exit(0);
    }

  filename = argv[1];
  n_states = atoi(argv[2]);

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
  data = (_fgmm_real *) malloc(sizeof(_fgmm_real) * dim * n_data);
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

  fclose(infile);
  fgmm_alloc(&GMM,n_states,dim);
  /* random initialization */ 
  fgmm_init_random(GMM,data,n_data);
  
  printf("end loading file\n");
  
  for(i=0;i<n_data*10;i++) // 10 pass , should be eventually randomised 
    fgmm_update(GMM,&data[(i%n_data)*dim]);
  fgmm_dump(GMM);
  
  infile = fopen("samples.txt","w");

  samp = (_fgmm_real *) malloc(sizeof(_fgmm_real) *dim);
  for(i=0;i<1000;i++)
    {
      fgmm_draw_sample(GMM,samp);
      for(j=0;j<dim;j++)
	fprintf(infile,"%f ",samp[j]);
      fprintf(infile,"\n");
    } 
  fclose(infile);
  fgmm_free(&GMM);
  free(data);
  free(samp);
  return 0;
}

