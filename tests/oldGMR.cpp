#include "gmr.h"
#include <sys/time.h>


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
  REALTYPE *  data;
  data = (REALTYPE *) malloc(sizeof(REALTYPE) * dim * n_data);
  REALTYPE * pdata = data;
  
  int j;
  for(i=0;i<n_data;i++)
    {
      for(j=0;j<dim;j++)
	{
	  fscanf(infile,"%f",pdata);
	  pdata++;
	}
    }

  Matrix mdata(data,n_data,dim);

  GaussianMixture gmm;
  gmm.initEM_random(n_states,mdata);

  struct timeval t1,t2;
  gettimeofday(&t1,NULL);
  gmm.doEM(mdata,0,10);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/10);
  //gmm.debug();
  return 0;
}
