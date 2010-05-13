#include "gmr.h"
#include "fgmm++.hpp"
#include <iostream>
#include <sys/time.h>

int main(int argc,char ** argv)
{
  if(argc < 3) 
    {
      std::cout << "usage ./bench nstate dim nsamples " << std::endl;
      exit(0);
    }
  int nstates = strtol(argv[1],NULL,10);
  int dim = strtol(argv[2],NULL,10);
  int nsamples = strtol(argv[3],NULL,10);
  Gmm mygmm(nstates,dim);
  for(int i=0;i<nstates;i++)
    {

      float mean[dim];
      float covar[dim*(dim+1)/2];
      int curs = 0;
      for(int j=0;j<dim;j++)
	{
	  if(j==(i%dim))
	    mean[j] = 1. + i/dim;
	  else 
	    mean[j] = 0.;
	  covar[curs] = 0.1;
	  curs++;
	  for(int k=j+1;k<dim;k++)
	    covar[curs++] = 0.;
	}
      mygmm.SetPrior(i,1./nstates);
      mygmm.SetMean(i,mean);
      mygmm.SetCovariance(i,covar);
    }
  float * data = new float[dim*nsamples];
  for(int i=0;i<nsamples;i++)
      mygmm.Draw(&data[i*dim]);
  
  struct timeval t1,t2;
   
  Gmm mygmm2(nstates,dim);
  mygmm2.init(data,nsamples);
  std::cout << "fgmm ::" << std::endl;
  gettimeofday(&t1,NULL);
  int niter = mygmm2.Em(data,nsamples);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/niter);
 
  Matrix mdata(data,nsamples,dim);

  GaussianMixture gmm;
  gmm.initEM_random(nstates,mdata);

  
  gettimeofday(&t1,NULL);
  gmm.doEM(mdata,0,10);
  gettimeofday(&t2,NULL);
  timersub(&t2,&t1,&t1);
  std::cout << "GMR ::" << std::endl;
  printf("%f ms / iterations\n",(t1.tv_sec*1000. + t1.tv_usec*.001)/10);
  //gmm.debug();
  return 0;
  
}
    
