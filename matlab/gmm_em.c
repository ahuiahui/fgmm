#include "mex.h"
#include "fgmm.h"

/* The gateway function */
/*
 * (pi,mean,sigma) = EM(data, nstates) 
 */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double * indata;
  if(nrhs != 2) 
    {
      mexErrMsgIdAndTxt("EM","Two inputs required. data and number of states");
    }
  
  if(nlhs != 3) 
    {
      mexErrMsgIdAndTxt("EM","Two outputs required. data and number of states");
    }

  
  indata = mxGetPr(prhs[0]);

  int nstates = (int) mxGetScalar(prhs[1]);
  mwSize dim = mxGetN(prhs[0]);
  mwSize len_data = mxGetM(prhs[0]);

  
  
  float * fdata = (float *) malloc( dim*len_data*sizeof(float));
  int i = 0;
  int j = 0;
  for(i=0;i<dim*len_data;i++)
    fdata[i] = (float) indata[i];

  struct gmm * GMM; 
  
  fgmm_alloc(&GMM,nstates,dim);
  fgmm_init_random(GMM,fdata,len_data);
  float like;
  fgmm_em(GMM,fdata,len_data,&like,1e-4);
  
  plhs[0] = mxCreateDoubleMatrix(1,nstates,mxREAL); 
  double * priors = mxGetPr(plhs[0]);
  
  for(i=0;i<nstates;i++)
    priors[i] = fgmm_get_prior(GMM,i);
      
  plhs[1] = mxCreateDoubleMatrix(nstates,dim ,mxREAL);
  
  double * means =mxGetPr(plhs[1]);
  float * cmean;
  for(i=0;i<nstates;i++)
    {
      cmean = fgmm_get_mean(GMM,i);
      for(j=0;j<dim;j++)
	  means[i*dim + j] = (double) cmean[j];
    }
  mwSize sig_dims[3];

  sig_dims[0] = dim;
  sig_dims[1] = dim;
  sig_dims[2] = nstates;

  plhs[2] = mxCreateNumericArray(3, sig_dims , mxDOUBLE_CLASS,mxREAL); 

  float fsigma[dim*dim];
  double * sigma = mxGetPr(plhs[2]);
  int size_sig = dim*dim;
  for(i=0;i<nstates;i++)
    {
      fgmm_get_covar(GMM,i,fsigma);
      for(j=0;j<dim*dim;j++)
	{
	  sigma[i*size_sig + j] = (double) fsigma[j];
	}
    }

  fgmm_free(&GMM);
  free(fdata);
}
  
  
  

