#include "mex.h"
#include "fgmm.h"

/* The gateway function */
/*
 * (pi,mean,sigma) = gmm_em(data, nstates)  -> random shitty init
 * (pi,mu,sigma) = gmm_em(data, nstates ,  Priors, Mu  , Sigma)
 */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double * indata;

  double * in_priors;
  double * in_means;
  double * in_sigmas;

  short initialize = 0;
  if(!(nrhs == 2 || nrhs == 5) )
    {
      mexErrMsgIdAndTxt("EM","Two inputs required. data number of states");
    }

  if(nrhs == 5)
    initialize = 1;

  if(nlhs != 3) 
    {
      mexErrMsgIdAndTxt("EM","Three outputs required. [Priors, Means , Sigmas]");
    }

  
  indata = mxGetPr(prhs[0]);

  int nstates = (int) mxGetScalar(prhs[1]);
  mwSize dim = mxGetM(prhs[0]);
  mwSize len_data = mxGetN(prhs[0]);

  
  
  float * fdata = (float *) malloc( dim*len_data*sizeof(float));
  int i = 0;
  int j = 0;
  int k = 0;
  for(i=0;i<dim*len_data;i++)
    fdata[i] = (float) indata[i];

  struct gmm * GMM; 
  
  fgmm_alloc(&GMM,nstates,dim);
  fgmm_init_random(GMM,fdata,len_data);

  float * _mu;
  float * _sigma;
  if(initialize)
    {

      in_priors = mxGetPr(prhs[2]);
      in_means = mxGetPr(prhs[3]);
      in_sigmas = mxGetPr(prhs[4]);
      
      _mu = (float *) malloc(sizeof(float) * dim);
      _sigma = (float *) malloc(sizeof(float) * dim * dim);

      for(i=0;i<nstates;i++)
	{
	  for(j=0;j<dim;j++)
	    {
	      _mu[j] = (float) in_means[i*dim + j];
	      for(k=0;k<dim;k++)
		_sigma[j*dim + k] = (float) in_sigmas[i*dim*dim + j*dim + k];
	    }

	  fgmm_set_prior(GMM,i, (float) in_priors[i]);
	  fgmm_set_mean(GMM,i, _mu);
	  fgmm_set_covar(GMM,i,_sigma);
	  
	}
      free(_mu);
      free(_sigma);
    }
	  
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
	{
	  means[i + j*nstates] = (double) cmean[j];
	}
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
  
  
  

