#include "fgmm.h"
#include "gaussian.h"
#include <stdlib.h>
#include <math.h> // isinf , isnan 
#include <float.h>
#include <stdio.h>

#define max_iter 100


/**
 * for all data compute p(x|i) prob that state i generated data point x
 * correspond to the E step of EM. 
 *
 * @param *pix is an alloc'd float table of dimension nstates*data_len 
 * returns total log_likelihood 
 */
float fgmm_e_step(struct gmm * GMM,
		   const float * data,
		  int data_len,
		  float * pix)
{
      /* E step */
  float log_lik=0;
  float like;
  float * pxi;
  int data_i=0;
  int state_i;
  
  pxi = (float *) malloc(sizeof(float) * GMM->nstates);
  for(;data_i<data_len;data_i++)
    {
      like=0;
      for(state_i=0;state_i<GMM->nstates;state_i++)
	{
	  pxi[state_i] = gaussian_pdf(&GMM->gauss[state_i], data + data_i*GMM->dim) ;
	      
	  //printf("state %d -> lik : %f\n",state_i,pxi[state_i]);
	  like += pxi[state_i]* GMM->gauss[state_i].prior;
	  /* pdata++;
	     ppxi++; */
	}
      if(like<= FLT_MIN)
	{
	  printf("too far from current distrib %d\n",data_i);
	  exit(0);
	}
      log_lik += log(like);
      /* if(isnan(log_lik) || isinf(log_lik))
	 exit(0); */
      for(state_i=0;state_i<GMM->nstates;state_i++)
	{
	  pix[data_i + state_i*data_len] = pxi[state_i] * GMM->gauss[state_i].prior / like;
	  if(pix[data_i + state_i*data_len] <= FLT_MIN) 
	    pix[data_i + state_i*data_len] = FLT_MIN;
	      
	}
	  
    }
  free(pxi);
  return log_lik;
}

void fgmm_m_step(struct gmm * GMM,
		 const float * data,
		 int data_len,
		 float * pix)
{
  int state_i,k;
  for(state_i=0;state_i<GMM->nstates;state_i++)
    {
      GMM->gauss[state_i].prior = 0;
      for(k=0;k<GMM->dim;k++)
	GMM->gauss[state_i].mean[k] = 0;

      GMM->gauss[state_i].prior = smat_covariance(GMM->gauss[state_i].covar,
						  data_len,
						  &pix[state_i*data_len],
						  data,
						  GMM->gauss[state_i].mean);
      GMM->gauss[state_i].prior /= data_len;
      invert_covar(&GMM->gauss[state_i]);
    }
  

}

/** perform em on the giver data
 * @param data : the given dataset (data_length*3 floats) 
 *               /!\ aligned malloc'd float *
 * @param num_states : number of states of the GMM
 * @return  # of iterations 
 */
int fgmm_em( struct gmm * GMM,
	     const float * data,
	     int data_length, 
	     float * end_loglikelihood,
	     float likelihood_epsilon,
	     const float * weights) // if not NULL, weighted version .. 
{
  float * pix;
  float log_lik;
  int niter=0;
  float oldlik=0;
  float deltalik=0;
  int state_i;
  int d=0;
  
  pix = (float *) malloc( sizeof(float) * data_length * GMM->nstates);

  for(state_i=0;state_i<GMM->nstates;state_i++)
    {
      invert_covar(&GMM->gauss[state_i]);
    }


  for(niter=0;niter<max_iter;niter++)
    {
      
      log_lik = fgmm_e_step(GMM,data,data_length,pix);
      log_lik/=data_length;
      #ifndef NDEBUG 
      printf("Log lik :: %f \n",log_lik);
      #endif
      // M step 
      deltalik = log_lik - oldlik;
      oldlik = log_lik;
      
      if(fabs(deltalik) < likelihood_epsilon)
	break;
      
      if(weights != NULL) 
	{
	  for(d=0;d<data_length;d++)
	    {
	      for(state_i=0;state_i< GMM->nstates; state_i++)
		pix[d*GMM->nstates + state_i] *= weights[d];
		
	    }
	}

      fgmm_m_step(GMM,data,data_length,pix);
      //      pdata = data;
    }
  if(end_loglikelihood != NULL)
    *end_loglikelihood = log_lik;
  
  free(pix);
  return niter; 
}

