#include "fgmm.h"
#include "gaussian.h"
//#include "vec3.h"
#include <stdlib.h>
#include <math.h> // isinf , isnan 
#include <float.h>
#include <stdio.h>

#define max_iter 100


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
	     float likelihood_epsilon)
{
  int data_i=0;
  int state_i=0;
  int k=0;
  float * pxi;
  float * pix;
  float like;
  float log_lik;
  int niter=0;

  pxi = (float *) malloc(sizeof(float)* GMM->nstates);
  pix = (float *) malloc( sizeof(float) * data_length * GMM->nstates);

  float oldlik=0;
  float deltalik=0;
  
  for(state_i=0;state_i<GMM->nstates;state_i++)
    {
      invert_covar(&GMM->gauss[state_i]);
    }


  for(niter=0;niter<max_iter;niter++)
    {
      /* E step */
      log_lik=0;
      for(data_i=0;data_i<data_length;data_i++)
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
	      pix[data_i + state_i*data_length] = pxi[state_i] * GMM->gauss[state_i].prior / like;
	    }
	  
	}
      
      log_lik/=data_length;
      printf("Log lik :: %f \n",log_lik);
      // M step 
      deltalik = log_lik - oldlik;
      oldlik = log_lik;
      
      if(fabs(deltalik) < likelihood_epsilon)
	break;
      
      //      pdata = data;
      for(state_i=0;state_i<GMM->nstates;state_i++)
	{
	  GMM->gauss[state_i].prior = 0;
	  for(k=0;k<GMM->dim;k++)
	    GMM->gauss[state_i].mean[k] = 0;

	  GMM->gauss[state_i].prior = smat_covariance(GMM->gauss[state_i].covar,
					       data_length,
					       &pix[state_i*data_length],
					       data,
					       GMM->gauss[state_i].mean);
	  GMM->gauss[state_i].prior /= data_length;
	  invert_covar(&GMM->gauss[state_i]);
	  /*printf("gauss : %d :: \n",state_i);
	    dump(&GMM[state_i]); 
	    printf("%f\n",GMM[state_i].nfactor); */	  
	}
    }
  if(end_loglikelihood != NULL)
    *end_loglikelihood = log_lik;
  
  free(pxi);
  free(pix);
  return niter; 
}

