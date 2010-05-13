#include "fgmm.h"
#include "gaussian.h"
#include <stdio.h>
#include <stdlib.h>

void fgmm_alloc(struct gmm * gm,int nstates,int dim)
{
  int i=0;
  gm->nstates = nstates;
  gm->dim = dim;
  gm->gauss = (struct gaussian *) malloc(sizeof(struct gaussian) * nstates );
  
  for(i=0;i<nstates;i++)
    gaussian_init(&gm->gauss[i],dim); // this alloc memory for each gaussian 
  /*
      GMM[state_i].prior = 1./N_STATES;
      for(j=0;j<DIM;j++)
	GMM[state_i].mean[j] = ((float)rand()/RAND_MAX)*2. - 1.;
	//dump(&GMM[state_i]); */
}

void fgmm_free(struct gmm * gmm)
{
  int i=0;
  for(i=0;i<gmm->nstates;i++)
    gaussian_free(&gmm->gauss[i]);
  free(gmm->gauss);
}
  
/* associate one random data point to 
   a gaussian */
void fgmm_init_random(struct gmm * gmm,
		     const float * data,
		     int data_len)
{
  int state_i =0;
  int i=0;
  int point_idx=0;
  float weights[data_len];
  for(i=0;i<data_len;i++)
    {
      weights[i] = 1.;
    }
 
  smat_covariance(gmm->gauss[0].covar,
		  data_len,
		  weights,
		  data,
		  gmm->gauss[0].mean);

  /*  float xx = 1./gmm->nstates;
      smat_multf(gmm->gauss[0].covar,&xx);*/

  for(;state_i < gmm->nstates;state_i++)
    {
      point_idx = rand()%data_len;
      fgmm_set_mean(gmm,state_i,&data[point_idx*gmm->dim]);
      if(state_i>0) 
	{
	  fgmm_set_covar(gmm,state_i,gmm->gauss[0].covar->_);
	}
      fgmm_set_prior(gmm,state_i,1./gmm->nstates);
    }
}

void fgmm_draw_sample(struct gmm * gmm, float * out)
{
  int st=-1;
  float cumprod=0.;
  float v = ((float)rand())/RAND_MAX;
  while((cumprod < v) && ( st<(gmm->nstates-1)))
    {
      st++;
      cumprod += gmm->gauss[st].prior;
    }
  //printf("%d %f %f\n",st,v,cumprod);
  gaussian_draw(&(gmm->gauss[st]),out);
}


void fgmm_set_prior(struct gmm * gmm,int state, float prior)
{
  gmm->gauss[state].prior = prior;
}

void fgmm_set_mean(struct gmm * gmm,int state, const float * mean)
{
  int i=0;
  for(;i<gmm->dim;i++)
    gmm->gauss[state].mean[i] = mean[i];
}

void fgmm_set_covar(struct gmm * gmm,int state, float * covar)
{
  int i=0;
  for(;i<gmm->gauss[state].covar->_size;i++)
    gmm->gauss[state].covar->_[i] = covar[i];
  invert_covar(&gmm->gauss[state]);
}

void fgmm_dump(struct gmm * gmm)
{
  int state_i=0;
  for(;state_i<gmm->nstates;state_i++)
    {
      printf("Gaussian %d ::\n",state_i);
      dump(&(gmm->gauss[state_i]));
    }
}
  
float fgmm_get_pdf( struct gmm * gmm,
		   float * point)
{
  int state_i = 0;
  float like=0;
  for(;state_i<gmm->nstates;state_i++)
    {
      like += gmm->gauss[state_i].prior * gaussian_pdf(&(gmm->gauss[state_i]),point);
    }
  return like;
}
