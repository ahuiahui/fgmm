#include "gmm.h"
#include <stdlib.h>

void gmm_alloc(struct gmm * gm,int nstates,int dim)
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

void gmm_free(struct gmm * gmm)
{
  int i=0;
  for(i=0;i<gmm->nstates;i++)
    gaussian_free(&gmm->gauss[i]);
  free(gmm->gauss);
}
    

void gmm_draw_sample(struct gmm * gmm, float * out)
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
