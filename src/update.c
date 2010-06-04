#include "fgmm.h"
#include "gaussian.h"

void fgmm_update(struct gmm * gmm, float * data_point)
{
  int state_i=0;
  float pxi[gmm->nstates];
  float nf = 0;
  float eta = .05;
  for(;state_i<gmm->nstates;state_i++)
    {
      pxi[state_i] = gaussian_pdf(&gmm->gauss[state_i], data_point);
      nf += pxi[state_i];
    }
  float weighted_lr;
  for(state_i=0;state_i<gmm->nstates;state_i++)
    {
      weighted_lr = eta * pxi[state_i]/nf;
      gaussian_update(&gmm->gauss[state_i],data_point,weighted_lr);
      invert_covar(&gmm->gauss[state_i]);
    }
}
  
  
