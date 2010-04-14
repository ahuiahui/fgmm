#include "gaussian.h"

struct gmm {
  struct gaussian * gauss;
  int nstates;
  int dim;
};

/** 
 * draw one sample from the gmm 
 */
void gmm_draw_sample(struct gmm *, float * out);

/**
 *  alloc all the memory needed by the model 
 *  
 *  all gaussians are init'd to zero mean, unity covariance 
 *  zero prior proba */
void gmm_alloc(struct gmm *, int nstates, int dim);

/** free everything allocated by the above function 
 */
void gmm_free(struct gmm *);
