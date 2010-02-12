#include "smat.h"

/** One gaussian distribution */

struct gaussian{
  float prior; /* prior probability */
  int dim;  /* dimensionality */ 
  float * mean ;
  struct smat * covar; /* covariance matrix */ 
  struct smat * covar_cholesky; /* cache for cholesky decomp of covar */ 
  float nfactor; /* cache for determinant of covar */ 
};   

/** compute the probability density at vector value 
    
    value should be at the same dimension than g->dim */
inline float gaussian_pdf(struct gaussian* g, const float* value);

/** randomize the gaussian */
void init_gaussian(struct gaussian* g,int dim);


void invert_covar(struct gaussian* g);

void dump(struct gaussian* g);
