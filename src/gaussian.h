#include "smat.h"

/** One gaussian distribution */

struct gaussian{
  float prior; /* prior probability */
  int dim;  /* dimensionality */ 
  float * mean ;
  struct smat * covar; /* covariance matrix */ 
  struct smat * covar_cholesky; /* cache for cholesky decomp of covar */ 
  struct smat * icovar_cholesky; /* cholesky matrix with inverse diagonal */
  float nfactor; /* cache for determinant of covar */ 
};   

/** compute the probability density at vector value 
    
    value should be at the same dimension than g->dim */
inline float gaussian_pdf(struct gaussian* g, const float* value);

/** alloc memory for the gaussian 
    and init it to zero with identity covariance matrix
*/
void gaussian_init(struct gaussian* g,int dim);
void gaussian_free(struct gaussian* g);

void invert_covar(struct gaussian* g);

void dump(struct gaussian* g);

/* draw one sample from the gaussian */
void gaussian_draw(struct gaussian* g, float * out);

/* get the projection of the gaussian on the given dimensions 
 * if result in NULL or wrong dimension .. is it (re) alloc'd */
void gaussian_get_subgauss(struct gaussian* g, struct gaussian* result,
			   int n_dim, int * dims);

/** random sample from normal law ( mu = 0, sigma = 1. ) **/ 
inline float randn_boxmuller();

/** incremental mean/var update */
void gaussian_update(struct gaussian * g, float * datapoint, float learning_rate);
