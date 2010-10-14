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
_minline float gaussian_pdf(struct gaussian* g, const float* x)
{
  float dist2;
  float dist = smat_sesq(g->icovar_cholesky,g->mean,x);
  dist *= .5;
  
  dist2 =  expf(-dist)/g->nfactor;
  //dist = 0.2;
  // returning zero here would give weird results for EM 
  /* if( isnan(dist2))
    {
      printf("NaN gaussian pdf .. ");
      printf(" %e %e \n ", dist, g->nfactor);
      }*/
  if(dist2 == 0)
    dist2 = FLT_MIN;
  return dist2;
};


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

#define ranf() ( (float) rand())/RAND_MAX

/** random sample from normal law ( mu = 0, sigma = 1. ) **/ 
_minline float randn_boxmuller( void )
{
  float x1, x2, w;
  do {
    x1 = 2.0 * ranf() - 1.0;
    x2 = 2.0 * ranf() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );

  w = sqrt( (-2.0 * log( w ) ) / w );
  x1 *= w;
  /* x2 *= w */   /* 2nd indpdt gaussian */
  return x1;
};



/** incremental mean/var update */
void gaussian_update(struct gaussian * g, 
		     const float * datapoint, 
		     float learning_rate);
