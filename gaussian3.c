#include "gaussian3.h"
#include "vec3.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>

/* check the inverse covariance computation */ 
/* #define CHECK_INVERSE  */ 

float gaussian_pdf(struct gaussian3d* g, const float* x)
{
  float dist = 0;
  /* dist = -.5 * (x - mu)^T Sigma^-1 (x-mu) */
  float cdata[3];

  diff3(x,g->mean,cdata);

  dist += cdata[0] * (g->icovar[0]*cdata[0] + g->icovar[3]*cdata[1] + g->icovar[5]*cdata[2]);
  dist += cdata[1] * (g->icovar[3]*cdata[0] + g->icovar[1]*cdata[1] + g->icovar[4]*cdata[2]);
  dist += cdata[2] * (g->icovar[5]*cdata[0] + g->icovar[4]*cdata[1] + g->icovar[2]*cdata[2]);

  dist *= .5;
  dist =  expf(-dist)/g->nfactor;
  //dist = 0.2;
  // returning zero here would give weird results for EM 
  if(dist == 0)
    dist = FLT_MIN;
  return dist;
}


void dump(struct gaussian3d* g)
{
  int k=0;
  printf("  prior : %f \n",g->prior);
  printf("  mean : ");
  for(k=0;k<3;k++)
    printf("%f  ",g->mean[k]);
  printf("\n");

  printf("  covariance : ");
  for(k=0;k<6;k++)
    printf("%f  ",g->covar[k]);
  printf("\n");
  /*
  for(k=0;k<6;k++)
  printf("%f  ",g->icovar[k]);*/
  printf("\n");

}

void invert_covar(struct gaussian3d* g)
{
  float det=0;
  det = g->covar[0]*g->covar[1]*g->covar[2] +	\
    g->covar[3]*g->covar[4]*g->covar[5]*2 -	\
    g->covar[0]*g->covar[4]*g->covar[4] -	\
    g->covar[3]*g->covar[3]*g->covar[2] -	\
    g->covar[1]*g->covar[5]*g->covar[5];
  
  assert(det!=0);
  
  g->icovar[0] = (g->covar[1]*g->covar[2] - g->covar[4]*g->covar[4])/det;
  g->icovar[1] = (g->covar[0]*g->covar[2] - g->covar[5]*g->covar[5])/det;
  g->icovar[2] = (g->covar[0]*g->covar[1] - g->covar[3]*g->covar[3])/det;
  g->icovar[3] = (g->covar[5]*g->covar[4] - g->covar[3]*g->covar[2])/det;
  g->icovar[4] = (g->covar[5]*g->covar[3] - g->covar[0]*g->covar[4])/det;
  g->icovar[5] = (g->covar[3]*g->covar[4] - g->covar[5]*g->covar[1])/det;
  
  g->nfactor = M_PI*sqrtf(M_PI*det);

#ifdef CHECK_INVERSE
  
  float val = 0;
  
  val = g->covar[0]*g->icovar[0] + g->covar[3]*g->icovar[3] + g->covar[5]*g->icovar[5];
  
  assert(fabs(val-1.) < 1e-6);
  val = g->covar[0]*g->icovar[3] + g->covar[3]*g->icovar[1] + g->covar[5]*g->icovar[4];
  
  assert(fabs(val) < 1e-6);
  
#endif // checkinverse values 
  
}

void init_random(struct gaussian3d* g)
{
  int k=0;
  for(k=0;k<3;k++)
    {
      g->mean[k] = (float)rand() / RAND_MAX;
      g->covar[k] = 1.;
      g->icovar[k] = 1.;
      g->covar[k+3] = 0.;
      g->icovar[k+3] = 0.;
    }
  invert_covar(g);
}
  
