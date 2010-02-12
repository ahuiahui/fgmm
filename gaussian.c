#include "gaussian.h"
//#include "smat.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>

/* check the inverse covariance computation */ 
/* #define CHECK_INVERSE  */ 

float gaussian_pdf(struct gaussian* g, const float* x)
{
  float dist = 0;
  /* dist = -.5 * (x - mu)^T Sigma^-1 (x-mu) */

  float cdata[g->dim];
  float tmp[g->dim];
  float ivect[g->dim];

  int i=0;
  for(i=0;i<g->dim;i++)
    {
      cdata[i] = x[i] - g->mean[i];
    }

  smat_tforward(g->covar_cholesky,cdata,tmp);
  smat_tforward(g->covar_cholesky,tmp,ivect);
  
  for(i=0;i<g->dim;i++)
    {
      dist += ivect[i]*cdata[i];
    }
  dist *= .5;
  dist =  expf(-dist)/g->nfactor;
  //dist = 0.2;
  // returning zero here would give weird results for EM 
  if(dist == 0)
    dist = FLT_MIN;
  return dist;
}


void dump(struct gaussian* g)
{
  int k=0;
  printf("  prior : %f \n",g->prior);
  printf("  mean : ");
  for(k=0;k<g->dim;k++)
    printf("%f  ",g->mean[k]);
  printf("\n");

  printf("  covariance : ");
  /*for(k=0;k<6;k++)
    printf("%f  ",g->covar[k]);*/
  smat_pmatrix(g->covar);
}

void invert_covar(struct gaussian* g)
{
  float det=1.;
  int i=0,diag=0;
  smat_cholesky(g->covar,g->covar_cholesky);
  for(i=0;i<g->dim;i++)
    {
      det *= g->covar_cholesky->_[diag];
      diag += g->dim - i;
      
    }
  det = det*det;
  g->nfactor = sqrtf( pow(M_PI,g->dim) * det);
}

void init_gaussian(struct gaussian * g,int dim)
{
  int i;
  g->dim = dim;
  g->mean = (float *) malloc(dim* sizeof(float));
  g->covar = NULL;
  g->covar_cholesky = NULL;
  for(i=0;i<dim;i++)
    g->mean[i] = 0.;
  smat_zero(&(g->covar),dim);
  smat_identity(g->covar); // just in case :) 
  smat_zero(&(g->covar_cholesky),dim);
}
/*
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
  }*/
  
  
