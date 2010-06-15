#include "gaussian.h"
//#include "smat.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <assert.h>


#define ranf() ( (float) rand())/RAND_MAX

/* fast normal law generator */
float randn_boxmuller()
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
}

/* check the inverse covariance computation */ 
/* #define CHECK_INVERSE  */ 

float gaussian_pdf(struct gaussian* g, const float* x)
{
  
  
  /* dist = -.5 * (x - mu)^T Sigma^-1 (x-mu) */
  /*
  float dist = 0;
  float cdata[g->dim];
  float tmp[g->dim];
  float ivect[g->dim];

  int i=0;
  for(i=0;i<g->dim;i++)
    {
      cdata[i] = x[i] - g->mean[i];
    }

  smat_tforward(g->covar_cholesky,cdata,tmp);
  smat_tbackward(g->covar_cholesky,tmp,ivect);
  
  for(i=0;i<g->dim;i++)
    {
      dist += ivect[i]*cdata[i];
    }
  */
  
  float dist = smat_sesq(g->icovar_cholesky,g->mean,x);
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
  int i=0,j=0;
  smat_cholesky(g->covar,g->covar_cholesky);
  float * pichol = g->icovar_cholesky->_;
  float * chol = g->covar_cholesky->_;

  for(i=0;i<g->dim;i++)
    {
      det *= *chol;
      *pichol = 1./(*chol);

      chol++;
      pichol++;

      for(j=i+1;j<g->dim;j++)
	{
	  *chol++ = *pichol++;
	}
    }
  det = det*det;
  g->nfactor = sqrtf( pow(M_PI,g->dim) * det);
}

void gaussian_init(struct gaussian * g,int dim)
{
  int i;
  g->dim = dim;
  g->mean = (float *) malloc(dim * sizeof(float));
  g->covar = NULL;
  g->covar_cholesky = NULL;
  g->icovar_cholesky = NULL;
  smat_zero(&(g->icovar_cholesky),dim);
  for(i=0;i<dim;i++)
    g->mean[i] = 0.;
  smat_zero(&(g->covar),dim);
  smat_identity(g->covar); // just in case :) 
  smat_zero(&(g->covar_cholesky),dim);
  invert_covar(g);
}

void gaussian_free(struct gaussian * g)
{
  free(g->mean);
  smat_free(&g->covar);
  smat_free(&g->covar_cholesky);
  smat_free(&g->icovar_cholesky);
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
  

void gaussian_draw(struct gaussian * g, float * out)
{
  int i=0;
  float tvec[g->dim];
  for(;i<g->dim;i++)
    tvec[i] = randn_boxmuller();
  smat_multv(g->covar_cholesky,tvec,out);
  for(i=0;i<g->dim;i++)
    out[i] += g->mean[i];
}

void gaussian_get_subgauss(struct gaussian* g, struct gaussian* result,
			   int n_dim, int * dims)
{
  if(result->dim != n_dim)
    {
      gaussian_free(result);
      gaussian_init(result,n_dim);
    }
  smat_get_submatrix(g->covar,result->covar,n_dim,dims);
  int i=0;
  for(;i<n_dim;i++)
    result->mean[i] = g->mean[dims[i]];
  invert_covar(result);
}
      
void gaussian_update(struct gaussian * g,
		     const float * data, 
		     float lr)
{
  int i=0;
  int j=0;
  int curs=0;
  for(;i<g->dim;i++)
    {
      g->mean[i] += lr*(data[i] - g->mean[i]);
      for(j=i;j<g->dim;j++)
	{
	  g->covar->_[curs] += lr*( (data[i]-g->mean[i])*(data[j] - g->mean[j]) - g->covar->_[curs]) ;
	  curs++;
	}
    }
}
	    
