#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>  // cholesky needs sqrtf

#include "smat.h"

/* Matrix - Vector multiplication 

 out = m * v 

*/
void smat_multv(const struct smat* m, const float * v,float * out)
{
  float * pcoef = m->_;
  int i,j;
  for(i=0;i<m->dim;i++)
    out[i] = 0;
  for(i=0;i<m->dim;i++)
    {
      for(j=i;j<m->dim;j++)
	{
	  out[i] += *pcoef * v[j];
	  if(j>i)
	    out[j] += *pcoef * v[i];
	  pcoef++;
	}
    }
}

/* scales Matrix with a float :: 

   m = m*f 
*/

void smat_multf(struct smat* m,const float *f)
{
  int i=0;
  for(i=0;i<m->_size;i++)
    m->_[i] *= *f;
}


/* zero the matrix and does the memory initialisation */
void smat_zero(struct smat ** mat,int dim)
{
  struct smat * m = *mat;
  if(m==NULL)
    {
      m = (struct smat *) malloc ( sizeof(struct smat));
      m->dim = dim;
      m->_size = dim*(dim+1)/2;
      m->_ = (float *) malloc(sizeof(float)*m->_size);
      *mat = m;
    }
  int i;
  for(i=0;i<m->_size;i++)
    m->_[i] = 0.;
}

float smat_get_value(struct smat * mat,int row, int col)
{
  int tmp;
  assert((row < mat->dim ) && (col < mat->dim));
  if(row > col)
    {
      tmp = row;
      row = col;
      col = tmp;
    }
  int i=0;
  int idx=0;
  for(;i<row;i++)
    {
      idx += mat->dim - i;
      col -= (i+1);
    }
  return mat->_[idx+col];
}

void smat_get_submatrix(struct smat * mat , struct smat * res,
			 int n_dims, int * dims )
{
  int i=0;
  int j;
  float *pres = res->_;

  for(;i<n_dims;i++)
    {
      for(j=i;j<n_dims;j++)
	{
	  *pres = smat_get_value(mat,dims[i],dims[j]);
	  pres++;
	}
    }
}
	    


void smat_free(struct smat ** mat)
{
  free( (*mat)->_ );
  free( *mat );
  *mat = NULL;
}

/* fill the mat with identity  /!\ you must first 
 allocate memory (wigh smat_zero for ex .. ) */
void smat_identity(struct smat * mat)
{
  int i=0;
  int j;
  float * pmat = mat->_;
  for(;i<mat->dim;i++)
    {
      (*pmat++) = 1.;
      for(j=i+1;j<mat->dim;j++)
	(*pmat++) = 0.;
    }
}

void smat_add_diagonal(struct smat * mat, float value)
{
  int i=0;
  int j;
  float * pmat = mat->_;
  for(;i<mat->dim;i++)
    {
      *pmat += value;      
      pmat += (mat->dim - i);
    }
}
  

/* print matrix to screen (for debug purposes ) */ 
void smat_pmatrix(const struct smat* mat)
{
  int i,j;
  float * pmat = mat->_;
  for(i=0;i<mat->dim;i++)
    {
      for(j=0;j<i;j++)
	printf("       ");
      for(;j<mat->dim;j++)
	{
	  printf("%e  ",*pmat);
	  pmat++;
	}
      printf("\n");
    }
}


/* Cholesky decomposition 
  
   out is a UPPER triang matrix such as  out^T * out = in 
   
   will fail if in is not SDP */

void smat_cholesky(const struct smat* in,struct smat* out)
{
  assert(in->dim == out->dim);
  float tmp[in->dim][in->dim];
  int line=0;
  int col=0;
  int i;
  float ts=0;
  float *pout = out->_;
  float *pin = in->_;
  for(line=0;line<in->dim;line++)
    {
      ts = 0;
      for(i=0;i<line;i++)
	  ts += tmp[i][line]*tmp[i][line];
      assert((*pin - ts) > 0.);
      tmp[line][line] = *pout = sqrtf( *pin - ts);
      pout++;
      pin++;
      for(col=line+1;col<in->dim;col++)
	{
	  ts = 0.;
	  for(i=0;i<line;i++)
	    ts += tmp[i][line]*tmp[i][col];
	  tmp[line][col] = *pout = (*pin - ts) / tmp[line][line];
	  pin++;
	  pout++;
	}    
    }
}

/* L^T * L  for a triang SUP matrix  */
void smat_ttmult(const struct smat* tri, struct smat* out)
{
  int uidx=0;
  int didx=0;
  int oidx=0;
  int line=0;
  smat_zero(&out,tri->dim);
  int linend = tri->dim - 1;

  for(didx=0;didx<tri->_size;didx++)
    {
      
      for(uidx=didx;uidx <= linend ;uidx++)
	{
	  //printf("%d %d %d\n",oidx,didx,uidx);
	  out->_[oidx] += tri->_[uidx] * tri->_[didx];
	  oidx++;
	}
      if(didx == linend)
	{
	  linend += (tri->dim - (line+1));
	  line++;
	  oidx = didx + 1; 
	  //printf("%d %d\n",oidx,linend);
	}
    }
}

/** returns (x - bias)' (ichol^T ichol)^-1 ( x - bias ) 
 *  ichol is the cholesky decomposition of Sigma, with inverted diagonal 
 */

float smat_sesq(struct smat * ichol,const float * bias,const float * x)
{
  float out = 0.;
  int i,j;
  float cdata[ichol->dim];
  float * pichol = ichol->_;
  for(i=0;i<ichol->dim;i++)
    cdata[i] = 0.;      
  for(i=0;i<ichol->dim;i++)
    {
      cdata[i] += x[i] - bias[i];
      cdata[i] *= *pichol++;
      for(j=i+1;j<ichol->dim;j++)
	{
	  cdata[j] -= (*pichol++)*cdata[i];
	}
      out += cdata[i]*cdata[i];
    }
  return out;
}
  

/* resolve  L*y = b 
   where L is * LOWER *  triangular */ 

void smat_tforward(struct smat * lower, float * b, float * y) 
{
  int i,j;
  float * pL = lower->_;
  for(i=0;i<lower->dim;i++)
    y[i] = b[i];
  for(i=0;i<lower->dim;i++)
    {
      //y[i] = b[i];
      y[i] /= (*pL++);
      for(j=i+1;j<lower->dim;j++)
	{
	  y[j] -= (*pL++)*y[i];
	}
      
    }
}

void smat_as_square(const struct smat * mat, float * square) 
{
  int i,j;
  float * pmat = mat->_;
  for(i=0;i<mat->dim;i++) 
    {
      square[i*mat->dim + i] = (*pmat++);
      for(j=i+1;j<mat->dim;j++)
	{
	  square[i*mat->dim + j] = *pmat;
	  square[j*mat->dim + i] = *pmat;
	  pmat++;
	}
    }
}

void smat_from_square(struct smat * mat, const float * square)
{
  int i,j;
  float * pmat = mat->_;
  for(i=0;i<mat->dim;i++) 
    {
      (*pmat++) = square[i*mat->dim + i];
      for(j=i+1;j<mat->dim;j++)
	{
	  *pmat = square[i*mat->dim + j];
	  pmat++;
	}
    }
}

/* resolve L*y = b 
   where L is upper triangular (like the result of cholesky ..) */

void smat_tbackward(const struct smat * upper, float * b, float * y)
{
  int i,j;
  float * pU = upper->_ + upper->_size -1; // points to the end  
  for (i = upper->dim - 1; i >= 0; i--)
  {
    y[i] = b[i];
    
    for (j = upper->dim -1 ; j > i ; j--)
    {
      y[i] -= (*pU--) * y[j];
    }
    assert(*pU != 0.);
    y[i] /= (*pU--);
  }
}

float smat_covariance(struct smat * cov, 
		     int ndata, 
		     const float * weight,
		     const float * data,
		     float * mean)
{
  const float * pdata = data;
  const float * pweight = weight;
  float * pcov = cov->_;
  float cdata[cov->dim]; // fixme
  int i=0,j=0,k=0;
  smat_zero(&cov,cov->dim);
  float norm=0;
  for(i=0;i<cov->dim;i++)
    {
      mean[i] = 0.;
    }
  for(i=0;i<ndata;i++)
    {
      for(j=0;j<cov->dim;j++)
	mean[j] += (*pweight)*(*pdata++);
      norm += *pweight;
      pweight++;
    }
  for(i=0;i<cov->dim;i++)
    {
      mean[i] /= norm;
    }
  pdata = data;
  pweight = weight;
  for(i=0;i<ndata;i++)
    {      
      pcov = cov->_;
      for(j=0;j<cov->dim;j++)
	{
	  cdata[j] = (*pdata++) - mean[j];
	}

      j=0;
      k=0;
      while(j<cov->dim)
	{
	  *pcov += (*pweight)*cdata[j]*cdata[k];
	  pcov++;
	  k++;
	  if(k==cov->dim) // next line
	    k=++j;
	}
      pweight++;
    }
  for(i=0;i<cov->_size;i++)
    cov->_[i] /= norm;
  return norm;
}

float smat_covariance_diag(struct smat * cov, 
			   int ndata, 
			   const float * weight,
			   const float * data,
			   float * mean)
{
  const float * pdata = data;
  const float * pweight = weight;
  float * pcov = (float *) malloc(sizeof(float) * cov->dim);
  float cdata[cov->dim];
  float *pmat = cov->_;

  int i=0,j=0,k=0;
  smat_zero(&cov,cov->dim);
  float norm=0;
  for(i=0;i<cov->dim;i++)
    {
      mean[i] = 0.;
      pcov[i] = 0.;
    }
  for(i=0;i<ndata;i++)
    {
      for(j=0;j<cov->dim;j++)
	mean[j] += (*pweight)*(*pdata++);
      norm += *pweight;
      pweight++;
    }
  for(i=0;i<cov->dim;i++)
    {
      mean[i] /= norm;
    }
  pdata = data;
  pweight = weight;
  for(i=0;i<ndata;i++)
    {      
      pcov = cov->_;
      for(j=0;j<cov->dim;j++)
	{
	  cdata[j] = (*pdata++) - mean[j];
	}

      j=0;
      /*k=0;*/
      for(j=0;j<cov->dim;j++)
	{
	  pcov[j] += (*pweight) * cdata[j]*cdata[j];
	  /* *pcov += (*pweight)*cdata[j]*cdata[k]; */
	  /* pcov++; */
	  /* k++; */
	  /* if(k==cov->dim) // next line */
	  /*   k=++j; */
	}
      pweight++;
    }

  for(k=0;k<cov->dim;k++)
    {
      (*pmat++) = pcov[k];
      for(j=k+1;j<cov->dim;j++)
	(*pmat++) = 0.;
    }

  for(i=0;i<cov->_size;i++)
    cov->_[i] /= norm;
  return norm;
}
