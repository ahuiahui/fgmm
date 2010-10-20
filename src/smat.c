/************************************************************************/
/* This file is part of libfgmm.				        */
/* 								        */
/* libfgmm is free software: you can redistribute it and/or modify      */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.				        */
/* 								        */
/* libfgmm is distributed in the hope that it will be useful,	        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU Lesser General Public License for more details.			        */
/* 								        */
/* You should have received a copy of the GNU Lesser General Public License    */
/* along with libfgmm.  If not, see <http://www.gnu.org/licenses/>.     */
/* 								        */
/* Copyright 2010        LASA  - EPFL   http://lasa.epfl.ch             */
/*                                                                      */
/*       Florent D'halluin   <florent.dhalluin@epfl.ch>		        */
/************************************************************************/

#include <stdio.h>
#include <assert.h>

#include "smat.h"

/* zero the matrix and does the memory initialisation */
void smat_zero(struct smat ** mat,int dim)
{
  struct smat * m = *mat;
  int i;
  if(m==NULL)
    {
      m = (struct smat *) malloc ( sizeof(struct smat));
      m->dim = dim;
      m->_size = dim*(dim+1)/2;
      m->_ = (float *) malloc(sizeof(float)*m->_size);
      *mat = m;
    }
  for(i=0;i<m->_size;i++)
    m->_[i] = 0.;
}

float smat_get_value(struct smat * mat,int row, int col)
{
  int tmp;
  int i=0;
  int idx=0;
  assert((row < mat->dim ) && (col < mat->dim));
  if(row > col)
    {
      tmp = row;
      row = col;
      col = tmp;
    }
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
   
   returns 1 on success, 0 on failure (like the matrix not 
   being strictly positive or full ranked */

int smat_cholesky(const struct smat* in,struct smat* out)
{
  float * tmp; //[in->dim][in->dim];
  int line=0;
  int col=0;
  int i;
  float ts=0;
  float *pout = out->_;
  float *pin = in->_;
  
  assert(in->dim == out->dim);
  tmp = (float *) malloc(sizeof(float) * in->dim * in->dim);
  
  for(line=0;line<in->dim;line++)
    {
      ts = 0;
      for(i=0;i<line;i++)
	  ts += tmp[i*in->dim + line]*tmp[i*in->dim + line];

      if((*pin - ts) <= 0.) // positiveness check
	{
	  free(tmp);
	  return 0.;
	}

      tmp[line*in->dim + line] = *pout = sqrtf( *pin - ts);
      pout++;
      pin++;
      for(col=line+1;col<in->dim;col++)
	{
	  ts = 0.;
	  for(i=0;i<line;i++)
	    ts += tmp[i*in->dim + line]*tmp[i*in->dim + col];
	  tmp[line*in->dim + col] = *pout = (*pin - ts) / tmp[line*in->dim + line];
	  pin++;
	  pout++;
	}    
    }
  free(tmp);
  return 1;
}

/* L^T * L  for a triang SUP matrix  */
void smat_ttmult(const struct smat* tri, struct smat* out)
{
  int uidx=0;
  int didx=0;
  int oidx=0;
  int line=0;
  int linend = tri->dim - 1;
  smat_zero(&out,tri->dim);
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

 /*
float smat_sesq(struct smat * ichol,const float * bias,const float * x)
{
  float out = 0.;
  int i,j;
  float * cdata; //[ichol->dim];
  float * pichol = ichol->_;
  
  cdata = (float *) malloc(sizeof(float) * ichol->dim);
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
  free(cdata);
  return out;
}
  */

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
  float * cdata;
  int i=0,j=0,k=0;
  float norm=0;
  
  smat_zero(&cov,cov->dim);
  cdata = (float *) malloc(sizeof(float) * cov->dim);
  
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
  free(cdata);
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
  float * pcov;
  float *pmat = cov->_;
  float tmp;

  int i=0,j=0,k=0;
  float norm=0;
  
  smat_zero(&cov,cov->dim);
  pcov = (float *) malloc(sizeof(float) * cov->dim);
  
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
      for(j=0;j<cov->dim;j++)
	{
	  tmp = (*pdata++) - mean[j];
	  pcov[j] += (*pweight) * tmp*tmp;
	}
      pweight++;
    }

  for(k=0;k<cov->dim;k++)
    {
      (*pmat++) = pcov[k] / norm;
      for(j=k+1;j<cov->dim;j++)
	(*pmat++) = 0.;
    }

  free(pcov);
  return norm;
}

float smat_covariance_single(struct smat * cov, 
			     int ndata, 
			     const float * weight,
			     const float * data,
			     float * mean)
{
  const float * pdata = data;
  const float * pweight = weight;
  float tmp;
  float *pmat = cov->_;
  float total_mean;
  float variance;
  int i=0,j=0,k=0;
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
  total_mean = 0.;
  for(i=0;i<cov->dim;i++)
    {
      mean[i] /= norm;
      total_mean += mean[i];
    }

  total_mean /= cov->dim;

  variance = 0.;

  pdata = data;
  pweight = weight;
  for(i=0;i<ndata;i++)
    {      
      for(j=0;j<cov->dim;j++)
	{
	  tmp = ((*pdata++) - total_mean);
	  variance += *pweight * tmp * tmp;
	}
      pweight++;
    }

  variance /= (norm * cov->dim);
  for(k=0;k<cov->dim;k++)
    {
      (*pmat++) = variance;
      for(j=k+1;j<cov->dim;j++)
	(*pmat++) = 0.;
    }
  return norm;
}
