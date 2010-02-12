#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>  // cholesky needs sqrtf

#include "smat.h"

/* Matrix - Vector multiplication 

 out = m * v 

*/
float smat_multv(const struct smat* m, const float * v,float * out)
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


/* fill the mat with identity  /!\ you must first 
 allocate memory (wigh smat_zero for ex .. ) */
void smat_identity(struct smat * mat)
{
  int i=0;
  float * pmat = mat->_;
  for(;i<mat->dim;i++)
    {
      *pmat = 1.;
      pmat += (mat->dim - i);
    }
}

/* print matrix to screen (for debug purposes ) */ 
void smat_pmatrix(const struct smat* mat)
{
  int i,j;
  int curs;
  float * pmat = mat->_;
  for(i=0;i<mat->dim;i++)
    {
      for(j=0;j<i;j++)
	printf("       ");
      for(;j<mat->dim;j++)
	{
	  printf("%.3f  ",*pmat);
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
  int curs=0;
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
      curs = line + in->dim;
      pout++;
      pin++;
      for(col=line+1;col<in->dim;col++)
	{
	  ts = 0.;
	  for(i=0;i<line;i++)
	    ts += tmp[i][line]*tmp[i][col];
	  tmp[line][col] = *pout = (*pin - ts) / tmp[line][line];
	  //curs += in->dim - (col - line);
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

