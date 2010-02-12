#include "smat.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

 /* tests .. */
#define DIM 100
#define INVERSE_PRECISION 1e-1 
int main(int argc, char ** argv)
{
  struct smat * m1 = NULL;
  struct smat * chol = NULL;
  struct smat * check = NULL;

  float *pm1;

  float a[DIM],b[DIM],c[DIM];
  int i=0,j=0;
  
  printf("testing smat_zero \n");
  smat_zero(&check,DIM);
  smat_zero(&m1,DIM);
  smat_zero(&chol,DIM);
  
  for(i=0;i++;i<check->_size)
    assert(check->_[i] == 0.);

  printf("..pass\n");

  printf("testing smat_identity \n");
  smat_identity(m1);
 
  assert(m1->_[0] == 1.);
  assert(m1->_[m1->_size-1] == 1.);
  assert(m1->_[1] == 0.);
  
  printf("..pass\n");
  
  printf("testing smat_ttmult \n");
  smat_ttmult(m1,check);
  pm1 = m1->_;
  for(i=0;i<m1->dim;i++)
    {
      assert(*(pm1++) == 1.);
      for(j=i+1;j<m1->dim;j++)
	assert((*pm1++) == 0.);
      
    }

  printf("..pass\n");

  printf("testing smat_multv \n");
  for(i=0;i<DIM;i++)
    a[i]= (float) i;
  smat_multv(m1,a,b);
  for(i=0;i<DIM;i++)
    {
      //printf("%f  -  %f\n",a[i],b[i]);
      assert(a[i] == b[i]);
    }
  printf("..pass\n");

  printf("testing cholesky \n");
  float * pmat = m1->_;
  for(i=0;i<DIM;i++)
    {
      for(j=i;j<DIM;j++)
	{
	  
	  *pmat = DIM - (j-i);
	  pmat++;
	}
    }
  //pmatrix(m1);
  smat_cholesky(m1,chol);
  //pmatrix(chol);
  
  smat_ttmult(chol,check);
  //pmatrix(check);

  for(i=0;i<m1->_size;i++)
    assert(check->_[i] == m1->_[i]);

  printf("..pass\n");

  printf("checking inversion :\n");

  smat_multv(m1,a,b);
 
  /*for(i=0;i<DIM;i++)
    printf("%f\n",b[i]);*/

  smat_tforward(chol,b,c);
  /*
  for(i=0;i<DIM;i++)
    printf("%f\n",c[i]); 
  */
  smat_tbackward(chol,c,b);
  for(i=0;i<DIM;i++)
    {
      //printf("%f %f\n",a[i],b[i]); 
      assert(fabs(a[i]-b[i]) < INVERSE_PRECISION);  // see the precision of inverse .. 
      }

   printf("..pass\n");
}
