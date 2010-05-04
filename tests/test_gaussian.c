#include "gaussian.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define randn_samples 1000000
#define one_sigma 0.682  /* percent of samples between -sigma and sigma */
#define three_sigma 0.044 /* percent of samples outside 2 sigma */

#define value 0.040071 /* pdf of 3d gaussian mean 0, sigma diagonal 1 @ (1.,1.,1.) */

int main(int argc,char ** argv)
{
  srand(time(NULL));
  printf("gaussian test suite, each test uses %d samples \n",randn_samples);
  int i;
  float vn=0;
  int sig_perc=0, sig3_perc=0;
  printf("box_muller testing :\n");
  for(i=0;i<randn_samples;i++)
    {
      vn = randn_boxmuller();
      if( fabs(vn) < 1.)
	sig_perc++;
      else if( fabs(vn) >= 2.)
	sig3_perc++;
    }
  //printf("%d %d\n",sig_perc,sig3_perc);
  assert( abs(sig_perc - one_sigma*randn_samples) < randn_samples / 400);
  assert( abs(sig3_perc - three_sigma*randn_samples) < randn_samples / 400);
  printf("..pass \n");

  /* ----------------------------------------------- */
  printf("simple gaussian pdf test :\n");
  struct gaussian g;
  
  gaussian_init(&g,3);
  invert_covar(&g);
  assert(g.dim == 3);
  float v[] = {1.,1.,1.};
  float pdf = gaussian_pdf(&g,v);
  assert(fabs(pdf - value) < 1e-5);
  printf("..pass \n");
  
  /* ----------------------------------------------- */
  printf("drawing sample from gaussian test :\n");
  float samp[3];
  float * data;
  float * weights;
  int k=0;
  data = (float *) malloc( sizeof(float) * randn_samples * 3);
  weights = (float *) malloc(sizeof(float) * randn_samples);
  for(i=0;i<randn_samples;i++)
    {
      weights[i] = 1.;
      gaussian_draw(&g,samp);
      //printf("%f %f %f\n",samp[0],samp[1],samp[2]);
      for(k=0;k<3;k++)
	{
	  data[3*i+k] = samp[k];
	}
    }
  struct smat * cv=NULL;
  float mean[3];
  smat_zero(&cv,3);
  smat_covariance(cv,randn_samples,weights,data,mean);
  //  smat_pmatrix(cv);
  /* checking covariance is identity */
  float * pcv = cv->_;
  for(i=0;i<3;i++)
    {
      assert(fabs(mean[i]) < 1e-2);
      assert( fabs(*(pcv++)- 1.) < 2e-2);
      for(k=i+1;k<3;k++)
	assert(fabs(*pcv++) < 1e-2);
      
    }
  printf("..pass\n");

  
  return 1;
}
  
