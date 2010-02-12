#include "gaussian.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define value 0.040071

int main(int argc,char ** argv)
{
  printf("simple gaussian pdf test :\n");
  struct gaussian g;
  
  init_gaussian(&g,3);
  invert_covar(&g);
  float v[] = {1.,1.,1.};
  float pdf = gaussian_pdf(&g,v);
  assert(fabs(pdf - value) < 1e-5);
  printf("..pass \n");
  return 1;
}
  
