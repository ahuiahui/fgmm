#include "vec3.h"

void diff3(const float* x,const float* y,float* r)
{
  r[0] = x[0] - y[0];
  r[1] = x[1] - y[1];  
  r[2] = x[2] - y[2];
}

void mmul3(float *x, float v)
{
  x[0]*=v;
  x[1]*=v;
  x[2]*=v;
}
 
