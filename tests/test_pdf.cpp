#include "fgmm++.hpp"
#include <iostream>
#include <assert.h>
#include <math.h>
#include <time.h>

/** 
 * test case for non diagonal covariance matrices, (and c++ ) 
 *
 * generates points using a 1 state gmm and look at their 
 * statistical properties 
 */

int main(int argc, char *argv[])
{
  srand(time(0));
  Gmm gmm(1,2); // one state 2 d 
  float cov[3] = { 10 , -4. , 2.};
  gmm.SetCovariance(0,cov);
  
  float sample[2];
  float mean = 0.;
  int count = 0;
   
  while(count < 500)
    {
      gmm.Draw(sample);
      if(sample[0] > 5.)
	{
	  mean += sample[1];
	  count++;
	}
      //std::cout << sample[0] << "  " << sample[1] << std::endl;
    }
  
  mean /= count;
  assert( fabs(mean + 2.55) < 1e-1);
  std::cout << "pdf test Pass : " << mean << std::endl;
  return 0;
}
