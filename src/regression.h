#include "gaussian.h"
#include "gmm.h"

struct gaussian_reg {
  struct gaussian * gauss;
  struct gaussian * subgauss; // input subgaussian Used to compute the weight of this
  int * input_dim;
  int * output_dim;
  int input_len;
  int output_len;
  float * reg_matrix; // store in->out A matrix 
};


struct fgmm_reg {
  struct gmm * model;
  int * input_dim;
  int * output_dim;
  int input_len;
  int output_len;
  struct gaussian_reg * subgauss;
};

