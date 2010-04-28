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


struct gmm_reg {
  struct gmm * model;
  int * input_dim;
  int * output_dim;
  int input_len;
  int output_len;
  struct gaussian_reg * subgauss;
};


/*
void gmm_regression_init(struct gmm_reg * reg)
{
}

void gmm_regression_init_g(struct gaussian_reg * gr)
{
  gaussian_get_subgauss(gr->gauss,&(gr->subgauss),
			input_dim,input_len);
  // reg_matrix = (Sigma^00)-1 * (Sigma^0i)
  if(reg_matrix != NULL)
    free(reg_matrix);
  reg_matrix =(float*)  malloc(sizeof(float) * input_len * output_len);
  struct smat * fcov = gr->gauss->covar;
  for(j=0;j<output_len;j++)
    {
      for(i=0;i<input_len;i++)
	{
	  reg_matrix[j * input_len + i] = smat_get_value(fcov,output_dim[j],input_dim[i]);
	}
    }
}
*/


/**
 * should be like gmm_regression_alloc() 
 *                gmm_regression_alloc_simple()
 *
 * then gmm_regression_init() 
 *                         _g() for single gaussian 
 */
void gmm_regression_alloc(struct gmm_reg ** regression,
			  struct gmm * gmm,
			  int input_len, int * input_dim,
			  int output_len, int * output_dim)
{

  struct gmm_reg * reg = *regression;
  reg = (struct gmm_reg*) malloc(sizeof(struct gmm_reg));

  int i = 0;
  reg->model = gmm;
  reg->input_len = input_len;
  reg->input_dim = (int*) malloc(sizeof(int)*input_len);
  for(;i<input_len;i++)
    reg->input_dim[i] = input_dim[i];
  reg->output_len = output_len;
  reg->output_dim = (int*) malloc(sizeof(int)*output_len);
  for(i=0;i<output_len;i++)
    reg->output_dim[i] = output_dim[i]; 
}

/**
 * initialise a regression structure , considering that 
 * first input_len dimensions are input and the rest the output
 */
void gmm_regression_alloc_simple(struct gmm_reg ** regression,
				 struct gmm * gmm,
				 int input_len)
{
  int output_len = gmm->dim - input_len;
  int inputs[input_len];
  int outs[output_len];
  int i;
  for(i=0;i<input_len;i++)
    {
      inputs[i] = i;
    }
  for(i=0;i<output_len;i++)
    {
      outs[i] = input_len + i;
    }
  gmm_regression_alloc(regression,gmm,input_len,inputs,output_len,outs);
}

/*
void gmm_regression(struct gmm * gmm, 
		    float * input,
		    float * output)
*/
