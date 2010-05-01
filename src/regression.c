#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "regression.h"

void gmm_regression_init_g(struct gaussian_reg * gr)
{
  int i,j;
  gr->subgauss = (struct gaussian *) malloc(sizeof(struct gaussian));
  gaussian_init(gr->subgauss,gr->input_len);
  gaussian_get_subgauss(gr->gauss,gr->subgauss,
			gr->input_len,gr->input_dim);
  // reg_matrix = (Sigma^00)-1 * (Sigma^0i)
  if(gr->reg_matrix != NULL)
    free(gr->reg_matrix);
  gr->reg_matrix =(float*)  malloc(sizeof(float) * gr->input_len * gr->output_len);
  struct smat * fcov = gr->gauss->covar;
  for(j=0;j<gr->output_len;j++)
    {
      for(i=0;i<gr->input_len;i++)
	{
	  gr->reg_matrix[j * gr->input_len + i] = smat_get_value(fcov,
							 gr->output_dim[j],
							 gr->input_dim[i]);
	}
    }
  dump(gr->subgauss);
}


void gmm_regression_init(struct gmm_reg * reg)
{
  int state=0;
  for(;state < reg->model->nstates ; state++)
    {
      gmm_regression_init_g(&reg->subgauss[state]);
    }
}



void gmm_regression_gaussian(struct gaussian_reg* gr, float * inputs,float * result)
{
  /*float result[gr->output_len];*/
  int j=0,i=0;
  float tmp[gr->input_len]; 
  float tmp2[gr->input_len]; 
  /* OPT : this computation is also done for the 
     subgauss pdf (ie weight of the gaussian in the regression .. */

  for(;i<gr->input_len;i++)
      tmp[i] = inputs[i] - gr->subgauss->mean[i];

  smat_tforward(gr->subgauss->covar_cholesky,tmp,tmp2);
  smat_tbackward(gr->subgauss->covar_cholesky,tmp2,tmp);

  for(i=0;i<gr->output_len;i++)
    {
      result[i] = gr->gauss->mean[ gr->output_dim[i]];
      for(;j<gr->input_len;j++)
	{
	  result[i] += gr->reg_matrix[i * gr->input_len + j]*tmp[j];
	}
    }
}


void gmm_regression(struct gmm_reg * reg, float * inputs,float * result)
{
  float weight = 0;
  /*float result[reg->output_len];*/
  float tmp[reg->output_len];
  float likelihood = 0;
  int state = 0;
  int i=0;
  for(i=0;i<reg->output_len;i++)
    result[i] = 0;
     
  for(;state<reg->model->nstates;state++)
    {
      weight = gaussian_pdf(reg->subgauss[state].subgauss,inputs);
      gmm_regression_gaussian(&reg->subgauss[state],inputs,tmp);
      for(i=0;i<reg->output_len;i++)
	result[i] += weight*tmp[i];
      likelihood += weight;
    }
  assert(likelihood > FLT_MIN);
  for(i=0;i<reg->output_len;i++)
    result[i] /= likelihood;
}


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

  struct gmm_reg * reg; 
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

  int state=0;
  reg->subgauss = (struct gaussian_reg*) malloc(sizeof(struct gaussian_reg) * reg->model->nstates);
  for(;state < reg->model->nstates ; state++)
    {
      reg->subgauss[state].gauss = &gmm->gauss[state];
      reg->subgauss[state].input_len = input_len;
      reg->subgauss[state].output_len = output_len;
      reg->subgauss[state].output_dim = reg->output_dim;
      reg->subgauss[state].input_dim = reg->input_dim;
      reg->subgauss[state].reg_matrix = NULL;
    }
  *regression = reg;
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
