#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include "regression.h"

void fgmm_regression_init_g(struct gaussian_reg * gr)
{
  int i,j;
  gr->subgauss = (struct gaussian *) malloc(sizeof(struct gaussian));
  gaussian_init(gr->subgauss,gr->input_len);
  smat_cholesky(gr->gauss->covar, gr->gauss->covar_cholesky);
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
  //dump(gr->subgauss);
}


void fgmm_regression_init(struct fgmm_reg * reg)
{
  int state=0;
  for(;state < reg->model->nstates ; state++)
    {
      fgmm_regression_init_g(&reg->subgauss[state]);
    }
}



void fgmm_regression_gaussian(struct gaussian_reg* gr, 
			      const float * inputs,
			      struct gaussian * result)
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
      result->mean[i] = gr->gauss->mean[ gr->output_dim[i]];
      for(;j<gr->input_len;j++)
	{
	  result->mean[i] += gr->reg_matrix[i * gr->input_len + j]*tmp[j];
	}
    }

  // result->covar = gr->gauss->covar - gr->reg_matrix *  (gr->subgauss->covar)^-1 gr->reg_matrix
  for(i=0;i<result->covar->_size;i++)
    {
      result->covar->_[i] = gr->gauss->covar->_[i];
    }

  float element;
  int off,k;
  for(i=0 ; i<gr->output_len ; i++)
  {
    
    for(j=0;j<gr->input_len;j++)
      tmp[j] = gr->reg_matrix[i*gr->input_len+j];

    smat_tforward(gr->subgauss->covar_cholesky,tmp,tmp2);
    smat_tbackward(gr->subgauss->covar_cholesky,tmp2,tmp);
    
    element = 0.;
    off = 0;
    
    for(j=0;j<(i+1);j++)
      {
	for(k=0;k<gr->input_len;k++) // scalar product here 
	  element += gr->reg_matrix[i*gr->input_len + k]*tmp[k];
	// column wise filling .. 
	result->covar->_[i+off] += element;
	off += (gr->input_len - j - 1); 
      }
  }
}

/** use a fgmm_ref struct to perform regression 
 * result and covare stores resulting mean and covariance (covariance 
 * is in magic smat order .. 
 * They must be alloc'd before calling this . 
 *
 * if covar == NULL , don't compute covariance 
 */
void fgmm_regression(struct fgmm_reg * reg, 
		     const float * inputs, // inputs dim (reg->input_len 
		     float * result, // outputs    (reg->output_len) /!\ alloc'd by user
		     float * covar)  // out covar  (reg->output_len ** 2/2)  /!\ alloc'd
{
  float weight = 0;
  float weight2 = 0;
  /*float result[reg->output_len];*/
  //float tmp[reg->output_len];

  struct gaussian loc_model ;
  gaussian_init(&loc_model,reg->output_len);

  float likelihood = 0;
  int state = 0;
  int i=0;
  
  for(i=0;i<reg->output_len;i++)
    result[i] = 0;
  
  if(covar != NULL)
    {
      for(i=0;i<loc_model.covar->_size;i++)
	covar[i] = 0.;
    }

  for(;state<reg->model->nstates;state++)
    {
      weight = gaussian_pdf(reg->subgauss[state].subgauss,inputs);
      fgmm_regression_gaussian(&reg->subgauss[state],inputs,&loc_model);

      for(i=0;i<reg->output_len;i++)
	result[i] += weight*loc_model.mean[i];

      weight2 = weight*weight;
      
      if(covar != NULL)
	{
	  for(i=0;i<loc_model.covar->_size;i++)
	    covar[i] += weight2*loc_model.covar->_[i];
	}

      likelihood += weight;
    }
  assert(likelihood > FLT_MIN);
  for(i=0;i<reg->output_len;i++)
    result[i] /= likelihood;

  gaussian_free(&loc_model);

}


/**
 * should be like fgmm_regression_alloc() 
 *                fgmm_regression_alloc_simple()
 *
 * then fgmm_regression_init() 
 *                         _g() for single gaussian 
 */
void fgmm_regression_alloc(struct fgmm_reg ** regression,
			  struct gmm * gmm,
			  int input_len, int * input_dim,
			  int output_len, int * output_dim)
{

  struct fgmm_reg * reg; 
  reg = (struct fgmm_reg*) malloc(sizeof(struct fgmm_reg)); 

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
void fgmm_regression_alloc_simple(struct fgmm_reg ** regression,
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
  fgmm_regression_alloc(regression,gmm,input_len,inputs,output_len,outs);
}

/*
void fgmm_regression(struct gmm * gmm, 
		    float * input,
		    float * output)
*/

void fgmm_regression_free(struct fgmm_reg ** regression)
{
  struct fgmm_reg * reg = *regression;
  free(reg->input_dim);
  free(reg->output_dim);
  int g=0;
  for(;g<reg->model->nstates;g++)
    {
      if(reg->subgauss[g].reg_matrix != NULL)
	free( reg->subgauss[g].reg_matrix );
    }
  free( reg->subgauss );
  free( reg );
  *regression = NULL;
}
  

/* conditionnal sampling */

void fgmm_regression_sampling(struct fgmm_reg * regression, 
			      const float * inputs,
			      float * output)
{
  float weights[regression->model->nstates];
  float nf=0;
  //float tmp[regression->output_len];
  //float likelihood = 0;
  float acc=0;
  int state = 0;
  //  int i=0;

  float picker = ((float) rand())/RAND_MAX;
  
  for(;state<regression->model->nstates;state++)
    {
      weights[state] = gaussian_pdf(regression->subgauss[state].subgauss,inputs);
      nf += weights[state];
    }

  state = 0;
  printf("%f %f \n",picker,acc);
  while(picker > acc)
    {
      acc += weights[state]/nf;
      state++;
    }
  state--;
  printf("rand state %d\n",state);
  struct gaussian * loc_model = (struct gaussian *) malloc(sizeof(struct gaussian));
  gaussian_init(loc_model,regression->output_len);
  
  fgmm_regression_gaussian(&regression->subgauss[state],inputs,loc_model);

  invert_covar(loc_model);
  gaussian_draw(loc_model,output);
  
  gaussian_free(loc_model);
  free(loc_model);
}
    
