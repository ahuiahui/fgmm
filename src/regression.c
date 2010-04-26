
struct gmm_reg {
  struct gmm * model;
  int * input_dim;
  int * output_dim;
  int input_len;
  int output_len;
};

void gmm_init_reg(struct gmm_reg * reg,
		  struct gmm * gmm,
		  int input_len, int * input_dim,
		  int output_len, int * output_dim)
{
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
void gmm_init_reg_simple(struct gmm_reg * reg,
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
  gmm_init_reg(reg,gmm,input_len,inputs,output_len,outs);
}

void gmm_regression(struct gmm * gmm, 
		    float * input,
		    float * output)
