#include <cstdlib>

extern "C" {

#include "fgmm.h"

}


class Gmm
{

private :
  struct gmm * c_gmm;
  struct fgmm_reg * c_reg;
  float likelihood;
  int ninput;

public :
  Gmm(int states, int dim)
  {
    //c_gmm = (struct gmm *) malloc(sizeof(struct gmm ));
    fgmm_alloc(&c_gmm,states,dim);
    c_reg = NULL;
  };

  ~Gmm()
  {
    if(c_reg != NULL) 
	fgmm_regression_free(&c_reg);
    fgmm_free(&c_gmm);
  };

  void init(float * data,int len)
  {
    fgmm_init_random(c_gmm,data,len);
  };

  void init(float * data,int len, int dim)
  {
    if(dim != c_gmm->dim) 
      {
	printf("Wrong dimension in init");
	return;
      }
    this->init(data,len);
  };

  void Dump()
  {
    fgmm_dump(this->c_gmm);
  };

  int Em(float * data,int len)
  {
    return fgmm_em(c_gmm,data,len,&likelihood,1e-4);
  };

  int Em(float * data,int len, int dim)
  {
    if(dim != c_gmm->dim) 
      {
	printf("Wrong dimension in EM");
	return -1;
      }
    return this->Em(data,len);
  };

  void SetPrior(int state, float val)
  {
    fgmm_set_prior(this->c_gmm,state,val);
  };
  
  void SetMean(int state, float * mean)
  {
    fgmm_set_mean(this->c_gmm,state,mean);
  };
 
  void SetCovariance(int state, float * covar)
  {
    fgmm_set_covar(this->c_gmm,state,covar);
  };

  void Draw(float * sample) 
  {
    fgmm_draw_sample(this->c_gmm,sample);
  };

  void Draw(float * sample,int dimo)
  {
    this->Draw(sample);
  };

  void InitRegression(int ninput){
    if( c_reg != NULL) 
      fgmm_regression_free(&c_reg);
    this->ninput = ninput;
    fgmm_regression_alloc_simple(&c_reg,c_gmm,ninput);
    fgmm_regression_init(c_reg);
  };

  void DoRegression(const float * input, float * output, float * covar=NULL)
  {
    fgmm_regression(c_reg,input,output,covar);
  };

  void DoSamplingRegression(const float * input, float * output)
  {
    fgmm_regression_sampling(c_reg,input,output);
  };

  void DoRegression(float * input, int dimi,
		    float * output, int dimo)
  {
    /*  if(dimi != this->c_reg->input_len or dimo != this->c_reg->output_len) 
	return;*/
    this->DoRegression(input,output);
  };


  void Update(const float * point,bool wta=false)
  {
    if(wta) 
      fgmm_update_wta(c_gmm,point);
    else 
      fgmm_update(c_gmm,point);
  };

};

