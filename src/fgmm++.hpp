/**
 * ------------
 * FGMM library 
 * ------------
 *
 * a fast(er) and light Gaussian mixture model implementation. 
 *        C++ bindings 
 *
 * Florent D'halluin <florent.dhalluin@epfl.ch> 
 */


#include <cstdlib>

extern "C" {

#include "fgmm.h"

}

/**
 * Gaussian Mixture Model class 
 */ 
class Gmm
{
public :

  /**
   * @param states : fixed number of states in the model 
   * @param dim    : dimensionnality of the space we are working in
   */

  int dim;
  int ninput;
  int nstates;

  Gmm(int states, int dim)
  {
    //c_gmm = (struct gmm *) malloc(sizeof(struct gmm ));
    fgmm_alloc(&c_gmm,states,dim);
    c_reg = NULL;
    this->dim = dim;
    this->ninput = 0;
    this->nstates = states;
  };

  ~Gmm()
  {
    if(c_reg != NULL) 
	fgmm_regression_free(&c_reg);
    fgmm_free(&c_gmm);
  };

  /**
   * call this before any kind of learning .. 
   * set means of gaussians by picking random points in 
   * the dataset, and set the variance using the variance 
   * of the dataset
   *
   * @param len : # of points in the dataset
   * @param data : dim*len array, datapoints. (row order) 
   */
  void init(float * data,int len)
  {
    fgmm_init_random(c_gmm,data,len);
  };

  void initKmeans(float * data,int len)
  {
    fgmm_init_kmeans(c_gmm,data,len);
  };


  /**
   * Just print the model's parameter on stdout
   */
  void Dump()
  {
    fgmm_dump(this->c_gmm);
  };

  /**
   * Expectation Maximization Algorithm. 
   */
  int Em(float * data,int len, float epsilon=1e-4)
  {
    return fgmm_em(c_gmm,data,len,&likelihood,epsilon,NULL);
  };


  float Pdf(float * obs, float * weights=NULL)
  {
    return fgmm_get_pdf(c_gmm,obs,weights);
  };

  /**
   * set Prior probability for the desired state
   */
  void SetPrior(int state, float val)
  {
    fgmm_set_prior(this->c_gmm,state,val);
  };
  
  /**
   * set the mean of the specified state, 
   *
   * @param state : the state index 
   * @param mean : an array of size dim, specify the mean
   */
  void SetMean(int state, float * mean)
  {
    fgmm_set_mean(this->c_gmm,state,mean);
  };

  /**
   * set the covariance of the specified state 
   *
   * @param  state : the state index
   * @param  covar : covariance matrix 
   * @param  AsSymetric : Using symetric matrix 
   *                 order .. dim*(dim+1)/2 
   * 
   *      Symetric matrix form : 
   *         [[ 1 2 3 4 ] 
   *          [ 2 5 6 7 ]
   *          [ 3 6 8 9 ]
   *          [ 4 7 9 10]]  
   *  
   * if not we are using a standart row order . 
   */
 
  void SetCovariance(int state, float * covar, bool AsSymetric=true)
  {
    if(AsSymetric) 
      fgmm_set_covar_smat(this->c_gmm,state,covar);
    else
      fgmm_set_covar(this->c_gmm,state,covar);
  };


  float GetPrior(int state)
  {
    return fgmm_get_prior(this->c_gmm,state);
  };

  void GetMean(int state, float * output)
  {
    float * pMean = fgmm_get_mean(this->c_gmm,state);
    for(int i=0;i<this->c_gmm->dim;i++)
      output[i] = pMean[i];
  }

  void GetCovariance(int state, float * out,bool AsSymetric=false)
  {
    if(!AsSymetric)
      {
	fgmm_get_covar(this->c_gmm,state,out);
      }
    else 
      {
	float * pC = fgmm_get_covar_smat(this->c_gmm,state);
	for(int i=0;i<this->c_gmm->dim*(this->c_gmm->dim+1)/2;
	    i++)
	  out[i] = pC[i];
      }
  }
  /**
   * draw a random sample from the model
   *
   * @param sample : the output sample, must be alloc'd of 
   *                 size dim. 
   */
  void Draw(float * sample) 
  {
    fgmm_draw_sample(this->c_gmm,sample);
  };


  /**
   * Initilization for Gaussian Mixture Regression 
   *
   * @param ninput : the first ninput dimension are the inputs, 
   *                 remaining dimensions are outputs. 
   */
  void InitRegression(int ninput){
    if( c_reg != NULL) 
      fgmm_regression_free(&c_reg);
    this->ninput = ninput;
    fgmm_regression_alloc_simple(&c_reg,c_gmm,ninput);
    fgmm_regression_init(c_reg);
  };

  /**
   * Perform the regression on one input point : 
   *
   * @param input : the input point (array of ninput) 
   * @param output : alloc'd array to store result. 
   * @param covar : eventually store resulting covariance is symetric matrix
   *                order (set SetCovariance ) 
   */
  void DoRegression(const float * input, float * output, float * covar=NULL)
  {
    fgmm_regression(c_reg,input,output,covar);
  };


  /**
   * Conditional sampling from the model : 
   * draw a sample in output subspace given the input point in 
   * input subspace. 
   * you must call InitRegression before. 
   */
   void DoSamplingRegression(const float * input, float * output)
   {
     fgmm_regression_sampling(c_reg,input,output);
   };


  /**
   * Online learning HIGHLY EXPERIMENTAL :: 
   * 
   * @param point : input point 
   * @param wta   : use winner take all update ( only update 
   *                 Most likely gaussian) 
   */
  void Update(const float * point,bool wta=false)
  {
    if(wta) 
      fgmm_update_wta(c_gmm,point);
    else 
      fgmm_update(c_gmm,point);
  };


  int GetLikelyState(const float * point)
  {
    return fgmm_most_likely_state(this->c_gmm,point);
  };
  
private :
  struct gmm * c_gmm;
  struct fgmm_reg * c_reg;
  float likelihood;

};

