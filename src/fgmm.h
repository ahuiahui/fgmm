
/**
 * structure holding one state of the  gmm 
 *
 * see gaussian.h for more details .. 
 */

struct gaussian;

/**
 * structure holding the whole Gaussian Mixture model 
 * this must be initialised by fgmm_alloc 
 *
 * short example  : 
 *   struct gmm mygmm;
 *   fgmm_alloc(&mygmm,4,1) // 1D 4 states GMM 
 *   int i = 0;
 *   for(i=0;i<4;i++)
 *     { 
 *       fgmm_set_prior(&mygmm,i,.25); 
 *        fgmm_set_mean(&mygmm,i,i*2.5):
 *        fgmm_set_covar(&mygmm,i,.1):
 *      }
 *    float sample; 
 *    fgmm_draw_sample(&mygmm,&sample);
 *    fgmm_free(&mygmm);
 */
       
struct gmm {
  struct gaussian * gauss;
  int nstates;
  int dim;
};

/**
 *  alloc all the memory needed by the model (alloc all gaussians struct and all ) 
 *  
 *  all gaussians are init'd to zero mean, unity covariance 
 *  zero prior proba
 *
 *  @param **gmm : pointer to an gmm struct  
 *  @param nstates : number of states 
 *  @param dim : dimensionnality of input space 
 */
void fgmm_alloc(struct gmm ** gmm, int nstates, int dim);

/**
 * free everything allocated by the above function 
 */
void fgmm_free(struct gmm ** gmm);

/**
 * initialize the model from the data by :
 *  for each gaussian :
 *     - pick one random data point from data 
 *     - set the mean to this point 
 *     - set covariance to data covariance/nstates
 *     - set prior to 1./nstates
 *
 * the model must be first alloc'd (with fgmm_alloc) 
 */
void fgmm_init_random(struct gmm * gmm,
		     const float * data,
		     int data_len);

/**
 * set the prior of a given state 
 *
 * @param state : index of the state
 * @param prior : the prior value (be careful to keep sum(priors) = 1. )
 */
void fgmm_set_prior(struct gmm *,int state, float prior);

/**
 * set the mean of a given state 
 *
 * @param state : index of the state
 * @param mean : table of dim float 
 */

void fgmm_set_mean(struct gmm *,int state, const float * mean);

/**
 * Set the covariance of a given state 
 *
 * @param state : index of the state 
 * @param covar : The covariance matrix 
 * 
 *      Symetric matrix form : 
 *         [[ 1 2 3 4 ]
 *          [ 2 5 6 7 ]
 *          [ 3 6 8 9 ]
 *          [ 4 7 9 10]]  
 *
 *  covar must be a float table of lenght (dim*dim+1)/2 
 */

void fgmm_set_covar(struct gmm *,int state, float * covar);

/**
 * print the gmm parameters to screen 
 */

void fgmm_dump(struct gmm * gmm);

/**
 * draw one sample from the gmm 
 *
 * @param out : *alloc'd * table of dim 
 *
 */
void fgmm_draw_sample(struct gmm *, float * out);


#define loglikelihood_eps 1e-4

/**
 * EM algorithm 
 *
 * @param GMM : the initialized gmm 
 * @param data : a dim x data_length float table holding training data (row order) 
 *               the k-th training point is data[dim*k]..data[dim*(k+1)] 
 * @param data_length : number of training points 
 * @param end_loglikelihood : will be set to the final loglikelihood 
 * @param likelihood_epsilon : if the loglikelihood variation is below this threshold, stops here
 * @return : The number of iterations 
 */
int fgmm_em( struct gmm * GMM,
	     const float * data,
	     int data_length, 
	     float * end_loglikelihood,
	     float likelihood_epsilon);

/**
 * return likelihood of point
 */
float fgmm_get_pdf( struct gmm * gmm,
		   float * point);


/**
 * Structure holding stuffs for the regression 
 * 
 * This must be alloc'd using fgmm_regression_alloc* on a 
 *  non allocated pointer. 
 */

struct fgmm_reg;

/**
 * init a fgmm_reg structure for a regression 
 * were the input_len first input dimenstion are 
 * the inputs and the remaining dimensions the output
 */
void fgmm_regression_alloc_simple(struct fgmm_reg ** regression,
				 struct gmm * gmm,
				 int input_len);

/**
 * init a fgmm_reg structure 
 * @param gmm : the model to perform regression on 
 * @param input_len is the number of input dimensions
 * @param input_dim are the indexes of the input dimensions
 * @param output_len is the number of outputs 
 * @param output_dim is their indexes .. 
 */ 
void fgmm_regression_alloc(struct fgmm_reg ** regression,
			  struct gmm * gmm,
			  int input_len, int * input_dim,
			  int output_len, int * output_dim);

/**
 * free all the memory allocated for the regression structure 
 */
void fgmm_regression_free(struct fgmm_reg ** regression);

/**
 * does all intermediate computation, call this 
 * between alloc and regression. Call this also if 
 * the gmm relative to the regression has changed
 *
 * this caches the inverse covariance matrices for instance
 */
void fgmm_regression_init(struct fgmm_reg * reg);

/**
 * does the regression 
 */
void fgmm_regression(struct fgmm_reg * reg, const float * inputs, 
		     float * outputs, float * covar);


/**
 * Conditional sampling
 *
 * draw a sample in the output subspace, given the input point 
 * out ~ p(x | input) 
 */
void fgmm_regression_sampling(struct fgmm_reg * reg,const float * inputs,
			      float * output);

/**
 * incremental update, update the model with a new datapoint
 * 
 * Highly experimental .. 
 */
void fgmm_update(struct gmm * gmm, const float * data_point);

