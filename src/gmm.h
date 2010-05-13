// #include "gaussian.h"

struct gaussian;

struct gmm {
  struct gaussian * gauss;
  int nstates;
  int dim;
};

/** 
 * draw one sample from the gmm 
 */
void fgmm_draw_sample(struct gmm *, float * out);

/**
 *  alloc all the memory needed by the model 
 *  
 *  all gaussians are init'd to zero mean, unity covariance 
 *  zero prior proba */
void fgmm_alloc(struct gmm *, int nstates, int dim);

/** free everything allocated by the above function 
 */
void fgmm_free(struct gmm *);

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


void fgmm_set_prior(struct gmm *,int state, float prior);
void fgmm_set_mean(struct gmm *,int state, const float * mean);
/**
 * Set the covariance of state # 
 *
 * Symetric matrix form : 
 * [[ 1 2 3 4 ]
 *  [ 2 5 6 7 ]
 *  [ 3 6 8 9 ]
 *  [ 4 7 9 10]] .. */
void fgmm_set_covar(struct gmm *,int state, float * covar);

/**
 * print the gmm parameters to screen 
 */
void fgmm_dump(struct gmm * gmm);

#define loglikelihood_eps 1e-4

/** EM algorithm */
int em( struct gmm * GMM,
	const float * data,
	int data_length, 
	float * end_loglikelihood,
	float likelihood_epsilon);

/**
 * return likelihood of point
 */
float fgmm_get_pdf( struct gmm * gmm,
		   float * point);


struct fgmm_reg;

void fgmm_regression_alloc_simple(struct fgmm_reg ** regression,
				 struct gmm * gmm,
				 int input_len);


void fgmm_regression_alloc(struct fgmm_reg ** regression,
			  struct gmm * gmm,
			  int input_len, int * input_dim,
			  int output_len, int * output_dim);

void fgmm_regression_free(struct fgmm_reg ** regression);

void fgmm_regression_init(struct fgmm_reg * reg);


void fgmm_regression(struct fgmm_reg * reg, float * inputs, float * outputs);
