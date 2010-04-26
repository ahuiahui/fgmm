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
void gmm_draw_sample(struct gmm *, float * out);

/**
 *  alloc all the memory needed by the model 
 *  
 *  all gaussians are init'd to zero mean, unity covariance 
 *  zero prior proba */
void gmm_alloc(struct gmm *, int nstates, int dim);

/** free everything allocated by the above function 
 */
void gmm_free(struct gmm *);

void gmm_set_prior(struct gmm *,int state, float prior);
void gmm_set_mean(struct gmm *,int state, float * mean);
/**
 * Set the covariance of state # 
 *
 * Symetric matrix form : 
 * [[ 1 2 3 4 ]
 *  [ 2 5 6 7 ]
 *  [ 3 6 8 9 ]
 *  [ 4 7 9 10]] .. */
void gmm_set_covar(struct gmm *,int state, float * covar);

/**
 * print the gmm parameters to screen 
 */
void gmm_dump(struct gmm * gmm);

/** EM algorithm */
int em( struct gmm * GMM,
	const float * data,
	int data_length, 
	float * end_loglikelihood);
