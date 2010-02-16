#include "gaussian.h"

int em( struct gaussian * GMM,
	const float * data,
	int dim,
	int data_length, 
	int num_states,
	float * end_loglikelihood);
