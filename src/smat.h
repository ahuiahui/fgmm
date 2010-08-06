/** 
    Awesome fast positive definite matrix computations .. 
    
    florent.dhalluin@epfl.ch */
    
/**
 * stores Symetric matrices in a efficient way : 
 *   data is stored row first order such as for a 3x3 symetrical matrix 
 *     1 2 3
 * m = 2 4 5    ->  m->_ = [ 1 2 3 4 5 6] 
 *     3 5 6
 */

struct smat {
  float * _; /* data is actually stored here */
  int dim;   /* dimensionnality of the matrix */
  int _size; /* actual size of the data pointer dim * (dim+1) /2 */
};



/**
 * allocate memory for smat if *mat == NULL 
 * and zero the matrix in all cases 
 *
 * Use this function to allocate memory for smat structures
 */
void smat_zero(struct smat ** mat,int dim);


/**
 *  free the memory used by the matrix 
 *  (allocated by smat_zero) 
 *  *mat is set to NULL after that
 */
void smat_free(struct smat ** mat);

/**
 *  Matrix x Vector multiplication 
 *
 *  out = m * v 
 */

// fixme sous windows :: __inline__ 
inline void smat_multv(const struct smat* m, const float * v,float * out);

/**
 *  _in place _ Matrix x float  multiplication 
 */
inline void smat_multf(struct smat* m,const float *f);

/**
 * get the value at row col if the matrix were stored as 
 * square matrix 
 */
float smat_get_value(struct smat * mat,int row,int col);

/**
 * return the symetrical matrix considering only 
 * a given subset of dimensions 
 */
void smat_get_submatrix(struct smat * mat , struct smat * res,
			int n_dims, int * dims);

/* fill the mat with identity  /!\ you must first 
 allocate memory (wigh smat_zero for ex .. ) */
void smat_identity(struct smat * mat);

/** 
 * add value on the diagonal .. (e.g. add diag noise on 
 * a gaussian .. ) 
 */
void smat_add_diagonal(struct smat * mat , float value );


/* print matrix to screen (for debug purposes ) */ 
void smat_pmatrix(const struct smat* mat);

/* transform a symetric ordered matrix to a square one .. */ 
void smat_as_square(const struct smat* mat, float *  square);
void smat_from_square(struct smat * mat, const float * square);

/* Cholesky decomposition 
  
   out is a UPPER triang matrix such as  out^T * out = in 
   
   will fail if in is not SDP */

void smat_cholesky(const struct smat* in,struct smat* out);


/* L^T * L  for a triang SUP matrix, such as cholesky results .. */
void smat_ttmult(const struct smat* tri, struct smat* out);

/* 
   Using cholesky to invert systems : 

   foward resolve L*y = b  where L is * LOWER *  triangular 

   backward resolve U*y = b where U is upper triangular (cholesky results ) 


   to solve A*y = b  :

   smat_cholesky(A,U);
   smat_tforward(U,b,tmp);
   smat_tbackward(U,tmp,y);
*/ 

inline void smat_tforward(struct smat * lower, float * b, float * y) ;
inline void smat_tbackward(const struct smat * upper, float * b, float * y);

/**
 * computes sesquilinear form :
 *   (x - bias)^T Sigma^-1 (x-bias) 
 * given two vectors x and bias and the cholesky decomposition of Sigma 
 * ichol is actualy the cholesky decomposition of Sigma, where its value 
 * on the diagonal are inverted ... 
 */
inline float smat_sesq(struct smat * ichol,const float * bias,const float * x);

/**
 * compute the weighted covariance matrix of a given dataset
 * 
 * mainly used in the maximisation step of EM 
 * weight are the weights of each element of data
 * ndata : number of elements in data (data's size is ndata*cov->dim ) 
 * mean : the computed weighted mean (must be alloc'd before , its size is cov->dim ) 
 */

float smat_covariance(struct smat * cov, 
		      int ndata, 
		      const float * weight,
		      const float * data,
		      float * mean);



float smat_covariance_diag(struct smat * cov, 
			   int ndata, 
			   const float * weight,
			   const float * data,
			   float * mean);
