/** 
    Awesome fast positive definite matrix computations .. 
    
    florent.dhalluin@epfl.ch */
    

struct smat {
  float * _;
  int dim;
  int _size;
};

inline void smat_multv(const struct smat* m, const float * v,float * out);
inline void smat_multf(struct smat* m,const float *f);

/* allocate memory for smat if mat == NULL 
   and zero the matrix in all cases */
void smat_zero(struct smat ** mat,int dim);

/* free the memory used by the matrix 
   (allocated by smat_zero) */
void smat_free(struct smat ** mat);

/* fill the mat with identity  /!\ you must first 
 allocate memory (wigh smat_zero for ex .. ) */
void smat_identity(struct smat * mat);


/* print matrix to screen (for debug purposes ) */ 
void smat_pmatrix(const struct smat* mat);


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



float smat_covariance(struct smat * cov, 
		      int ndata, 
		      const float * weight,
		      const float * data,
		      float * mean);



