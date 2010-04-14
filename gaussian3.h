
struct gaussian3d {
  float prior;
  float mean[3];
  float covar[6];  /* speaking of 3x3 **symetric** matrix .. */
  float icovar[6];
  float nfactor;
  /*  [ 0 3 5 
        - 1 4 
	- - 2 ] */
};   

inline float gaussian_pdf(struct gaussian3d* g, const float* value);
void init_random(struct gaussian3d* g);
void invert_covar(struct gaussian3d* g);
void dump(struct gaussian3d* g);
