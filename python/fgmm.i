%module fgmm
%{
  #define SWIG_FILE_WITH_INIT
  #include "fgmm++.hpp"
%}

%include "numpy.i"

%init %{
  import_array();
  %}

// This means that on the python side, the input array 
// must be explicitely using type np.dtype('f4')  -> 32 bits floats .. 

//%apply(float * IN_ARRAY1, int DIM1) {(float * data, int size)}


/* Typemap suite for (DATA_TYPE* IN_ARRAY2, DIM_TYPE DIM1, DIM_TYPE DIM2)
 */

// ---------------

class Gmm {
 public :
  Gmm(int states, int dim);
  ~Gmm();

  void Dump();

  void InitRegression(int ninputs);
  

  %apply(float * IN_ARRAY1, int DIM1){(float * input, int dimi)}
  %apply(float * ARGOUT_ARRAY1, int DIM1){(float * output, int dimo)}

  void DoRegression(float * input, int dimi, 
		    float * output, int dimo);

  void Draw(float * output, int dimo);
  
  %clear (float * input, int dimi , float * output, int dimo);
  
  %apply(float * IN_ARRAY2, int DIM1, int DIM2){(float * data, int size, int dim)}
  
  
  void init(float * data, int size,int dim);
  int Em(float * data, int size, int dim);

};
