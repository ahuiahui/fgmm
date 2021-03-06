libfgmm : a fast(er) and lighter gaussian mixture models library. 

(c) Florent D'halluin  <florent.dhalluin@epfl.ch> , LASA laboratory. 

A pure C implementation of GMM, using a custom made structure to store 
covariance matrix. This structure exploits the symetric definite positiveness
of covariance matrix to be more computationnaly and memory efficient. 

This library is also intended to have no dependency. It actually 
haven't any, except a C compiler, and eventually python for the build
script. 

Compilation
===========

the build system here is waf ( http://code.google.com/p/waf/) an
awesome building system written in python ( http://python.org ) .  You
don't need to install it, since it is provided with the source code,
but you do need python, so go install it, and become a really cool
computer user.  If you don't want to, write your own build script,
this should not be a big deal anyway.

To compile run : 

> ./waf configure (-h shows you compilation options ) 
> ./waf 
> ./waf test
> ./waf install 

Documentation 
=============

this time using doxygen, a doxyfile is provided

> doxygen Doxyfile  

generates you html doc in the doc dir. 

see also the tests examples.  


Python wrappers 
===============

call ./waf configure with --python option. (you must have the python-dev and numpy-dev packages
somewhere ) 


Matlab/Octave wrappers 
======================

/!\ Highly untested .. 

For now I just wrapped the EM algorithm 

* requirements : 

tested with octave 3.2.3 and matlab  7.10.0
with octave you'll need the octave-dev packages 

* compile with 

> ./waf configure [--matlab | --octave]  
> ./waf

get the mex(glx) files in the _build_/default/matlab/ 
and put them wherever you want matlab/octave to find
them. (ie somewhere in the path : /usr/share/octave/packages/
the pwd .. ) 

* usage :  (see also matlab/test_em.m) 

(my data is  n datapointx x dim , you want a 8 state GMM) 

[priors,mu,Sigmas,likelihood] = fgmm_em(data, number_of_state, weights, 
			                initial_priors, mu , sigmas)  

if weights is empty ( [] ) assume all all points have 1. for weights.
initial parameters are optionals, in that case a random initialization 
is performed and may lead to strange results. 
