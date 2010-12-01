/************************************************************************/
/* This file is part of libfgmm.				        */
/* 								        */
/* libfgmm is free software: you can redistribute it and/or modify      */
/* it under the terms of the GNU Lesser General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.				        */
/* 								        */
/* libfgmm is distributed in the hope that it will be useful,	        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU Lesser General Public License for more details.			        */
/* 								        */
/* You should have received a copy of the GNU Lesser General Public License    */
/* along with libfgmm.  If not, see <http://www.gnu.org/licenses/>.     */
/* 								        */
/* Copyright 2010        LASA  - EPFL   http://lasa.epfl.ch             */
/*                                                                      */
/*       Florent D'halluin   <florent.dhalluin@epfl.ch>		        */
/************************************************************************/

/* FOR INTERNAL USE ONLY , defines structures for GMR */ 

#include "gaussian.h"
#include "fgmm.h"

struct fgmm_reg;

struct gaussian_reg {
  struct gaussian * gauss;    // the gaussian we will do regression 
  struct gaussian * subgauss; // input subgaussian Used to compute the weight of this
  struct fgmm_reg * reg;             // pointer to reg structure holding info on in/out dimensions ..
  _fgmm_real * reg_matrix; // store in->out A matrix 
};


struct fgmm_reg {
  struct gmm * model; // the actual GMM model
  int * input_dim;    // input dimension indexes
  int * output_dim;   // output dimension indexes
  int input_len;      
  int output_len;
  struct gaussian_reg * subgauss; // array of gaussian_reg
};

