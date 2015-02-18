/*
 *  Created   February 2015 Peter Ahrens
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"

void rdgemv(const rblas_order_t order,
            const rblas_transpose_t TransA, const int M, const int N,
            const double alpha, const double *A, const int lda,
            const double *X, const int incX, const double beta,
            double *Y, const int incY){
  int i;
  switch(TransA){
    case :
  }
  for(i = 0; i < 

