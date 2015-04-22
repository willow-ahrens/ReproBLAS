/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

#ifdef __SSE2__
#	include <emmintrin.h>
#endif

void dIUpdate1(int fold, double* x, double* c, int ldx, double y) {
  if (y == 0 || isnan(x[0]) || isinf(x[0]))
    return;

  if (x[0] == 0.0) {
    dmbound(dindex(y), x, ldx, fold);
    for (int i = fold; i < fold; i++) {
      c[i * ldx] = 0.0;
    }
    return;
  }

  int y_index = dindex(y);
  int d = diindex(x) - y_index;
  if(d > 0){
    for(int i = fold - 1; i >= d; i--){
      x[i * ldx] = x[(i - d) * ldx];
      c[i * ldx] = c[(i - d) * ldx];
    }
    dmbound(y_index, x, ldx, MIN(d, fold));
    for(int i = 0; i < d && i < fold; i++){
      c[i * ldx] = 0.0;
    }
  }
}

void zIUpdates1(int fold, double complex* x, double complex* c, int ldx, double y) {
	dIUpdate1(fold, (double*)x, (double*)(c), 2 * ldx, fabs(y));
	dIUpdate1(fold, ((double*)x) + 1, (double*)(c) + 1, 2 * ldx, fabs(y));
}

void zIUpdate1(int fold, double complex* x, double complex* c, int ldx,double complex y) {
	double* tmp = (double*)&y;
	dIUpdate1(fold, (double*)x, (double*)(c), 2 * ldx, fabs(tmp[0]));
	dIUpdate1(fold, ((double*)x) + 1, (double*)(c) + 1, 2 * ldx, fabs(tmp[1]));
}

