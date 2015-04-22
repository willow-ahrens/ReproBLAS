/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "indexed.h"
#include "../Common/Common.h"

void sIUpdate1(int fold, float y, float* x, float* c, int ldx) {
  if (y == 0 || isnan(x[0]) || isinf(x[0]))
    return;

  if (x[0] == 0.0) {
    smbound(sindex(y), x, ldx, fold);
    for (int i = fold; i < fold; i++) {
      c[i * ldx] = 0.0;
    }
    return;
  }

  int y_index = sindex(y);
  int d = siindex(x) - y_index;
  if(d > 0){
    for(int i = fold - 1; i >= d; i--){
      x[i * ldx] = x[(i - d) * ldx];
      c[i * ldx] = c[(i - d) * ldx];
    }
    smbound(y_index, x, ldx, MIN(d, fold));
    for(int i = 0; i < d && i < fold; i++){
      c[i * ldx] = 0.0;
    }
  }
}

void cIUpdates1(int K, float complex* X, float* C, int INC, float Y) {
	Y = fabs(Y);
	sIUpdate1(K,Y,(float*)X  , C  , 2*INC);
	sIUpdate1(K,Y,(float*)X+1, C+1, 2*INC);
}

void cIUpdate1(int K, float complex* X, float* C,int INC, float complex Y) {
	float* tmp = (float*)&Y;
	sIUpdate1(K, fabs(tmp[0]), (float*)X  , C  , 2*INC);
	sIUpdate1(K, fabs(tmp[1]), (float*)X+1, C+1, 2*INC);
}

