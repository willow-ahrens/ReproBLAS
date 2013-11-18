/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rblas1.h"
#include "IndexedFP/dIndexed.h"
#include <omp.h>

#define MACHEPS DBL_EPSILON
#define PREC    53
#define OVERFLOW_THRES DBL_MAX_EXP
#define CHECK_NAN_INF

I_double ddotI(
	int N,
	double* x, int incx,
	double* y, int incy
) {
	I_double ret;
	dISetZero(ret);
#pragma omp parallel
	{
		int nthreads = omp_get_num_threads();
		int me = omp_get_thread_num();
		int q = N / nthreads;
		int r = N % nthreads;
		int B = (me < r) ? q + 1 : q;
		int start = (me < r) ? (q + 1) * me : (q+1) * r + (me - r) * q;
		Idouble tmp;
		dISetZero(tmp);

		ddotI1(B, x + start * incx, incx, y + start * incy, incy,
			DI_DEFAULT_FOLD, 0, tmp);

		#pragma omp critical
		{
			dIAdd(ret, tmp);
		}
	}

	return ret;
}

double rddot(
	int N,
	double* x, int incx,
	double* y, int incy
) {

	Idouble dot;

	dot = ddotI(N, x, incx, y, incy);

	return Iconv2d(dot);
}
