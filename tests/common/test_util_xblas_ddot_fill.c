/*
Copyright (c) 2008-2009 The University of California Berkeley.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/
#include <stdlib.h>
#ifndef BLAS_FPU_H
#define BLAS_FPU_H

/* Contains code to set up the FPU control registers on x86
   systems.  The current double-double code requires that 
   all arithmetic is done in double precision (as opposed to
   double-extended).                                         */

#ifdef x86
#ifdef _WIN32

#include <float.h>
#define FPU_FIX_DECL unsigned int __old_cw, __new_cw;
#define FPU_FIX_START \
  __old_cw = _control87(0, 0);  \
  __new_cw = _control87(0x00010000, 0x00030000);
#define FPU_FIX_STOP \
  _control87(*_old_cw, 0xFFFFFFFF);
#else  /* _WIN32 */

#if HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#define FPU_FIX_DECL unsigned short __old_cw, __new_cw;
#define FPU_FIX_START \
  _FPU_GETCW(__old_cw); \
  __new_cw = (__old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE; \
  _FPU_SETCW(__new_cw); 
#define FPU_FIX_STOP \
  _FPU_SETCW(__old_cw);
#endif  /* else _WIN32 */

#else   /* x86 */
#define FPU_FIX_DECL
#define FPU_FIX_START
#define FPU_FIX_STOP
#endif  /* else x86 */

#endif /* BLAS_FPU_H */

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

static void BLAS_error(const char *rname, int iflag, int ival, char *form, ...)
/*
 * Argument
 * ========
 * rname     (input) routine name
 *
 * iflag     (input) a negative value indicates that parameter number -IFLAG
 *                   caused the error; a nonnegative value is an
 *                   implementation-specific error code.
 *
 * ival      (input) the value of parameter number -IFLAG.
 */
{
  {
    va_list argptr;
    va_start(argptr, form);
    fprintf(stderr, "Error #%d from routine %s:\n", iflag, rname);
    if (form)
      vfprintf(stderr, form, argptr);
    else if (iflag < 0)
      fprintf(stderr,
	      "  Parameter number %d to routine %s had the illegal value %d\n",
	      -iflag, rname, ival);
    else
      fprintf(stderr, "  Unknown error code %d from routine %s\n",
	      iflag, rname);
    exit(iflag);
  }
}

enum blas_prec_type {
            blas_prec_single     = 211,
            blas_prec_double     = 212,
            blas_prec_indigenous = 213,
            blas_prec_extra      = 214 };

enum blas_conj_type {
            blas_conj    = 191,
            blas_no_conj = 192 };

/* constants */
#define BITS_S  24
#define BITS_D  53
#define BITS_E  106

/* Split a double into 2 parts with at most 26 bits each. (2^27 + 1) */
#define split 	(134217729.0)

/* macros */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

static double power(int i1, int i2)
{
  int i, j;
  double r = 1.0;

  if (i2 < 0)
    j = -i2;
  else
    j = i2;
  for (i = 0; i < j; ++i)
    r *= i1;
  if (i2 < 0)
    r = 1. / r;

  return r;
}

static double xrand(int *is)
/*
 * XRAND returns a uniformly distributed pseudorandom number in (0, 1).
 * IS is the seed, and is changed with each call.  The period of this
 * linear congruential generator is 2^26, according to Knuth vol. 2.
 */
{
  double s1, s2, ret_val;
#define f7    78125.0           /* 5.d0 ** 7 */
#define r26   1.4901161193847656e-8     /* 2^(-26) */
#define r28   3.7252902984619141e-9     /* 2^(-28) */
#define t28   268435456.0       /* 2^28 */

  s1 = *is;
  s2 = fmod(f7 * s1, t28);
  ret_val = (s1 + r26 * s2) * r28;
  *is = s2;

  return ret_val;
}

static int FixedBits(double r_true_l, double r_true_t)
/*
 * Purpose
 * =======
 *
 * Compute the number of fixed bits in r_true.
 * (r_true_l, r_true_t) is double-double representation of r_true.
 *
 */
{

  int b;                        /* Number of fixed bits in r_true */
  int i, k;
  double tmp_l, tmp_t;
  double res[5], t, temp;

  b = k = 0;
  res[0] = r_true_l;
  while (res[k] != 0.0) {       /* Each time cancel 53 bits */
    tmp_l = r_true_l;
    tmp_t = r_true_t;
    for (i = 0; i <= k; ++i) {  /* tmp = tmp - res[i] */
      t = -res[i];
      {
        /* Compute double-double = double-double + double. */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = tmp_l + t;
        e = t1 - tmp_l;
        t2 = ((t - e) + (tmp_l - (t1 - e))) + tmp_t;

        /* The result is t1 + t2, after normalization. */
        tmp_l = t1 + t2;
        tmp_t = t2 - (tmp_l - t1);
      }
    }
    b += BITS_D;
    ++k;
    res[k] = tmp_l;
  }

  if (k > 0) {
    int e;
    double f;

    /* Use Kahan's trick for the last nonzero residual. */
    b -= BITS_D;
    --k;
    t = fabs(res[k]);
    f = frexp(t, &e);           /* 1/2 <= f < 1 */
    /* g = f * 2;              1 <= g < 2 */
    for (i = 1;; ++i) {         /* Compute number of bits in g */
      t = f * power(2, i);      /* Shift left i bits */
      temp = floor(t);
      t -= temp;
      if (t == 0.0)
        break;
      ++b;
    }
  }
  return b;
}                               /* FixedBits */

static void BLAS_ddot_x(enum blas_conj_type conj, int n, double alpha,
                        const double *x, int incx, double beta,
                        const double *y, int incy,
                        double *r, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *     r <- beta * r + alpha * SUM_{i=0, n-1} x[i] * y[i].
 * 
 * Arguments
 * =========
 *  
 * conj   (input) enum blas_conj_type
 *        When x and y are complex vectors, specifies whether vector
 *        components x[i] are used unconjugated or conjugated. 
 * 
 * n      (input) int
 *        The length of vectors x and y.
 * 
 * alpha  (input) double
 * 
 * x      (input) const double*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) double
 *
 * y      (input) const double*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) double*
 * 
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  static const char routine_name[] = "BLAS_ddot_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      double *r_i = r;
      const double *x_i = x;
      const double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      double x_ii;
      double y_ii;
      double r_v;
      double prod;
      double sum;
      double tmp1;
      double tmp2;


      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if ((beta_i == 1.0) && (n == 0 || (alpha_i == 0.0)))
	return;



      r_v = r_i[0];
      sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	prod = x_ii * y_ii;	/* prod = x[i]*y[i] */
	sum = sum + prod;	/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      tmp1 = sum * alpha_i;	/* tmp1 = sum*alpha */
      tmp2 = r_v * beta_i;	/* tmp2 = r*beta */
      tmp1 = tmp1 + tmp2;	/* tmp1 = tmp1+tmp2 */
      *r = tmp1;		/* r = tmp1 */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      double *r_i = r;
      const double *x_i = x;
      const double *y_i = y;
      double alpha_i = alpha;
      double beta_i = beta;
      double x_ii;
      double y_ii;
      double r_v;
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1, tail_tmp1;
      double head_tmp2, tail_tmp2;
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if ((beta_i == 1.0) && (n == 0 || (alpha_i == 0.0)))
	return;

      FPU_FIX_START;

      r_v = r_i[0];
      head_sum = tail_sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x_ii * split;
	  a1 = con - x_ii;
	  a1 = con - a1;
	  a2 = x_ii - a1;
	  con = y_ii * split;
	  b1 = con - y_ii;
	  b1 = con - b1;
	  b2 = y_ii - b1;

	  head_prod = x_ii * y_ii;
	  tail_prod = (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	}			/* prod = x[i]*y[i] */
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_sum + head_prod;
	  bv = s1 - head_sum;
	  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_sum + tail_prod;
	  bv = t1 - tail_sum;
	  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_sum = t1 + t2;
	  tail_sum = t2 - (head_sum - t1);
	}			/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      {
	/* Compute double-double = double-double * double. */
	double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	con = head_sum * split;
	a11 = con - head_sum;
	a11 = con - a11;
	a21 = head_sum - a11;
	con = alpha_i * split;
	b1 = con - alpha_i;
	b1 = con - b1;
	b2 = alpha_i - b1;

	c11 = head_sum * alpha_i;
	c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	c2 = tail_sum * alpha_i;
	t1 = c11 + c2;
	t2 = (c2 - (t1 - c11)) + c21;

	head_tmp1 = t1 + t2;
	tail_tmp1 = t2 - (head_tmp1 - t1);
      }				/* tmp1 = sum*alpha */
      {
	/* Compute double_double = double * double. */
	double a1, a2, b1, b2, con;

	con = r_v * split;
	a1 = con - r_v;
	a1 = con - a1;
	a2 = r_v - a1;
	con = beta_i * split;
	b1 = con - beta_i;
	b1 = con - b1;
	b2 = beta_i - b1;

	head_tmp2 = r_v * beta_i;
	tail_tmp2 = (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
      }				/* tmp2 = r*beta */
      {
	/* Compute double-double = double-double + double-double. */
	double bv;
	double s1, s2, t1, t2;

	/* Add two hi words. */
	s1 = head_tmp1 + head_tmp2;
	bv = s1 - head_tmp1;
	s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

	/* Add two lo words. */
	t1 = tail_tmp1 + tail_tmp2;
	bv = t1 - tail_tmp1;
	t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

	s2 += t1;

	/* Renormalize (s1, s2)  to  (t1, s2) */
	t1 = s1 + s2;
	s2 = s2 - (t1 - s1);

	t2 += s2;

	/* Renormalize (t1, t2)  */
	head_tmp1 = t1 + t2;
	tail_tmp1 = t2 - (head_tmp1 - t1);
      }				/* tmp1 = tmp1+tmp2 */
      *r = head_tmp1;		/* r = tmp1 */

      FPU_FIX_STOP;
    }
    break;
  }
}

static double rand_half_1(int, int *);
static double ulp(double);
static void
r_truth(enum blas_conj_type, int, double, const double *, int,
        double, const double *, int, double *, double *, double *);
static void
gen_y_to_cancel(int, int, enum blas_conj_type, double, double *, double *);
static void
gen_r_to_cancel(int, enum blas_conj_type, double, double,
                double *, double *, double *, int *);

void
util_xblas_ddot_fill(int n, int n_fix2, int n_mix, int norm,
                  char charconj,
                  double *alpha, int alpha_flag, double *beta, int beta_flag,
                  double *x, double *y, int *seed,
                  double *r, double *r_true_l, double *r_true_t)
/*
 * Purpose
 * =======
 *
 * This routine generates the test vectors X and Y for C_DDOT.
 *
 * Arguments
 * =========
 * 
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int* 
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int B, frees, y_free, i, k, m, s;
  double a, rtmps;
  double rtmpd, eps, eps_out;
  double f;

  enum blas_conj_type conj = blas_no_conj;
  if(charconj == 'c' || charconj == 'C'){
    conj = blas_conj;
  }

  if (alpha_flag == 0)
    *alpha = xrand(seed);
  if (beta_flag == 0)
    *beta = xrand(seed);

  y_free = n - n_fix2;
  k = n_fix2;
  eps_out = power(2, -BITS_D);

  /*
   * Compute the number of bits in the prefix sum:
   *     alpha * SUM_{i=0,n_fix2-1}(x[i] * y[i])
   */
  *r = 0.0;
  r_truth(conj, n_fix2, *alpha, x, 1, 0.0, y, 1, r, r_true_l, r_true_t);
  B = FixedBits(*r_true_l, *r_true_t);

  /* Pick r at random */
  *r = xrand(seed);

  /* Pick the free X(i)'s at random. */
  for (i = n_fix2 + n_mix; i < n; ++i)
    x[i] = xrand(seed);

  if (alpha_flag == 1 && *alpha == 0.0) {
    /* Pick the free Y(i)'s at random. */
    for (i = n_fix2; i < n; ++i)
      y[i] = xrand(seed);
    /* Compute r_truth in double-double */
    r_truth(conj, n, *alpha, x, 1, *beta, y, 1, r, r_true_l, r_true_t);
    return;
  }

  if (beta_flag == 1 && *beta == 0.0) {
    if (B == 0) {               /* Assume alpha is not very big . */
      switch (y_free) {
      case 0:
        break;
      case 1:
        y[n_fix2] = xrand(seed);
        break;
      case 2:
        /*
         * Make SUM_{i=0,1}(x[k+i] * y[k+i]) small ...
         * ... = alpha * eps^2
         */
        if (n_mix == 0) {       /* Both x[k] and x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = eps_out^2, small */
          a = rand_half_1(26, seed);    /* strictly < 1 */
          x[k] = a;
          y[k] = a;
          x[k + 1] = a + eps_out;       /* exact */
          y[k + 1] = -a + eps_out;      /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed, x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = eps_out^2 */
          a = x[k];
          y[k] = a;
          eps = ulp(a);
          x[k + 1] = a + eps;   /* exact */
          y[k + 1] = -a + eps;  /* exact */
        } else {                /* Both x[k] and x[k+1] fixed; cancel 53
                                 * bits. */
          y[k] = xrand(seed);
          gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
        }
        break;
      case 3:
        /*
         * Make SUM_{i=0,2}(x[k+i] * y[k+i]) small
         * ... x[k]*y[k] = -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small
         */
        y[k] = -x[k + 2];
        y[k + 2] = x[k];
        f = x[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (x[k + 1] == 0.)
          y[k + 1] = 0.;
        else
          y[k + 1] = f * power(2, -s) / x[k + 1];
        break;
      case 4:
        /*
         * Make SUM_{i=0,3}(x[k+i] * y[k+i]) small
         * ... x[k]*y[k] = -x[k+3]*y[k+3] exact, x[k+1]*y[k+1] small,
         *     x[k+2]*y[k+2] = 0.
         */
        y[k] = -x[k + 3];
        y[k + 3] = x[k];
        f = x[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (x[k + 1] == 0.)
          y[k + 1] = 0.;
        else
          y[k + 1] = f * power(2, -s) / x[k + 1];
        y[k + 2] = 0.0;
        break;
      default:                 /* y_free >= 5 */
        /*
         * Make SUM_{i=0,n-1}(x[k+i] * y[k+i]) small
         * ... Cancel >= 106 bits.
         */
        a = xrand(seed);        /* Flip a coin */
        if (a < 0.5) {
          /*
           * Use last 2 to cancel bits, and leading ones to add bits.
           */
          y[k] = xrand(seed);
          rtmpd = *alpha * x[k] * y[k];
          s = 30;
          for (i = k + 1; i < n - 2; ++i) {
            rtmpd *= power(2, -s);
            if (*alpha * x[i] == 0.)
              y[i] = 0.;
            else
              y[i] = rtmpd / (*alpha * x[i]);
          }
          gen_y_to_cancel(n - 2, n, conj, *alpha, x, y);

        } else {                /* Use the scheme similar to case 3. */
          m = (y_free - 1) / 2;
          rtmpd = 0.0;
          for (i = 0; i < m; ++i) {
            y[k + i] = -x[n - i - 1];
            rtmpd = MAX(rtmpd, fabs(x[k + i] * y[k + i]));
          }
          for (i = 0; i < m; ++i)
            y[n - i - 1] = x[k + i];
          s = 100;              /* Should be random between 40 and 100. */
          if (x[k + m] == 0.)
            y[k + m] = 0.;
          else
            y[k + m] = rtmpd * power(2, -s) / x[k + m];
          if (y_free % 2 == 0)
            y[k + m + 1] = 0.0;
        }
        break;
      }                         /* end switch */
    } else {                    /* B > 0 */
      if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
        gen_y_to_cancel(k, n, conj, *alpha, x, y);
      } else {                  /* At least y[n-1] is free. */
        if (y_free == 1) {
          /* Cancel min(B,53) bits. */
          gen_y_to_cancel(k, n, conj, *alpha, x, y);
        } else {                /* >= 2 frees. */
          /*
           * There are 2 possibilities:
           * (1) use 1 to add bits, and y_free-1 to cancel 53*(y_free-1) bits
           * (2) use all to cancel min(B, 53*y_free) bits
           * Goal is to maximize the # of bits cancelled. By equating (1)
           * and (2), we find the crossover point is y_free = B/53 + 1.
           */
          if (y_free > B / (double) BITS_D + 1) {       /* Use scheme (1) */
            rtmps = 0.0;
            BLAS_ddot_x(conj, k, *alpha, x, 1, 0.0, y, 1,
                        &rtmps, blas_prec_extra);
            rtmpd = rtmps;
            s = 100;            /* Should be random between 40 and 100. */
            if (*alpha * x[k] == 0.)
              y[k] = 0.;
            else
              y[k] = rtmpd * power(2, -s) / (*alpha * x[k]);
            gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
          } else {              /* Use scheme (2) */
            gen_y_to_cancel(k, n, conj, *alpha, x, y);
          }
        }
      }                         /* end else B < 106 */
    }                           /* end else B > 0 */

    /* Compute r_truth in double-double */
    r_truth(conj, n, *alpha, x, 1, *beta, y, 1, r, r_true_l, r_true_t);
    return;
  }

  /* if beta == 0 */
  /* Now, beta is non-zero. */
  if (B == 0) {
    switch (y_free) {
    case 0:
      break;
    case 1:
      /* Make alpha*x[k]*y[k] + beta*r small. */
      /* Count number of frees in alpha, x[k], and beta. */
      frees = 0;
      if (alpha_flag == 0)
        ++frees;
      if (beta_flag == 0)
        ++frees;
      if (n_mix == 0)
        ++frees;
      if (frees >= 2) {
        /* alpha*x[k]*y[k] + beta*r = -alpha * eps_out^2 */
        a = rand_half_1(26, seed);      /* strictly < 1, only leading 26 bits */
        *r = -a * a;            /* exact */
        if (alpha_flag == 1) {  /* alpha fixed */
          *beta = *alpha;
          x[k] = a + eps_out;   /* exact */
          y[k] = a - eps_out;   /* exact */
        } else if (beta_flag == 1) {    /* beta fixed */
          *alpha = *beta;
          x[k] = a + eps_out;   /* exact */
          y[k] = a - eps_out;   /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed */
          *beta = x[k];
          *alpha = a + eps_out; /* exact */
          y[k] = a - eps_out;   /* exact */
        } else {                /* alpha, x[k], and beta all free */
          *beta = *alpha;
          x[k] = a + eps_out;   /* exact */
          y[k] = a - eps_out;   /* exact */
        }
      } else {                  /* Cancel 53 bits. */
        y[k] = xrand(seed);
        gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
      }
      break;
    case 2:
      /* Make SUM_{i=0,1}(alpha * x[k+i] * y[k+i]) + beta*r small. */
      /* Count number of frees in alpha, x[k], and beta. */
      frees = 0;
      if (alpha_flag == 0)
        ++frees;
      if (beta_flag == 0)
        ++frees;
      if (n_mix == 0)
        ++frees;
      if (frees > 0) {
        /*
         * Make alpha*x[k]*y[k] = -beta*r exact, alpha*x[k+1]*y[k+1] small.
         */
        y[k] = -1.0;
        if (alpha_flag == 0) {  /* alpha free */
          *alpha = *beta;
          *r = x[k];
        } else if (beta_flag == 0) {    /* beta free */
          *beta = *alpha;
          *r = x[k];
        } else {                /* x[k] free */
          x[k] = *beta;
          *r = *alpha;
        }
        f = *alpha * x[k] * y[k];
        s = 100;                /* Should be random between 40 and 100. */
        if (*alpha * x[k + 1] == 0.)
          y[k + 1] = 0.;
        else
          y[k + 1] = f * power(2, -s) / (*alpha * x[k + 1]);
      } else {                  /* Cancel 53 bits. */
        y[k] = xrand(seed);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
      }
      break;
    default:                   /* Actual frees >= 4 */
      /*
       * Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
       * ... Cancel >= 106 bits.
       */
      a = xrand(seed);          /* Flip a coin */
      if (a < 0.5) {
        /*
         * Use last 2 ( Y(n) and r) to cancel bits, and
         * leading ones to add bits.
         */
        y[k] = xrand(seed);
        rtmpd = *alpha * x[k] * y[k];
        s = 30;
        for (i = k + 1; i < n - 1; ++i) {
          rtmpd *= power(2, -s);
          if (*alpha * x[i] == 0.)
            y[i] = 0.;
          else
            y[i] = rtmpd / (*alpha * x[i]);
        }
        gen_y_to_cancel(n - 1, n, conj, *alpha, x, y);
        gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
      } else {                  /* Use the scheme similar to case 2. */
        /*
         * Make the middle one tiny, and the rest cancel each other
         * from both ends.
         */
        *r = 0.0;
        m = (y_free - 1) / 2;
        rtmpd = 0.0;
        for (i = 0; i < m; ++i) {
          y[k + i] = -x[n - i - 1];
          rtmpd = MAX(rtmpd, fabs(x[k + i] * y[k + i]));
        }
        for (i = 0; i < m; ++i)
          y[n - i - 1] = x[k + i];
        s = 100;                /* Should be random between 40 and 100. */
        if (x[k + m] == 0.)
          y[k + m] = 0.;
        else
          y[k + m] = rtmpd * power(2, -s) / x[k + m];
        if (y_free % 2 == 0)
          y[k + m + 1] = 0.0;
      }
      break;
    }                           /* end switch */
  } else {                      /* B > 0 */
    if (B >= BITS_E) {          /* Choose Y(i)'s and r to cancel. */
      gen_y_to_cancel(k, n, conj, *alpha, x, y);
      gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
    } else {                    /* >= 2 frees. Use y[k] to add bits. */
      frees = y_free + 1;
      /*
       * There are 2 possibilities:
       * (1) use 1 to add bits, and y_free-1 to cancel 53*(y_free-1) bits
       * (2) use all to cancel min(B, 53*y_free) bits
       * Goal is to maximize the # of bits cancelled. By equating (1)
       * and (2), we find the crossover point is y_free = B/53 + 1.
       */
      if (frees > B / (double) BITS_D + 1) {    /* Use scheme (1) */
        rtmps = 0.0;
        BLAS_ddot_x(conj, k, *alpha, x, 1, 0.0, y, 1,
                    &rtmps, blas_prec_extra);
        s = 100;                /* Should be random between 40 and 100. */
        if (*alpha * x[k] == 0.)
          y[k] = 0.;
        else
          y[k] = rtmps * power(2, -s) / (*alpha * x[k]);
        gen_y_to_cancel(k + 1, n, conj, *alpha, x, y);
      } else {                  /* Use scheme (2) */
        gen_y_to_cancel(k, n, conj, *alpha, x, y);
      }
      gen_r_to_cancel(n, conj, *alpha, *beta, x, y, r, seed);
    }
  }

  /* Compute r_truth in double-double */
  r_truth(conj, n, *alpha, x, 1, *beta, y, 1, r, r_true_l, r_true_t);
}                               /* testgen_BLAS_ddot */


static double ulp(double a)
/*
 * Purpose 
 * =======
 * 
 * Compute the unit last place of a double precision number.
 */
{
  int e;
  frexp(a, &e);
  return power(2, e - BITS_D);
}


static double rand_half_1(int l_bits, int *seed)
/*
 * Purpose =======
 * 
 * Generate random number in the interval [0.5, 1). l_bits specifies that only
 * the leading l_bits are nonzero.
 * 
 */
{
  double a = xrand(seed);       /* [0,1] */
  a /= 2.;
  a += 0.5;
  if (l_bits < BITS_D) {
    double s = power(2, l_bits);
    double t = a / s;           /* shift right l_bits */
    t = (t + a) - a;            /* cancel trailing bits */
    a = t * s;                  /* shift back */
  }
  return a;
}

static void
gen_y_to_cancel(int k, int n, enum blas_conj_type conj,
                double alpha, double *x, double *y)
/*
 * Purpose =======
 * 
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 * 
 */
{
  int i;
  double rtmp;

  for (i = k; i < n; ++i) {
    rtmp = 0.0;
    BLAS_ddot_x(conj, i, alpha, x, 1, 0.0, y, 1, &rtmp, blas_prec_extra);
    if (alpha * x[i] == 0.)
      y[i] = 0.;
    else
      y[i] = -rtmp / (alpha * x[i]);
  }
}


static void
gen_r_to_cancel(int n, enum blas_conj_type conj,
                double alpha, double beta,
                double *x, double *y, double *r, int *seed)
/*
 * Purpose =======
 * 
 * Generate r to cancel as much as possible.
 * 
 */
{
  double rtmp;

  if (beta == 0.0)
    *r = xrand(seed);
  else {
    rtmp = 0.0;
    BLAS_ddot_x(conj, n, alpha, x, 1, 0.0, y, 1, &rtmp, blas_prec_extra);
    *r = -rtmp / beta;
  }
}


static void r_truth(enum blas_conj_type conj, int n, double alpha, const double *x, int incx, double beta, const double *y, int incy, double *r,        /* input */
                    double *r_true_l, double *r_true_t)
{
  int i, ix = 0, iy = 0;
  double *r_i = r;
  const double *x_i = x;
  const double *y_i = y;
  double alpha_i = alpha;
  double beta_i = beta;
  double x_ii;
  double y_ii;
  double r_v;
  double prod_l, prod_t;
  double sum_l, sum_t;
  double tmp1_l, tmp1_t;
  double tmp2_l, tmp2_t;

  /* Immediate return */
  if (n < 0) {
    *r_true_l = *r_true_t = 0.0;
    return;
  }

  r_v = r_i[0];
  sum_l = sum_t = 0.0;          /* sum = 0 */


  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incx;
  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = y_i[iy];
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = x_ii * split;
      a1 = con - x_ii;
      a1 = con - a1;
      a2 = x_ii - a1;
      con = y_ii * split;
      b1 = con - y_ii;
      b1 = con - b1;
      b2 = y_ii - b1;

      prod_l = x_ii * y_ii;
      prod_t = (((a1 * b1 - prod_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }                           /* prod = x[i]*y[i] */
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = sum_l + prod_l;
      e = t1 - sum_l;
      t2 = ((prod_l - e) + (sum_l - (t1 - e))) + sum_t + prod_t;

      /* The result is t1 + t2, after normalization. */
      sum_l = t1 + t2;
      sum_t = t2 - (sum_l - t1);
    }                           /* sum = sum+prod */
    ix += incx;
    iy += incy;
  }                             /* endfor */
  {
    /* Compute double-double = double-double * double. */
    double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

    con = sum_l * split;
    a11 = con - sum_l;
    a11 = con - a11;
    a21 = sum_l - a11;
    con = alpha_i * split;
    b1 = con - alpha_i;
    b1 = con - b1;
    b2 = alpha_i - b1;

    c11 = sum_l * alpha_i;
    c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

    c2 = sum_t * alpha_i;
    t1 = c11 + c2;
    e = t1 - c11;
    t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

    tmp1_l = t1 + t2;
    tmp1_t = t2 - (tmp1_l - t1);
  }                             /* tmp1 = sum*alpha */
  {
    /* Compute double_double = double * double. */
    double a1, a2, b1, b2, con;

    con = r_v * split;
    a1 = con - r_v;
    a1 = con - a1;
    a2 = r_v - a1;
    con = beta_i * split;
    b1 = con - beta_i;
    b1 = con - b1;
    b2 = beta_i - b1;

    tmp2_l = r_v * beta_i;
    tmp2_t = (((a1 * b1 - tmp2_l) + a1 * b2) + a2 * b1) + a2 * b2;
  }                             /* tmp2 = r*beta */
  {
    /* Compute double-double = double-double + double-double. */
    double e, t1, t2;

    /* Knuth trick. */
    t1 = tmp1_l + tmp2_l;
    e = t1 - tmp1_l;
    t2 = ((tmp2_l - e) + (tmp1_l - (t1 - e))) + tmp1_t + tmp2_t;

    /* The result is t1 + t2, after normalization. */
    tmp1_l = t1 + t2;
    tmp1_t = t2 - (tmp1_l - t1);
  }                             /* tmp1 = tmp1+tmp2 */

  /* Return r_truth */
  *r_true_l = tmp1_l;
  *r_true_t = tmp1_t;
}                               /* end r_truth */
