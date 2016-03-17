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
#include <math.h>
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


static void ddadd(double dda_l, double dda_t, double ddb_l, double ddb_t,
                  double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 *
 * This subroutine computes ddc = dda + ddb.
 *
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
  double e, t1, t2;

  /* Compute dda + ddb using Knuth's trick. */
  t1 = dda_l + ddb_l;
  e = t1 - dda_l;
  t2 = ((ddb_l - e) + (dda_l - (t1 - e))) + dda_t + ddb_t;

  /* The result is t1 + t2, after normalization. */
  *ddc_l = t1 + t2;
  *ddc_t = t2 - (*ddc_l - t1);

}                               /* end ddadd */

static void ddmuld(double dda_l, double dda_t, double db,
                   double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 *
 * This routine multiplies the DD number DDA by the DP number DB to yield
 * the DD product DDC.
 *
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
  double a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2;

  /* This splits dda(1) and db into high-order and low-order words. */
  cona = dda_l * split;
  conb = db * split;
  a1 = cona - (cona - dda_l);
  b1 = conb - (conb - db);
  a2 = dda_l - a1;
  b2 = db - b1;

  /* Multilply dda(1) * db using Dekker's method. */
  c11 = dda_l * db;
  c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

  /* Compute dda(2) * db (only high-order word is needed). */
  c2 = dda_t * db;

  /* Compute (c11, c21) + c2 using Knuth's trick. */
  t1 = c11 + c2;
  e = t1 - c11;
  t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

  /* The result is t1 + t2, after normalization. */
  *ddc_l = t1 + t2;
  *ddc_t = t2 - (*ddc_l - t1);

}                               /* end ddmuld */

static void dddiv(double dda_l, double dda_t,
                  double ddb_l, double ddb_t, double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 * 
 * This divides the DD number DDA by the DD number DDB to yield the DD
 * quotient DDC.
 * 
 * Taken from D. H. Bailey's ddfun90.f.
 *
 */
  double a1, a2, b1, b2, cona, conb, c11, c2, c21, e,
    s1, s2, t1, t2, t11, t12, t21, t22;

  /* Compute a DP approximation to the quotient. */
  s1 = dda_l / ddb_l;

  /* This splits s1 and ddb(1) into high-order and low-order words. */
  cona = s1 * split;
  conb = ddb_l * split;
  a1 = cona - (cona - s1);
  b1 = conb - (conb - ddb_l);
  a2 = s1 - a1;
  b2 = ddb_l - b1;

  /* Multiply s1 * ddb(1) using Dekker's method. */
  c11 = s1 * ddb_l;
  c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

  /* Compute s1 * ddb(2) (only high-order word is needed). */
  c2 = s1 * ddb_t;

  /* Compute (c11, c21) + c2 using Knuth's trick. */
  t1 = c11 + c2;
  e = t1 - c11;
  t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

  /* The result is t1 + t2, after normalization. */
  t12 = t1 + t2;
  t22 = t2 - (t12 - t1);

  /* Compute dda - (t12, t22) using Knuth's trick. */
  t11 = dda_l - t12;
  e = t11 - dda_l;
  t21 = ((-t12 - e) + (dda_l - (t11 - e))) + dda_t - t22;

  /* Compute high-order word of (t11, t21) and divide by ddb(1). */
  s2 = (t11 + t21) / ddb_l;

  /* The result is s1 + s2, after normalization. */
  *ddc_l = s1 + s2;
  *ddc_t = s2 - (*ddc_l - s1);

}                               /* end dddiv */


static void z_ddmuld(double *dda_l, double *dda_t, double *db,
                     double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 * 
 * This multiplies the complex DD number DDA by the complex D number DB to
 * yield the complex DD product DDC.
 * 
 */
  double t_l, t_t, t1_l, t1_t, t2_l, t2_t;
  /* real part */
  ddmuld(dda_l[0], dda_t[0], db[0], &t1_l, &t1_t);
  ddmuld(dda_l[1], dda_t[1], -db[1], &t2_l, &t2_t);
  ddadd(t1_l, t1_t, t2_l, t2_t, &t_l, &t_t);
  ddc_l[0] = t_l;
  ddc_t[0] = t_t;
  /* imaginary part */
  ddmuld(dda_l[1], dda_t[1], db[0], &t1_l, &t1_t);
  ddmuld(dda_l[0], dda_t[0], db[1], &t2_l, &t2_t);
  ddadd(t1_l, t1_t, t2_l, t2_t, &t_l, &t_t);
  ddc_l[1] = t_l;
  ddc_t[1] = t_t;
}

static void z_dddivd(double *dda_l, double *dda_t, double *db,
                     double *ddc_l, double *ddc_t)
{
/* Purpose
 * =======
 * 
 * This divides the complex DD number DDA by the complex DP number DB to
 * yield the complex DD quotient DDC.
 * 
 */
  double db_conj[2];
  double t_l[2], t_t[2];
  double d_l, d_t;

  /* b_r^2 + b_i^2 in double-double */
  d_l = db[0];
  d_t = 0.0;
  ddmuld(d_l, d_t, d_l, &t_l[0], &t_t[0]);
  d_l = db[1];
  d_t = 0.0;
  ddmuld(d_l, d_t, d_l, &t_l[1], &t_t[1]);
  ddadd(t_l[0], t_t[0], t_l[1], t_t[1], &d_l, &d_t);

  db_conj[0] = db[0];
  db_conj[1] = -db[1];
  z_ddmuld(dda_l, dda_t, db_conj, t_l, t_t);
  /*printf("\tz_dddivd() -> a * conj(b) = %10.8e + %10.8ei\n",
     t_l[0],t_l[1]); */
  dddiv(t_l[0], t_t[0], d_l, d_t, &ddc_l[0], &ddc_t[0]);
  dddiv(t_l[1], t_t[1], d_l, d_t, &ddc_l[1], &ddc_t[1]);
}

static void BLAS_zdot_x(enum blas_conj_type conj, int n, const void *alpha,
                        const void *x, int incx, const void *beta,
                        const void *y, int incy, void *r, enum blas_prec_type prec)

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
 * alpha  (input) const void*
 * 
 * x      (input) const void*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) const void*
 *
 * y      (input) const void*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) void*
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
  static const char routine_name[] = "BLAS_zdot_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      double *r_i = (double *) r;
      const double *x_i = (double *) x;
      const double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double x_ii[2];
      double y_ii[2];
      double r_v[2];
      double prod[2];
      double sum[2];
      double tmp1[2];
      double tmp2[2];


      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if (((beta_i[0] == 1.0 && beta_i[1] == 0.0))
	  && (n == 0 || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)))
	return;



      r_v[0] = r_i[0];
      r_v[1] = r_i[0 + 1];
      sum[0] = sum[1] = 0.0;
      incx *= 2;
      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      if (conj == blas_conj) {
	for (i = 0; i < n; ++i) {
	  x_ii[0] = x_i[ix];
	  x_ii[1] = x_i[ix + 1];
	  y_ii[0] = y_i[iy];
	  y_ii[1] = y_i[iy + 1];
	  x_ii[1] = -x_ii[1];
	  {
	    prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
	    prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
	  }			/* prod = x[i]*y[i] */
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];	/* sum = sum+prod */
	  ix += incx;
	  iy += incy;
	}			/* endfor */
      } else {
	/* do not conjugate */

	for (i = 0; i < n; ++i) {
	  x_ii[0] = x_i[ix];
	  x_ii[1] = x_i[ix + 1];
	  y_ii[0] = y_i[iy];
	  y_ii[1] = y_i[iy + 1];

	  {
	    prod[0] = (double) x_ii[0] * y_ii[0] - (double) x_ii[1] * y_ii[1];
	    prod[1] = (double) x_ii[0] * y_ii[1] + (double) x_ii[1] * y_ii[0];
	  }			/* prod = x[i]*y[i] */
	  sum[0] = sum[0] + prod[0];
	  sum[1] = sum[1] + prod[1];	/* sum = sum+prod */
	  ix += incx;
	  iy += incy;
	}			/* endfor */
      }

      {
	tmp1[0] = (double) sum[0] * alpha_i[0] - (double) sum[1] * alpha_i[1];
	tmp1[1] = (double) sum[0] * alpha_i[1] + (double) sum[1] * alpha_i[0];
      }				/* tmp1 = sum*alpha */
      {
	tmp2[0] = (double) r_v[0] * beta_i[0] - (double) r_v[1] * beta_i[1];
	tmp2[1] = (double) r_v[0] * beta_i[1] + (double) r_v[1] * beta_i[0];
      }				/* tmp2 = r*beta */
      tmp1[0] = tmp1[0] + tmp2[0];
      tmp1[1] = tmp1[1] + tmp2[1];	/* tmp1 = tmp1+tmp2 */
      ((double *) r)[0] = tmp1[0];
      ((double *) r)[1] = tmp1[1];	/* r = tmp1 */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      double *r_i = (double *) r;
      const double *x_i = (double *) x;
      const double *y_i = (double *) y;
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;
      double x_ii[2];
      double y_ii[2];
      double r_v[2];
      double head_prod[2], tail_prod[2];
      double head_sum[2], tail_sum[2];
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if (((beta_i[0] == 1.0 && beta_i[1] == 0.0))
	  && (n == 0 || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)))
	return;

      FPU_FIX_START;

      r_v[0] = r_i[0];
      r_v[1] = r_i[0 + 1];
      head_sum[0] = head_sum[1] = tail_sum[0] = tail_sum[1] = 0.0;
      incx *= 2;
      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      if (conj == blas_conj) {
	for (i = 0; i < n; ++i) {
	  x_ii[0] = x_i[ix];
	  x_ii[1] = x_i[ix + 1];
	  y_ii[0] = y_i[iy];
	  y_ii[1] = y_i[iy + 1];
	  x_ii[1] = -x_ii[1];
	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[0] * split;
	      a1 = con - x_ii[0];
	      a1 = con - a1;
	      a2 = x_ii[0] - a1;
	      con = y_ii[0] * split;
	      b1 = con - y_ii[0];
	      b1 = con - b1;
	      b2 = y_ii[0] - b1;

	      head_t1 = x_ii[0] * y_ii[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[1] * split;
	      a1 = con - x_ii[1];
	      a1 = con - a1;
	      a2 = x_ii[1] - a1;
	      con = y_ii[1] * split;
	      b1 = con - y_ii[1];
	      b1 = con - b1;
	      b2 = y_ii[1] - b1;

	      head_t2 = x_ii[1] * y_ii[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_prod[0] = head_t1;
	    tail_prod[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[1] * split;
	      a1 = con - x_ii[1];
	      a1 = con - a1;
	      a2 = x_ii[1] - a1;
	      con = y_ii[0] * split;
	      b1 = con - y_ii[0];
	      b1 = con - b1;
	      b2 = y_ii[0] - b1;

	      head_t1 = x_ii[1] * y_ii[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[0] * split;
	      a1 = con - x_ii[0];
	      a1 = con - a1;
	      a2 = x_ii[0] - a1;
	      con = y_ii[1] * split;
	      b1 = con - y_ii[1];
	      b1 = con - b1;
	      b2 = y_ii[1] - b1;

	      head_t2 = x_ii[0] * y_ii[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_prod[1] = head_t1;
	    tail_prod[1] = tail_t1;
	  }			/* prod = x[i]*y[i] */
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum[0];
	    tail_a = tail_sum[0];
	    head_b = head_prod[0];
	    tail_b = tail_prod[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum[0] = head_t;
	    tail_sum[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum[1];
	    tail_a = tail_sum[1];
	    head_b = head_prod[1];
	    tail_b = tail_prod[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum[1] = head_t;
	    tail_sum[1] = tail_t;
	  }			/* sum = sum+prod */
	  ix += incx;
	  iy += incy;
	}			/* endfor */
      } else {
	/* do not conjugate */

	for (i = 0; i < n; ++i) {
	  x_ii[0] = x_i[ix];
	  x_ii[1] = x_i[ix + 1];
	  y_ii[0] = y_i[iy];
	  y_ii[1] = y_i[iy + 1];

	  {
	    /* Compute complex-extra = complex-double * complex-double. */
	    double head_t1, tail_t1;
	    double head_t2, tail_t2;
	    /* Real part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[0] * split;
	      a1 = con - x_ii[0];
	      a1 = con - a1;
	      a2 = x_ii[0] - a1;
	      con = y_ii[0] * split;
	      b1 = con - y_ii[0];
	      b1 = con - b1;
	      b2 = y_ii[0] - b1;

	      head_t1 = x_ii[0] * y_ii[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[1] * split;
	      a1 = con - x_ii[1];
	      a1 = con - a1;
	      a2 = x_ii[1] - a1;
	      con = y_ii[1] * split;
	      b1 = con - y_ii[1];
	      b1 = con - b1;
	      b2 = y_ii[1] - b1;

	      head_t2 = x_ii[1] * y_ii[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    head_t2 = -head_t2;
	    tail_t2 = -tail_t2;
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_prod[0] = head_t1;
	    tail_prod[0] = tail_t1;
	    /* Imaginary part */
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[1] * split;
	      a1 = con - x_ii[1];
	      a1 = con - a1;
	      a2 = x_ii[1] - a1;
	      con = y_ii[0] * split;
	      b1 = con - y_ii[0];
	      b1 = con - b1;
	      b2 = y_ii[0] - b1;

	      head_t1 = x_ii[1] * y_ii[0];
	      tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double_double = double * double. */
	      double a1, a2, b1, b2, con;

	      con = x_ii[0] * split;
	      a1 = con - x_ii[0];
	      a1 = con - a1;
	      a2 = x_ii[0] - a1;
	      con = y_ii[1] * split;
	      b1 = con - y_ii[1];
	      b1 = con - b1;
	      b2 = y_ii[1] - b1;

	      head_t2 = x_ii[0] * y_ii[1];
	      tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_t1 + head_t2;
	      bv = s1 - head_t1;
	      s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_t1 + tail_t2;
	      bv = t1 - tail_t1;
	      t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t1 = t1 + t2;
	      tail_t1 = t2 - (head_t1 - t1);
	    }
	    head_prod[1] = head_t1;
	    tail_prod[1] = tail_t1;
	  }			/* prod = x[i]*y[i] */
	  {
	    double head_t, tail_t;
	    double head_a, tail_a;
	    double head_b, tail_b;
	    /* Real part */
	    head_a = head_sum[0];
	    tail_a = tail_sum[0];
	    head_b = head_prod[0];
	    tail_b = tail_prod[0];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum[0] = head_t;
	    tail_sum[0] = tail_t;
	    /* Imaginary part */
	    head_a = head_sum[1];
	    tail_a = tail_sum[1];
	    head_b = head_prod[1];
	    tail_b = tail_prod[1];
	    {
	      /* Compute double-double = double-double + double-double. */
	      double bv;
	      double s1, s2, t1, t2;

	      /* Add two hi words. */
	      s1 = head_a + head_b;
	      bv = s1 - head_a;
	      s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	      /* Add two lo words. */
	      t1 = tail_a + tail_b;
	      bv = t1 - tail_a;
	      t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	      s2 += t1;

	      /* Renormalize (s1, s2)  to  (t1, s2) */
	      t1 = s1 + s2;
	      s2 = s2 - (t1 - s1);

	      t2 += s2;

	      /* Renormalize (t1, t2)  */
	      head_t = t1 + t2;
	      tail_t = t2 - (head_t - t1);
	    }
	    head_sum[1] = head_t;
	    tail_sum[1] = tail_t;
	  }			/* sum = sum+prod */
	  ix += incx;
	  iy += incy;
	}			/* endfor */
      }

      {
	/* Compute complex-extra = complex-extra * complex-double. */
	double head_a0, tail_a0;
	double head_a1, tail_a1;
	double head_t1, tail_t1;
	double head_t2, tail_t2;
	head_a0 = head_sum[0];
	tail_a0 = tail_sum[0];
	head_a1 = head_sum[1];
	tail_a1 = tail_sum[1];
	/* real part */
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_a0 * split;
	  a11 = con - head_a0;
	  a11 = con - a11;
	  a21 = head_a0 - a11;
	  con = alpha_i[0] * split;
	  b1 = con - alpha_i[0];
	  b1 = con - b1;
	  b2 = alpha_i[0] - b1;

	  c11 = head_a0 * alpha_i[0];
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_a0 * alpha_i[0];
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_a1 * split;
	  a11 = con - head_a1;
	  a11 = con - a11;
	  a21 = head_a1 - a11;
	  con = alpha_i[1] * split;
	  b1 = con - alpha_i[1];
	  b1 = con - b1;
	  b2 = alpha_i[1] - b1;

	  c11 = head_a1 * alpha_i[1];
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_a1 * alpha_i[1];
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_t2 = t1 + t2;
	  tail_t2 = t2 - (head_t2 - t1);
	}
	head_t2 = -head_t2;
	tail_t2 = -tail_t2;
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_tmp1[0] = head_t1;
	tail_tmp1[0] = tail_t1;
	/* imaginary part */
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_a1 * split;
	  a11 = con - head_a1;
	  a11 = con - a11;
	  a21 = head_a1 - a11;
	  con = alpha_i[0] * split;
	  b1 = con - alpha_i[0];
	  b1 = con - b1;
	  b2 = alpha_i[0] - b1;

	  c11 = head_a1 * alpha_i[0];
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_a1 * alpha_i[0];
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_a0 * split;
	  a11 = con - head_a0;
	  a11 = con - a11;
	  a21 = head_a0 - a11;
	  con = alpha_i[1] * split;
	  b1 = con - alpha_i[1];
	  b1 = con - b1;
	  b2 = alpha_i[1] - b1;

	  c11 = head_a0 * alpha_i[1];
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_a0 * alpha_i[1];
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_t2 = t1 + t2;
	  tail_t2 = t2 - (head_t2 - t1);
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_tmp1[1] = head_t1;
	tail_tmp1[1] = tail_t1;
      }
      /* tmp1 = sum*alpha */
      {
	/* Compute complex-extra = complex-double * complex-double. */
	double head_t1, tail_t1;
	double head_t2, tail_t2;
	/* Real part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = r_v[0] * split;
	  a1 = con - r_v[0];
	  a1 = con - a1;
	  a2 = r_v[0] - a1;
	  con = beta_i[0] * split;
	  b1 = con - beta_i[0];
	  b1 = con - b1;
	  b2 = beta_i[0] - b1;

	  head_t1 = r_v[0] * beta_i[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = r_v[1] * split;
	  a1 = con - r_v[1];
	  a1 = con - a1;
	  a2 = r_v[1] - a1;
	  con = beta_i[1] * split;
	  b1 = con - beta_i[1];
	  b1 = con - b1;
	  b2 = beta_i[1] - b1;

	  head_t2 = r_v[1] * beta_i[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	head_t2 = -head_t2;
	tail_t2 = -tail_t2;
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_tmp2[0] = head_t1;
	tail_tmp2[0] = tail_t1;
	/* Imaginary part */
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = r_v[1] * split;
	  a1 = con - r_v[1];
	  a1 = con - a1;
	  a2 = r_v[1] - a1;
	  con = beta_i[0] * split;
	  b1 = con - beta_i[0];
	  b1 = con - b1;
	  b2 = beta_i[0] - b1;

	  head_t1 = r_v[1] * beta_i[0];
	  tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = r_v[0] * split;
	  a1 = con - r_v[0];
	  a1 = con - a1;
	  a2 = r_v[0] - a1;
	  con = beta_i[1] * split;
	  b1 = con - beta_i[1];
	  b1 = con - b1;
	  b2 = beta_i[1] - b1;

	  head_t2 = r_v[0] * beta_i[1];
	  tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_t1 + head_t2;
	  bv = s1 - head_t1;
	  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_t1 + tail_t2;
	  bv = t1 - tail_t1;
	  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t1 = t1 + t2;
	  tail_t1 = t2 - (head_t1 - t1);
	}
	head_tmp2[1] = head_t1;
	tail_tmp2[1] = tail_t1;
      }				/* tmp2 = r*beta */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_tmp1[0];
	tail_a = tail_tmp1[0];
	head_b = head_tmp2[0];
	tail_b = tail_tmp2[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_tmp1[0] = head_t;
	tail_tmp1[0] = tail_t;
	/* Imaginary part */
	head_a = head_tmp1[1];
	tail_a = tail_tmp1[1];
	head_b = head_tmp2[1];
	tail_b = tail_tmp2[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_tmp1[1] = head_t;
	tail_tmp1[1] = tail_t;
      }				/* tmp1 = tmp1+tmp2 */
      ((double *) r)[0] = head_tmp1[0];
      ((double *) r)[1] = head_tmp1[1];	/* r = tmp1 */

      FPU_FIX_STOP;
    }
    break;
  }
}


/* Complex-Complex Multiplication */
static void z_mul(double a[], double b[], double c[])
{
  double cr, ci;
  cr = a[0] * b[0] - a[1] * b[1];
  ci = a[1] * b[0] + a[0] * b[1];
  c[0] = cr;
  c[1] = ci;
}

/* Complex Division c = a/b */
static void z_div(double a[], double b[], double c[])
{
  double ratio, den;
  double abr, abi, cr, ci;

  if ((abr = b[0]) < 0.)
    abr = -abr;
  if ((abi = b[1]) < 0.)
    abi = -abi;
  if (abr <= abi) {
    if (abi == 0) {
      BLAS_error("z_div: division by zero", 0, 0, NULL);
    }
    ratio = b[0] / b[1];
    den = b[1] * (1 + ratio * ratio);
    cr = (a[0] * ratio + a[1]) / den;
    ci = (a[1] * ratio - a[0]) / den;
  } else {
    ratio = b[1] / b[0];
    den = b[0] * (1 + ratio * ratio);
    cr = (a[0] + a[1] * ratio) / den;
    ci = (a[1] - a[0] * ratio) / den;
  }
  c[0] = cr;
  c[1] = ci;
}

static double ulp(double a)
/*
 * Purpose 
 * =======
 * 
 * Compute the unit last place of a double precision number.
 */
{
  double f;
  int e;
  f = frexp(a, &e);
  return power(2, e - BITS_D);
}


static double rand_half_1(int l_bits, int *seed)
/*
 * Purpose
 * =======
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


static void r_truth(enum blas_conj_type conj, int n, void *alpha, const void *x, int incx, void *beta, const void *y, int incy, void *r,        /* input */
                    double *r_true_l, double *r_true_t)
{
  int i, ix = 0, iy = 0;
  double *r_i = (double *) r;
  const double *x_i = (double *) x;
  const double *y_i = (double *) y;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double x_ii[2];
  double y_ii[2];
  double r_v[2];
  double prod_l[2], prod_t[2];
  double sum_l[2], sum_t[2];
  double tmp1_l[2], tmp1_t[2];
  double tmp2_l[2], tmp2_t[2];

  /* Immediate return */
  if (n < 0) {
    r_true_l[0] = r_true_l[1] = r_true_t[0] = r_true_t[1] = 0.0;
    return;
  }

  r_v[0] = r_i[0];
  r_v[1] = r_i[0 + 1];
  sum_l[0] = sum_l[1] = sum_t[0] = sum_t[1] = 0.0;      /* sum = 0 */
  for (i = 0; i < n; ++i) {
    x_ii[0] = x_i[ix];
    x_ii[1] = x_i[ix + 1];
    if (conj == blas_conj)
      x_ii[1] = -x_ii[1];
    y_ii[0] = y_i[iy];
    y_ii[1] = y_i[iy + 1];
    {
      /*
       * Compute complex-extra = complex-double * complex-double.
       */
      double t1_l, t1_t;
      double t2_l, t2_t;
      /* Real part */
      {
        /* Compute double_double = double * double. */
        double a1, a2, b1, b2, con;

        con = x_ii[0] * split;
        a1 = con - x_ii[0];
        a1 = con - a1;
        a2 = x_ii[0] - a1;
        con = y_ii[0] * split;
        b1 = con - y_ii[0];
        b1 = con - b1;
        b2 = y_ii[0] - b1;

        t1_l = x_ii[0] * y_ii[0];
        t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      {
        /* Compute double_double = double * double. */
        double a1, a2, b1, b2, con;

        con = x_ii[1] * split;
        a1 = con - x_ii[1];
        a1 = con - a1;
        a2 = x_ii[1] - a1;
        con = y_ii[1] * split;
        b1 = con - y_ii[1];
        b1 = con - b1;
        b2 = y_ii[1] - b1;

        t2_l = x_ii[1] * y_ii[1];
        t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      t2_l = -t2_l;
      t2_t = -t2_t;
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = t1_l + t2_l;
        e = t1 - t1_l;
        t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

        /* The result is t1 + t2, after normalization. */
        t1_l = t1 + t2;
        t1_t = t2 - (t1_l - t1);
      }
      prod_l[0] = t1_l;
      prod_t[0] = t1_t;
      /* Imaginary part */
      {
        /* Compute double_double = double * double. */
        double a1, a2, b1, b2, con;

        con = x_ii[1] * split;
        a1 = con - x_ii[1];
        a1 = con - a1;
        a2 = x_ii[1] - a1;
        con = y_ii[0] * split;
        b1 = con - y_ii[0];
        b1 = con - b1;
        b2 = y_ii[0] - b1;

        t1_l = x_ii[1] * y_ii[0];
        t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      {
        /* Compute double_double = double * double. */
        double a1, a2, b1, b2, con;

        con = x_ii[0] * split;
        a1 = con - x_ii[0];
        a1 = con - a1;
        a2 = x_ii[0] - a1;
        con = y_ii[1] * split;
        b1 = con - y_ii[1];
        b1 = con - b1;
        b2 = y_ii[1] - b1;

        t2_l = x_ii[0] * y_ii[1];
        t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = t1_l + t2_l;
        e = t1 - t1_l;
        t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

        /* The result is t1 + t2, after normalization. */
        t1_l = t1 + t2;
        t1_t = t2 - (t1_l - t1);
      }
      prod_l[1] = t1_l;
      prod_t[1] = t1_t;
    }                           /* prod = x[i]*y[i] */
    {
      double t_l, t_t;
      double a_l, a_t;
      double b_l, b_t;
      /* Real part */
      a_l = sum_l[0];
      a_t = sum_t[0];
      b_l = prod_l[0];
      b_t = prod_t[0];
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = a_l + b_l;
        e = t1 - a_l;
        t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

        /* The result is t1 + t2, after normalization. */
        t_l = t1 + t2;
        t_t = t2 - (t_l - t1);
      }
      sum_l[0] = t_l;
      sum_t[0] = t_t;
      /* Imaginary part */
      a_l = sum_l[1];
      a_t = sum_t[1];
      b_l = prod_l[1];
      b_t = prod_t[1];
      {
        /*
         * Compute double-double = double-double + double-double.
         */
        double e, t1, t2;

        /* Knuth trick. */
        t1 = a_l + b_l;
        e = t1 - a_l;
        t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

        /* The result is t1 + t2, after normalization. */
        t_l = t1 + t2;
        t_t = t2 - (t_l - t1);
      }
      sum_l[1] = t_l;
      sum_t[1] = t_t;
    }                           /* sum = sum+prod */
    ix += 2;
    iy += 2;
  }                             /* endfor */
  {
    /* Compute complex-extra = complex-extra * complex-double. */
    double a0_l, a0_t;
    double a1_l, a1_t;
    double t1_l, t1_t;
    double t2_l, t2_t;
    a0_l = sum_l[0];
    a0_t = sum_t[0];
    a1_l = sum_l[1];
    a1_t = sum_t[1];
    /* real part */
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a0_l * split;
      a11 = con - a0_l;
      a11 = con - a11;
      a21 = a0_l - a11;
      con = alpha_i[0] * split;
      b1 = con - alpha_i[0];
      b1 = con - b1;
      b2 = alpha_i[0] - b1;

      c11 = a0_l * alpha_i[0];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a0_t * alpha_i[0];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a1_l * split;
      a11 = con - a1_l;
      a11 = con - a11;
      a21 = a1_l - a11;
      con = alpha_i[1] * split;
      b1 = con - alpha_i[1];
      b1 = con - b1;
      b2 = alpha_i[1] - b1;

      c11 = a1_l * alpha_i[1];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a1_t * alpha_i[1];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t2_l = t1 + t2;
      t2_t = t2 - (t2_l - t1);
    }
    t2_l = -t2_l;
    t2_t = -t2_t;
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp1_l[0] = t1_l;
    tmp1_t[0] = t1_t;
    /* imaginary part */
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a1_l * split;
      a11 = con - a1_l;
      a11 = con - a11;
      a21 = a1_l - a11;
      con = alpha_i[0] * split;
      b1 = con - alpha_i[0];
      b1 = con - b1;
      b2 = alpha_i[0] - b1;

      c11 = a1_l * alpha_i[0];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a1_t * alpha_i[0];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    {
      /* Compute double-double = double-double * double. */
      double a11, a21, b1, b2, c11, c21, c2, con, e, t1, t2;

      con = a0_l * split;
      a11 = con - a0_l;
      a11 = con - a11;
      a21 = a0_l - a11;
      con = alpha_i[1] * split;
      b1 = con - alpha_i[1];
      b1 = con - b1;
      b2 = alpha_i[1] - b1;

      c11 = a0_l * alpha_i[1];
      c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

      c2 = a0_t * alpha_i[1];
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      t2_l = t1 + t2;
      t2_t = t2 - (t2_l - t1);
    }
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp1_l[1] = t1_l;
    tmp1_t[1] = t1_t;
  }
  /* tmp1 = sum*alpha */
  {
    /* Compute complex-extra = complex-double * complex-double. */
    double t1_l, t1_t;
    double t2_l, t2_t;
    /* Real part */
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[0] * split;
      a1 = con - r_v[0];
      a1 = con - a1;
      a2 = r_v[0] - a1;
      con = beta_i[0] * split;
      b1 = con - beta_i[0];
      b1 = con - b1;
      b2 = beta_i[0] - b1;

      t1_l = r_v[0] * beta_i[0];
      t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[1] * split;
      a1 = con - r_v[1];
      a1 = con - a1;
      a2 = r_v[1] - a1;
      con = beta_i[1] * split;
      b1 = con - beta_i[1];
      b1 = con - b1;
      b2 = beta_i[1] - b1;

      t2_l = r_v[1] * beta_i[1];
      t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    t2_l = -t2_l;
    t2_t = -t2_t;
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp2_l[0] = t1_l;
    tmp2_t[0] = t1_t;
    /* Imaginary part */
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[1] * split;
      a1 = con - r_v[1];
      a1 = con - a1;
      a2 = r_v[1] - a1;
      con = beta_i[0] * split;
      b1 = con - beta_i[0];
      b1 = con - b1;
      b2 = beta_i[0] - b1;

      t1_l = r_v[1] * beta_i[0];
      t1_t = (((a1 * b1 - t1_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double_double = double * double. */
      double a1, a2, b1, b2, con;

      con = r_v[0] * split;
      a1 = con - r_v[0];
      a1 = con - a1;
      a2 = r_v[0] - a1;
      con = beta_i[1] * split;
      b1 = con - beta_i[1];
      b1 = con - b1;
      b2 = beta_i[1] - b1;

      t2_l = r_v[0] * beta_i[1];
      t2_t = (((a1 * b1 - t2_l) + a1 * b2) + a2 * b1) + a2 * b2;
    }
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = t1_l + t2_l;
      e = t1 - t1_l;
      t2 = ((t2_l - e) + (t1_l - (t1 - e))) + t1_t + t2_t;

      /* The result is t1 + t2, after normalization. */
      t1_l = t1 + t2;
      t1_t = t2 - (t1_l - t1);
    }
    tmp2_l[1] = t1_l;
    tmp2_t[1] = t1_t;
  }                             /* tmp2 = r*beta */
  {
    double t_l, t_t;
    double a_l, a_t;
    double b_l, b_t;
    /* Real part */
    a_l = tmp1_l[0];
    a_t = tmp1_t[0];
    b_l = tmp2_l[0];
    b_t = tmp2_t[0];
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = a_l + b_l;
      e = t1 - a_l;
      t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

      /* The result is t1 + t2, after normalization. */
      t_l = t1 + t2;
      t_t = t2 - (t_l - t1);
    }
    tmp1_l[0] = t_l;
    tmp1_t[0] = t_t;
    /* Imaginary part */
    a_l = tmp1_l[1];
    a_t = tmp1_t[1];
    b_l = tmp2_l[1];
    b_t = tmp2_t[1];
    {
      /* Compute double-double = double-double + double-double. */
      double e, t1, t2;

      /* Knuth trick. */
      t1 = a_l + b_l;
      e = t1 - a_l;
      t2 = ((b_l - e) + (a_l - (t1 - e))) + a_t + b_t;

      /* The result is t1 + t2, after normalization. */
      t_l = t1 + t2;
      t_t = t2 - (t_l - t1);
    }
    tmp1_l[1] = t_l;
    tmp1_t[1] = t_t;
  }                             /* tmp1 = tmp1+tmp2 */

  /* Return r_truth = tmp1 */
  r_true_l[0] = tmp1_l[0];
  r_true_l[1] = tmp1_l[1];
  r_true_t[0] = tmp1_t[0];
  r_true_t[1] = tmp1_t[1];
}                               /* end r_truth */


static void
gen_y_to_cancel(int k, int n, enum blas_conj_type conj,
                void *alpha, void *x, void *y)
/*
 * Purpose
 * =======
 * 
 * Generate Y(i)'s from k to n-1 to cancel as much as possible.
 * 
 */
{
  int i, ii;
  double zero[2] = { 0.0, 0.0 };
  double r[2] = { 0.0, 0.0 };
  double tmp[2], tmp_l[2], tmp_t[2];
  double *x_i = x, *y_i = y;
  double r_true_l[2], r_true_t[2];

  for (i = k; i < n; ++i) {
    /* y[i] = -rtmp / (alpha * x[i]); */
    r_truth(conj, i, alpha, x, 1, zero, y, 1, r, r_true_l, r_true_t);
    ii = 2 * i;
    if (conj == blas_conj) {
      tmp[0] = x_i[ii];
      tmp[1] = -x_i[ii + 1];
      z_mul((double *) alpha, tmp, tmp);
    } else {
      z_mul((double *) alpha, &x_i[ii], tmp);
    }
    if (tmp[0] == 0. && tmp[1] == 0.)
      y_i[ii] = y_i[ii + 1] = 0.;
    else {
      z_dddivd(r_true_l, r_true_t, tmp, tmp_l, tmp_t);
      y_i[ii] = -tmp_l[0];
      y_i[ii + 1] = -tmp_l[1];
    }
  }
}


static void
gen_r_to_cancel(int n, enum blas_conj_type conj,
                void *alpha, void *beta, void *x, void *y, void *r, int *seed)
/*
 * Purpose
 * =======
 * 
 * Generate r to cancel as much as possible.
 * 
 */
{
  double zero[2] = { 0.0, 0.0 };
  double rtmp[2] = { 0.0, 0.0 };
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double r_true_l[2], r_true_t[2];

  if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
    r_i[0] = xrand(seed);
    r_i[1] = xrand(seed);
  } else {
    r_truth(conj, n, alpha, x, 1, zero, y, 1, rtmp, r_true_l, r_true_t);
    z_dddivd(r_true_l, r_true_t, beta_i, r_true_l, r_true_t);
    r_i[0] = -r_true_l[0];
    r_i[1] = -r_true_l[1];
  }
}



void
util_xblas_zdot_fill(int n, int n_fix2, int n_mix, int norm,
                  char charconj,
                  void *alpha, int alpha_flag, void *beta, int beta_flag,
                  void *x, void *y, int *seed,
                  void *r, double r_true_l[], double r_true_t[])
/*
 * Purpose
 * =======
 *
 * This routine generates the test vectors X and Y for C_ZDOT.
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
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int* 
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double[]
 *         The leading (real,imaginary) parts of the truth in double-double.
 *
 * r_true_t (output) double[]
 *         The trailing (real,imaginary) parts of the truth in double-double.
 *
 */
{
  int B, frees, y_free, i, ii, k, s;
  double zero[2] = { 0.0, 0.0 };
  double a, b, eps, eps_out;
  double f[2], rtmp[2];
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;
  double *r_i = (double *) r;
  double *x_i = (double *) x, *y_i = (double *) y;

  enum blas_conj_type conj = blas_no_conj;
  if(charconj == 'c' || charconj == 'C'){
    conj = blas_conj;
  }

  if (alpha_flag == 0) {
    alpha_i[0] = xrand(seed);
    alpha_i[1] = xrand(seed);
  }
  if (beta_flag == 0) {
    beta_i[0] = xrand(seed);
    beta_i[1] = xrand(seed);
  }

  y_free = n - n_fix2;
  k = 2 * n_fix2;
  eps_out = power(2, -BITS_D);

  /*
   * Compute the number of bits in the prefix sum:
   *     alpha * SUM_{i=0,n_fix2-1}(x[i] * y[i])
   */
  r_i[0] = r_i[1] = 0.0;
  r_truth(conj, n_fix2, alpha, x, 1, zero, y, 1, r, r_true_l, r_true_t);
  B = FixedBits(r_true_l[0], r_true_t[0]);      /* real */
  B = MAX(B, FixedBits(r_true_l[1], r_true_t[1]));      /* imag */

  /* Pick r at random */
  r_i[0] = xrand(seed);
  r_i[1] = xrand(seed);

  /* Pick the free X(i)'s at random. */
  for (i = n_fix2 + n_mix; i < n; ++i) {
    ii = 2 * i;
    x_i[ii] = xrand(seed);
    x_i[ii + 1] = xrand(seed);
  }

  if (alpha_flag == 1 && alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    /* Pick the free Y(i)'s at random. */
    for (i = n_fix2; i < n; ++i) {
      ii = 2 * i;
      y_i[ii] = xrand(seed);
      y_i[ii + 1] = xrand(seed);
    }
    /* Compute r_truth in double-double */
    r_truth(conj, n, alpha, x, 1, beta, y, 1, r, r_true_l, r_true_t);
    return;
  }

  if (beta_flag == 1 && beta_i[0] == 0.0 && beta_i[1] == 0.0) {
    if (B == 0) {               /* Assume alpha is not very big . */
      switch (y_free) {
      case 0:
        break;
      case 1:
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        break;
      case 2:
        /*
         * Make SUM_{i=0,1}(x[k+i] * y[k+i]) small ... alpha * eps^2
         */
        if (n_mix == 0) {       /* Both x[k] and x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = eps_out^2, small.
           * For complex, each real number is multiplied by (i+1),
           * the result is 2i * eps_out^2.
           */
          a = rand_half_1(26, seed);    /* strictly < 1 */
          x_i[k] = a;           /* real */
          x_i[k + 1] = a;       /* imag */
          y_i[k] = a;
          y_i[k + 1] = a;
          x_i[k + 2] = a + eps_out;     /* exact */
          x_i[k + 3] = a + eps_out;     /* exact */
          y_i[k + 2] = -a + eps_out;    /* exact */
          y_i[k + 3] = -a + eps_out;    /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed, x[k+1] free. */
          /* x[k]*y[k] + x[k+1]*y[k+1] = (eps1 + eps2*i)^2 */
          a = x_i[k];           /* real */
          b = x_i[k + 1];       /* imag */
          if (conj == blas_conj)
            b = -b;
          y_i[k] = a;
          y_i[k + 1] = b;
          eps = ulp(a);
          x_i[k + 2] = a + eps; /* exact */
          y_i[k + 2] = -a + eps;        /* exact */
          eps = ulp(b);
          x_i[k + 3] = b + eps; /* exact */
          if (conj == blas_conj)
            x_i[k + 3] = -x_i[k + 3];
          y_i[k + 3] = -b + eps;        /* exact */
        } else {                /* Both x[k] and x[k+1] fixed; cancel 53 bits. */
          y_i[k] = xrand(seed);
          y_i[k + 1] = xrand(seed);
          gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
        }
        break;
      case 3:
        /*
         * Make SUM_{i=0,2}(x[k+i] * y[k+i]) small
         * ... x[k]*y[k] = -x[k+2]*y[k+2] exact, x[k+1]*y[k+1] small
         */
        y_i[k] = -x_i[k + 4];
        y_i[k + 1] = -x_i[k + 5];
        y_i[k + 4] = x_i[k];
        y_i[k + 5] = x_i[k + 1];
        rtmp[0] = x_i[k];
        if (conj == blas_conj) {
          rtmp[1] = -x_i[k + 1];
          y_i[k + 1] = -y_i[k + 1];
          y_i[k + 5] = -y_i[k + 5];
        } else {
          rtmp[1] = x_i[k + 1];
        }
        z_mul(rtmp, &y_i[k], f);
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k + 2];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 3];
        else
          rtmp[1] = x_i[k + 3];
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k + 2] = y_i[k + 3] = 0.;
        else
          z_div(f, rtmp, &y_i[k + 2]);
        break;
      case 4:
        /*
         * Make SUM_{i=0,3}(x[k+i] * y[k+i]) small
         * ... x[k]*y[k] = -x[k+3]*y[k+3] exact, x[k+1]*y[k+1] small,
         *     x[k+2]*y[k+2] = 0.
         */
        y_i[k] = -x_i[k + 6];
        y_i[k + 1] = -x_i[k + 7];
        y_i[k + 6] = x_i[k];
        y_i[k + 7] = x_i[k + 1];
        rtmp[0] = x_i[k];
        if (conj == blas_conj) {
          rtmp[1] = -x_i[k + 1];
          y_i[k + 1] = -y_i[k + 1];
          y_i[k + 7] = -y_i[k + 7];
        } else {
          rtmp[1] = x_i[k + 1];
        }
        z_mul(rtmp, &y_i[k], f);
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k + 2];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 3];
        else
          rtmp[1] = x_i[k + 3];
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k + 2] = y_i[k + 3] = 0.;
        else
          z_div(f, rtmp, &y_i[k + 2]);
        y_i[k + 4] = 0.0;
        y_i[k + 5] = 0.0;
        break;
      default:                 /* y_free >= 5 */
        /*
         * Make SUM_{i=0,n-1}(x[k+i] * y[k+i]) small.
         * Use last 2 to cancel bits, and leading ones to add bits.
         * ... Cancel >= 106 bits.
         */
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        rtmp[0] = x_i[k];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 1];
        else
          rtmp[1] = x_i[k + 1];
        z_mul(rtmp, &y_i[k], rtmp);
        z_mul(alpha_i, rtmp, rtmp);
        s = 30;
        b = power(2, -s);
        for (i = n_fix2 + 1; i < n - 2; ++i) {
          rtmp[0] *= b;
          rtmp[1] *= b;
          ii = 2 * i;
          if (conj == blas_conj) {
            f[0] = x_i[ii];
            f[1] = -x_i[ii + 1];
          } else {
            f[0] = x_i[ii];
            f[1] = x_i[ii + 1];
          }
          z_mul(alpha_i, f, f);
          if (f[0] == 0. && f[1] == 0.)
            y_i[ii] = y_i[ii + 1] = 0.;
          else
            z_div(rtmp, f, &y_i[ii]);
        }
        gen_y_to_cancel(n - 2, n, conj, alpha, x, y);
      }                         /* end switch */
    } else {                    /* B > 0 */
      if (B >= BITS_E) {        /* Choose Y(i)'s to cancel. */
        gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
      } else {                  /* At least y[n-1] is free. */
        if (y_free == 1) {
          /* Cancel min(B,53) bits. */
          gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
        } else {                /* >= 2 frees. */
          /*
           * There are 2 possibilities:
           * (1) use 1 to add bits, and y_free-1 to cancel 53*(y_free-1) bits
           * (2) use all to cancel min(B, 53*y_free) bits
           * Goal is to maximize the # of bits cancelled. By equating (1)
           * and (2), we find the crossover point is y_free = B/53 + 1.
           */
          if (y_free > B / (double) BITS_D + 1) {       /* Use scheme (1) */
            f[0] = f[1] = 0.0;
            r_truth(conj, n_fix2, alpha, x, 1, zero, y, 1, f,
                    r_true_l, r_true_t);
            f[0] = r_true_l[0];
            f[1] = r_true_l[1];
            s = 100;            /* Should be random between 40 and 100. */
            b = power(2, -s);
            f[0] *= b;
            f[1] *= b;
            rtmp[0] = x_i[k];
            if (conj == blas_conj)
              rtmp[1] = -x_i[k + 1];
            else
              rtmp[1] = x_i[k + 1];
            z_mul(alpha_i, rtmp, rtmp);
            if (rtmp[0] == 0. && rtmp[1] == 0.)
              y_i[k] = y_i[k + 1] = 0.;
            else
              z_div(f, rtmp, &y_i[k]);
            gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
          } else {              /* Use scheme (2) */
            gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
          }
        }
      }                         /* end else B < 106 */
    }                           /* end else B > 0 */

    /* Compute r_truth in double-double */
    r_truth(conj, n, alpha, x, 1, zero, y, 1, r, r_true_l, r_true_t);
    return;
  }

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
        /* alpha*x[k]*y[k] + beta*r = -alpha * eps_out^2
         * For complex, each real number is multiplied by (i+1) to
         * yield final result -2i * alpha * eps_out^2.
         */
        a = rand_half_1(26, seed);      /* strictly < 1, only leading 26 bits */
        r_i[0] = 0.0;           /* real */
        r_i[1] = -a * a * 2;    /* imag, exact */
        if (beta_flag == 1) {   /* beta fixed, alpha must be free */
          alpha_i[0] = beta_i[0];
          alpha_i[1] = beta_i[1];
          x_i[k] = a + eps_out; /* real, exact */
          x_i[k + 1] = a + eps_out;     /* imag, exact */
          if (conj == blas_conj)
            x_i[k + 1] = -x_i[k + 1];
          y_i[k] = a - eps_out; /* exact */
          y_i[k + 1] = a - eps_out;     /* exact */
        } else if (n_mix == 1) {        /* x[k] fixed */
          beta_i[0] = x_i[k];
          beta_i[1] = x_i[k + 1];
          if (conj == blas_conj)
            beta_i[1] = -beta_i[1];
          alpha_i[0] = a + eps_out;     /* real, exact */
          alpha_i[1] = a + eps_out;     /* imag, exact */
          y_i[k] = a - eps_out; /* exact */
          y_i[k + 1] = a - eps_out;     /* exact */
        } else {                /* alpha fixed or free, x[k] and beta free */
          beta_i[0] = alpha_i[0];
          beta_i[1] = alpha_i[1];
          x_i[k] = a + eps_out; /* real, exact */
          x_i[k + 1] = a + eps_out;     /* imag, exact */
          if (conj == blas_conj)
            x_i[k + 1] = -x_i[k + 1];
          y_i[k] = a - eps_out; /* exact */
          y_i[k + 1] = a - eps_out;     /* exact */
        }
      } else {                  /* Cancel 53 bits. */
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
      }
      break;
    case 2:                    /* Actual frees = 3 */
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
        y_i[k] = -1.0;
        y_i[k + 1] = 0.0;
        if (alpha_flag == 0) {  /* alpha free */
          alpha_i[0] = beta_i[0];
          alpha_i[1] = beta_i[1];
          r_i[0] = x_i[k];
          r_i[1] = x_i[k + 1];
          if (conj == blas_conj)
            r_i[1] = -r_i[1];
        } else if (beta_flag == 0) {    /* beta free */
          beta_i[0] = alpha_i[0];
          beta_i[1] = alpha_i[1];
          r_i[0] = x_i[k];
          r_i[1] = x_i[k + 1];
          if (conj == blas_conj)
            r_i[1] = -r_i[1];
        } else {                /* x[k] free */
          x_i[k] = beta_i[0];
          x_i[k + 1] = beta_i[1];
          if (conj == blas_conj)
            x_i[k + 1] = -x_i[k + 1];
          r_i[0] = alpha_i[0];
          r_i[1] = alpha_i[1];
        }
        rtmp[0] = x_i[k];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 1];
        else
          rtmp[1] = x_i[k + 1];
        z_mul(rtmp, &y_i[k], f);
        z_mul(alpha_i, f, f);
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k + 2];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 3];
        else
          rtmp[1] = x_i[k + 3];
        z_mul(alpha_i, rtmp, rtmp);
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k + 2] = y_i[k + 3] = 0.;
        else
          z_div(f, rtmp, &y_i[k + 2]);
      } else {                  /* Cancel 53 bits. */
        y_i[k] = xrand(seed);
        y_i[k + 1] = xrand(seed);
        gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
        gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
      }
      break;
    default:                   /* Actual frees >= 4 */
      /*
       * Make SUM_{i=0,n-1}(alpha * x[k+i] * y[k+i]) + beta*r small.
       * Use last 2 ( Y(n) and r) to cancel bits, and leading ones to
       * add bits ... Cancel >= 106 bits.
       */
      y_i[k] = xrand(seed);
      y_i[k + 1] = xrand(seed);
      rtmp[0] = x_i[k];
      if (conj == blas_conj)
        rtmp[1] = -x_i[k + 1];
      else
        rtmp[1] = x_i[k + 1];
      z_mul(rtmp, &y_i[k], f);
      z_mul(alpha_i, f, f);
      s = 30;
      b = power(2, -s);
      for (i = n_fix2 + 1; i < n - 1; ++i) {
        f[0] *= b;
        f[1] *= b;
        ii = 2 * i;
        rtmp[0] = x_i[ii];
        if (conj == blas_conj)
          rtmp[1] = -x_i[ii + 1];
        else
          rtmp[1] = x_i[ii + 1];
        z_mul(alpha_i, rtmp, rtmp);
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[ii] = y_i[ii + 1] = 0.;
        else
          z_div(f, rtmp, &y_i[ii]);
      }
      gen_y_to_cancel(n - 1, n, conj, alpha, x, y);
      gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
      break;
    }                           /* end switch */
  } else {                      /* B > 0 */
    if (B >= BITS_E) {          /* Choose Y(i)'s and r to cancel. */
      gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
      gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
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
        f[0] = f[1] = 0.0;
        r_truth(conj, n_fix2, alpha, x, 1, zero, y, 1, f, r_true_l, r_true_t);
        f[0] = r_true_l[0];
        f[1] = r_true_l[1];
        s = 100;                /* Should be random between 40 and 100. */
        b = power(2, -s);
        f[0] *= b;
        f[1] *= b;
        rtmp[0] = x_i[k];
        if (conj == blas_conj)
          rtmp[1] = -x_i[k + 1];
        else
          rtmp[1] = x_i[k + 1];
        z_mul(alpha_i, rtmp, rtmp);
        if (rtmp[0] == 0. && rtmp[1] == 0.)
          y_i[k] = y_i[k + 1] = 0.;
        else
          z_div(f, rtmp, &y_i[k]);
        gen_y_to_cancel(n_fix2 + 1, n, conj, alpha, x, y);
      } else {                  /* Use scheme (2) */
        gen_y_to_cancel(n_fix2, n, conj, alpha, x, y);
      }
      gen_r_to_cancel(n, conj, alpha, beta, x, y, r, seed);
    }
  }

  /* Compute r_truth in double-double */
  r_truth(conj, n, alpha, x, 1, beta, y, 1, r, r_true_l, r_true_t);
}                               /* testgen_BLAS_zdot */
