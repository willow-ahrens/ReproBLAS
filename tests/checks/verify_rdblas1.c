#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "rdblas1_wrapper.h"

#include "../common/test_vecvec_fill_header.h"

#define MAX_NAME 100

int verify_rdblas1_reproducibility(int N, double* x, int incx, double* y, int incy, int func, double ref, Idouble Iref, int max_num_blocks) {
  // GENERATE DATA
  int i, j;
  double res;
  Idouble Ires;
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rdblas1_func(func))(N, x, incx, y, incy);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      dISetZero(Ires);
      for (j = 0; j < N; j += block_N) {
        block_N = block_N < N - j ? block_N : (N-j);
        dIAdd(&Ires, (wrap_Idblas1_func(func))(block_N, x + j * incx, incx, y + j * incy, incy));
      }
      res = Iconv2d(Ires);
    }
    if (res != ref) {
      printf("%s(x, y)[num_blocks=%d,block_N=%d] = %g != %g\n", wrap_rdblas1_name(func), num_blocks, block_N, res, ref);
      if (num_blocks == 1) {
        Ires = (wrap_Idblas1_func(func))(N, x, incx, y, incy);
      }
      printf("Ref I_double:\n");
      dIprint(Iref);
      printf("\nRes I_double:\n");
      dIprint(Ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  return 0;
}

extern const char* vecvec_fill_name(int argc, char** argv){
  static char namebuf[MAX_NAME];
  int func = opt_read_int(argc, argv, "-f", 0);
  snprintf(namebuf, MAX_NAME * sizeof(char), "Verify %s reproducibility", wrap_rdblas1_name(func));
  return namebuf;
}

extern int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type){
  int func = opt_read_int(argc, argv, "-f", 0);
  int rc = 0;
  double ref;
  Idouble Iref;
  int max_num_blocks = 1024;
  double *x = dvec_alloc(N, incx);
  double *y = dvec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  dvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  dvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  dvec_fill(N, x, incx, type, 1.0, opt_read_double(argc, argv, "-c", 1.0));

  //fill y with 1 where necessary
  dvec_fill(N, y, incy, vec_fill_CONSTANT, 1.0, 1.0);

  //nrm2 doesn't make sense with more than 1 block.
  if(func == verify_RDNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rdblas1_func(func))(N, x, incx, y, incy);
  Iref = (wrap_Idblas1_func(func))(N, x, incx, y, incy);

  dvec_reverse(N, x, incx);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_sort(N, x, incx, vec_order_INCREASING);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_sort(N, x, incx, vec_order_DECREASING);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_sort(N, x, incx, vec_order_INCREASING_MAGNITUDE);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_sort(N, x, incx, vec_order_DECREASING_MAGNITUDE);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_shuffle(N, x, incx);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_shuffle(N, x, incx);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_shuffle(N, x, incx);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  dvec_shuffle(N, x, incx);

  rc = verify_rdblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  return rc;
}
