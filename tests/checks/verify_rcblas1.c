#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "rcblas1_wrapper.h"

#include "../common/test_vecvec_fill_header.h"

#define NAME_SIZE 100

int verify_rcblas1_reproducibility(int N, float complex* x, int incx, float complex* y, int incy, int func, float complex ref, I_float_Complex Iref, int max_num_blocks) {
  // GENERATE DATA
  int i, j;
  float complex res;
  I_float_Complex Ires;
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rcblas1_func(func))(N, x, incx, y, incy);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      cISetZero(Ires);
      for (j = 0; j < N; j += block_N) {
        block_N = block_N < N - j ? block_N : (N-j);
        cIAdd(&Ires, (wrap_Icblas1_func(func))(block_N, x + j * incx, incx, y + j * incy, incy));
      }
      res = Iconv2c(Ires);
    }
    if (res != ref) {
      printf("%s(x, y)[num_blocks=%d,block_N=%d] = %g + %gi != %g + %gi\n", wrap_rcblas1_name(func), num_blocks, block_N, CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
      if (num_blocks == 1) {
        Ires = (wrap_Icblas1_func(func))(N, x, incx, y, incy);
      }
      printf("Ref I_float_Complex:\n");
      cIprint(Iref);
      printf("\nRes I_float_Complex:\n");
      cIprint(Ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  return 0;
}

extern const char* vecvec_fill_name(int argc, char** argv){
  static char namebuf[NAME_SIZE];
  int func = opt_read_int(argc, argv, "-f", 0);
  snprintf(namebuf, NAME_SIZE * sizeof(char), "Verify %s reproducibility", wrap_rcblas1_name(func));
  return namebuf;
}

extern int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type){
  int func = opt_read_int(argc, argv, "-f", 0);
  int rc = 0;
  float complex ref;
  I_float_Complex Iref;
  int max_num_blocks = 1024;
  float complex *x = cvec_alloc(N, incx);
  float complex *y = cvec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  cvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  cvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  cvec_fill(N, x, incx, type, 1.0, opt_read_float(argc, argv, "-c", 1.0));

  //fill y with -i where necessary
  cvec_fill(N, y, incy, vec_fill_CONSTANT, -_Complex_I, 1.0);

  //nrm2 doesn't make sense with more than 1 block.
  if(func == verify_RSCNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rcblas1_func(func))(N, x, incx, y, incy);
  Iref = (wrap_Icblas1_func(func))(N, x, incx, y, incy);

  cvec_reverse(N, x, incx);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_sort(N, x, incx, vec_order_INCREASING);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_sort(N, x, incx, vec_order_DECREASING);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_sort(N, x, incx, vec_order_INCREASING_MAGNITUDE);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_sort(N, x, incx, vec_order_DECREASING_MAGNITUDE);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_shuffle(N, x, incx);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_shuffle(N, x, incx);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_shuffle(N, x, incx);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  cvec_shuffle(N, x, incx);

  rc = verify_rcblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  return rc;
}
