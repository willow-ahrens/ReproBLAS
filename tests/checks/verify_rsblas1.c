#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "rsblas1_wrapper.h"
#include "test_vecvec_fill.h"

#define NAME_SIZE 100

int verify_rsblas1_reproducibility(int N, float* x, int incx, float* y, int incy, int func, float ref, Ifloat Iref, int max_num_blocks) {
  // GENERATE DATA
  int i, j;
  float res;
  Ifloat Ires;
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rsblas1_func(func))(N, x, incx, y, incy);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      sISetZero(Ires);
      for (j = 0; j < N; j += block_N) {
        block_N = block_N < N - j ? block_N : (N-j);
        sIAdd(&Ires, (wrap_Isblas1_func(func))(block_N, x + j * incx, incx, y + j * incy, incy));
      }
      res = Iconv2f(Ires);
    }
    if (res != ref) {
      printf("%s(x, y)[num_blocks=%d,block_N=%d] = %g != %g\n", wrap_rsblas1_name(func), num_blocks, block_N, res, ref);
      if (num_blocks == 1) {
        Ires = (wrap_Isblas1_func(func))(N, x, incx, y, incy);
      }
      printf("Ref I_float:\n");
      sIprint(Iref);
      printf("\nRes I_float:\n");
      sIprint(Ires);
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
  snprintf(namebuf, NAME_SIZE * sizeof(char), "Verify %s reproducibility", wrap_rsblas1_name(func));
  return namebuf;
}

extern int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type){
  int func = opt_read_int(argc, argv, "-f", 0);
  int rc = 0;
  float ref;
  Ifloat Iref;
  int max_num_blocks = 1024;
  float *x = svec_alloc(N, incx);
  float *y = svec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  svec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  svec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  svec_fill(N, x, incx, type, 1.0, opt_read_float(argc, argv, "-c", 1.0));

  //fill y with 1 where necessary
  svec_fill(N, y, incy, vec_fill_CONSTANT, 1.0, 1.0);

  //nrm2 doesn't make sense with more than 1 block.
  if(func == verify_RSNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rsblas1_func(func))(N, x, incx, y, incy);
  Iref = (wrap_Isblas1_func(func))(N, x, incx, y, incy);

  svec_reverse(N, x, incx);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_sort(N, x, incx, vec_order_INCREASING);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_sort(N, x, incx, vec_order_DECREASING);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_sort(N, x, incx, vec_order_INCREASING_MAGNITUDE);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_sort(N, x, incx, vec_order_DECREASING_MAGNITUDE);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_shuffle(N, x, incx);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_shuffle(N, x, incx);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_shuffle(N, x, incx);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  svec_shuffle(N, x, incx);

  rc = verify_rsblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  return rc;
}
