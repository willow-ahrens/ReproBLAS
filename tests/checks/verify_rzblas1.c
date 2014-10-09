#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "rzblas1_wrapper.h"

#define NAME_SIZE 100

static char namebuf[NAME_SIZE];

int verify_rzblas1_reproducibility(int N, double complex* x, int incx, double complex* y, int incy, int func, double complex ref, I_double_Complex Iref, int max_num_blocks) {
  // GENERATE DATA
  int i, j;
  double complex res;
  I_double_Complex Ires;
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rzblas1_func(func))(N, x, incx, y, incy);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      zISetZero(Ires);
      for (j = 0; j < N; j += block_N) {
        block_N = block_N < N - j ? block_N : (N-j);
        zIAdd(&Ires, (wrap_Izblas1_func(func))(block_N, x + j * incx, incx, y + j * incy, incy));
      }
      res = Iconv2z(Ires);
    }
    if (res != ref) {
      printf("%s(x, y)[num_blocks=%d,block_N=%d] = %g + %gi != %g + %gi\n", wrap_rzblas1_name(func), num_blocks, block_N, CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
      if (num_blocks == 1) {
        Ires = (wrap_Izblas1_func(func))(N, x, incx, y, incy);
      }
      printf("Ref I_double_Complex:\n");
      zIprint(Iref);
      printf("\nRes I_double_Complex:\n");
      zIprint(Ires);
      printf("\n");
      return 1;
    }
    num_blocks *= 2;
  }
  return 0;
}


extern const char* vecvec_fill_name(int argc, char** argv){
  int func = opt_read_int(argc, argv, "-f", 0);
  snprintf(namebuf, NAME_SIZE * sizeof(char), "Verify %s reproducibility", wrap_rzblas1_name(func));
  return namebuf;
}

extern int vecvec_fill_check(int argc, char** argv, int N, int incx, int incy, int type){
  int func = opt_read_int(argc, argv, "-f", 0);
  int rc = 0;
  double complex ref;
  I_double_Complex Iref;
  int max_num_blocks = 1024;
  double complex *x = zvec_alloc(N, incx);
  double complex *y = zvec_alloc(N, incy);

  vec_random_seed();

  //fill empty space with random data to check increment
  zvec_fill(N * incx, x, 1, vec_fill_RAND, 1.0, 1.0);
  zvec_fill(N * incy, y, 1, vec_fill_RAND, 1.0, 1.0);

  //fill x
  zvec_fill(N, x, incx, type, 1.0, opt_read_double(argc, argv, "-c", 1.0));

  //fill y with -i where necessary
  zvec_fill(N, y, incy, vec_fill_CONSTANT, -_Complex_I, 1.0);

  //nrm2 doesn't make sense with more than 1 block.
  if(func == verify_RDZNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rzblas1_func(func))(N, x, incx, y, incy);
  Iref = (wrap_Izblas1_func(func))(N, x, incx, y, incy);

  zvec_reverse(N, x, incx);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_sort(N, x, incx, vec_order_INCREASING);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_sort(N, x, incx, vec_order_DECREASING);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_sort(N, x, incx, vec_order_INCREASING_MAGNITUDE);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_sort(N, x, incx, vec_order_DECREASING_MAGNITUDE);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_shuffle(N, x, incx);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_shuffle(N, x, incx);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_shuffle(N, x, incx);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  zvec_shuffle(N, x, incx);

  rc = verify_rzblas1_reproducibility(N, x, incx, y, incy, func, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  return rc;
}
