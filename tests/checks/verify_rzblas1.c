#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/test_opt.h"
#include "rzblas1_wrapper.h"

#include "../common/test_vecvec_fill_header.h"

static opt_option func_type = {._named.header.type       = opt_named,
                               ._named.header.short_name = 'w',
                               ._named.header.long_name  = "w_type",
                               ._named.header.help       = "wrapped function type",
                               ._named.required          = 1,
                               ._named.n_names           = wrap_rzblas1_n_names,
                               ._named.names             = (char**)wrap_rzblas1_names,
                               ._named.descs             = (char**)wrap_rzblas1_descs,
                               ._named.value             = wrap_RZSUM};

int verify_rzblas1_reproducibility(int N, double complex* x, int incX, double complex* y, int incY, int func, double complex ref, I_double_Complex Iref, int max_num_blocks) {
  // GENERATE DATA
  int i, j;
  double complex res;
  I_double_Complex Ires;
  int num_blocks = 1;

  int block_N = (N + num_blocks - 1) / num_blocks;

  num_blocks = 1;
  while (num_blocks < N && num_blocks <= max_num_blocks) {
    if (num_blocks == 1)
      res = (wrap_rzblas1_func(func))(N, x, incX, y, incY);
    else {
      block_N =  (N + num_blocks - 1) / num_blocks;
      zISetZero(Ires);
      for (j = 0; j < N; j += block_N) {
        block_N = block_N < N - j ? block_N : (N-j);
        zIAdd(&Ires, (wrap_Izblas1_func(func))(block_N, x + j * incX, incX, y + j * incY, incY));
      }
      res = Iconv2z(Ires);
    }
    if (res != ref) {
      printf("%s(x, y)[num_blocks=%d,block_N=%d] = %g + %gi != %g + %gi\n", wrap_rzblas1_names[func], num_blocks, block_N, CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
      if (num_blocks == 1) {
        Ires = (wrap_Izblas1_func(func))(N, x, incX, y, incY);
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

int vecvec_fill_show_help(void){
  opt_show_option(func_type);
  return 0;
}

extern const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "Verify %s reproducibility", wrap_rzblas1_names[func_type._named.value]);
  return name_buffer;
}

extern int vecvec_fill_test(int argc, char** argv, int N, int incX, int incY, int type, double scale, double cond){
  int rc = 0;
  double complex ref;
  I_double_Complex Iref;
  int max_num_blocks = 1024;

  util_random_seed();

  double complex *x = util_zvec_alloc(N, incX);
  double complex *y = util_zvec_alloc(N, incY);

  opt_eval_option(argc, argv, &func_type);

  //fill x
  util_zvec_fill(N, x, incX, type, scale, cond);

  //fill y with -i where necessary
  util_zvec_fill(N, y, incY, util_Vec_Constant, -_Complex_I, 1.0);

  //nrm2 doesn't make sense with more than 1 block.
  if(func_type._named.value == wrap_RDZNRM2){
    max_num_blocks = 1;
  }

  //compute with unpermuted data
  ref  = (wrap_rzblas1_func(func_type._named.value))(N, x, incX, y, incY);
  Iref = (wrap_Izblas1_func(func_type._named.value))(N, x, incX, y, incY);

  util_zvec_reverse(N, x, incX, NULL, 1);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_sort(N, x, incX, NULL, 1, util_Increasing);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_sort(N, x, incX, NULL, 1, util_Decreasing);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_sort(N, x, incX, NULL, 1, util_Increasing_Magnitude);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_sort(N, x, incX, NULL, 1, util_Decreasing_Magnitude);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_shuffle(N, x, incX, NULL, 1);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_shuffle(N, x, incX, NULL, 1);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_shuffle(N, x, incX, NULL, 1);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  util_zvec_shuffle(N, x, incX, NULL, 1);

  rc = verify_rzblas1_reproducibility(N, x, incX, y, incY, func_type._named.value, ref, Iref, max_num_blocks);
  if(rc != 0){
    return rc;
  }

  return rc;
}
