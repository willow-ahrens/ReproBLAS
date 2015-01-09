#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "../common/test_perf.h"

#include "../common/test_vecvec_fill_header.h"

const char* vecvec_fill_bench_name(int argc, char** argv);
int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, double scale, double cond, int trials);

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];
  opt_option trials;

  trials.header.type       = opt_int;
  trials.header.short_name = 'a';
  trials.header.long_name  = "trials";
  trials.header.help       = "number of trials";
  trials._int.required     = 0;
  trials._int.min          = 1;
  trials._int.max          = INT_MAX;
  trials._int.value        = 1000;

  if(help._flag.exists){
    opt_show_option(trials);
    vecvec_fill_bench_name(argc, argv);
    return "";
  }

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", vecvec_fill_bench_name(argc, argv), trials);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type, double scale, double cond){
  opt_option trials;

  trials.header.type       = opt_int;
  trials.header.short_name = 'a';
  trials.header.long_name  = "trials";
  trials.header.help       = "number of trials";
  trials._int.required     = 0;
  trials._int.min          = 1;
  trials._int.max          = INT_MAX;
  trials._int.value        = 1000;

  opt_eval_option(argc, argv, &trials);
  int rc = vecvec_fill_bench_test(argc, argv, N, incx, incy, type, scale, cond, trials._int.value);
  return rc;
}
