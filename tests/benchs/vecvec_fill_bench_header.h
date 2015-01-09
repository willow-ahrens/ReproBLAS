#include <stdio.h>
#include "../common/test_opt.h"
#include "../common/test_vec.h"
#include "../common/test_perf.h"

#include "../common/test_vecvec_fill_header.h"

int vecvec_fill_bench_desc(void);
int vecvec_fill_bench_show_help(void);
const char* vecvec_fill_bench_name(int argc, char** argv);
int vecvec_fill_bench_test(int argc, char** argv, int N, int incx, int incy, int type, double scale, double cond, int trials);

static opt_option trials = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'a',
                            ._int.header.long_name  = "trials",
                            ._int.header.help       = "number of trials",
                            ._int.required          = 0,
                            ._int.min               = 1,
                            ._int.max               = INT_MAX,
                            ._int.value             = 1000};

static opt_option desc   = {._flag.header.type       = opt_int,
                            ._flag.header.short_name = 'd',
                            ._flag.header.long_name  = "desc",
                            ._flag.header.help       = "show benchmark description"};

int vecvec_fill_show_help(void){
  opt_show_option(trials);
  return vecvec_fill_bench_show_help();
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", vecvec_fill_bench_name(argc, argv), trials._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int incx, int incy, int type, double scale, double cond){
  opt_eval_option(argc, argv, &desc);
  if(desc._flag.exists){
    return vecvec_fill_bench_desc();
  }

  opt_eval_option(argc, argv, &trials);
  int rc = vecvec_fill_bench_test(argc, argv, N, incx, incy, type, scale, cond, trials._int.value);
  return rc;
}
