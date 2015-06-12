#include <stdio.h>
#include "../common/test_opt.h"

#include "../common/test_vecvec_fill_header.h"

int acc_vecvec_fill_show_help(void);
const char* acc_vecvec_fill_name(int argc, char** argv);
int acc_vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY, int trials);

static opt_option trials;

static void acc_vecvec_fill_options_initialize(void){
  trials._int.header.type       = opt_int;
  trials._int.header.short_name = 'a';
  trials._int.header.long_name  = "trials";
  trials._int.header.help       = "number of trials";
  trials._int.required          = 0;
  trials._int.min               = 1;
  trials._int.max               = INT_MAX;
  trials._int.value             = 10000;
}

int vecvec_fill_show_help(void){
  acc_vecvec_fill_options_initialize();

  opt_show_option(trials);
  return acc_vecvec_fill_show_help();
}

const char* vecvec_fill_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  acc_vecvec_fill_options_initialize();

  opt_eval_option(argc, argv, &trials);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s (%d trials)", acc_vecvec_fill_name(argc, argv), trials._int.value);
  return name_buffer;
}

int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY){
  acc_vecvec_fill_options_initialize();

  opt_eval_option(argc, argv, &trials);
  int rc = acc_vecvec_fill_test(argc, argv, N, FillX, RealScaleX, ImagScaleX, incX, FillY, RealScaleY, ImagScaleY, incY, trials._int.value);
  return rc;
}
