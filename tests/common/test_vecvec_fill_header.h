#include <stdio.h>
#include "test_opt.h"

#include "test_vecvec_header.h"

int vecvec_fill_show_help(void);
const char* vecvec_fill_name(int argc, char** argv);
int vecvec_fill_test(int argc, char** argv, int N, int FillX, double RealScaleX, double ImagScaleX, int incX, int FillY, double RealScaleY, double ImagScaleY, int incY);

static opt_option FillX;
static opt_option RealScaleX;
static opt_option ImagScaleX;
static opt_option FillY;
static opt_option RealScaleY;
static opt_option ImagScaleY;

static void vecvec_fill_options_initialize(void){
  FillX._named.header.type       = opt_named;
  FillX._named.header.short_name = '\0';
  FillX._named.header.long_name  = "FillX";
  FillX._named.header.help       = "X fill type";
  FillX._named.required          = 0;
  FillX._named.n_names           = (int)util_vec_fill_n_names;
  FillX._named.names             = (char**)util_vec_fill_names;
  FillX._named.descs             = (char**)util_vec_fill_descs;
  FillX._named.value             = 0;

  RealScaleX._double.header.type       = opt_double;
  RealScaleX._double.header.short_name = '\0';
  RealScaleX._double.header.long_name  = "RealScaleX";
  RealScaleX._double.header.help       = "X scale (real)";
  RealScaleX._double.required          = 0;
  RealScaleX._double.min               = -1 * DBL_MAX;
  RealScaleX._double.max               = DBL_MAX;
  RealScaleX._double.value             = 1.0;

  ImagScaleX._double.header.type       = opt_double;
  ImagScaleX._double.header.short_name = '\0';
  ImagScaleX._double.header.long_name  = "ImagScaleX";
  ImagScaleX._double.header.help       = "X scale (imaginary)";
  ImagScaleX._double.required          = 0;
  ImagScaleX._double.min               = -1 * DBL_MAX;
  ImagScaleX._double.max               = DBL_MAX;
  ImagScaleX._double.value             = 0.0;

  FillY._named.header.type       = opt_named;
  FillY._named.header.short_name = '\0';
  FillY._named.header.long_name  = "FillY";
  FillY._named.header.help       = "Y fill type";
  FillY._named.required          = 0;
  FillY._named.n_names           = (int)util_vec_fill_n_names;
  FillY._named.names             = (char**)util_vec_fill_names;
  FillY._named.descs             = (char**)util_vec_fill_descs;
  FillY._named.value             = 0;

  RealScaleY._double.header.type       = opt_double;
  RealScaleY._double.header.short_name = '\0';
  RealScaleY._double.header.long_name  = "RealScaleY";
  RealScaleY._double.header.help       = "Y scale (real)";
  RealScaleY._double.required          = 0;
  RealScaleY._double.min               = -1 * DBL_MAX;
  RealScaleY._double.max               = DBL_MAX;
  RealScaleY._double.value             = 1.0;

  ImagScaleY._double.header.type       = opt_double;
  ImagScaleY._double.header.short_name = '\0';
  ImagScaleY._double.header.long_name  = "ImagScaleY";
  ImagScaleY._double.header.help       = "Y scale (imaginary)";
  ImagScaleY._double.required          = 0;
  ImagScaleY._double.min               = -1 * DBL_MAX;
  ImagScaleY._double.max               = DBL_MAX;
  ImagScaleY._double.value             = 0.0;
}

int vecvec_show_help(void){
  vecvec_fill_options_initialize();

  opt_show_option(FillX);
  opt_show_option(RealScaleX);
  opt_show_option(ImagScaleX);
  opt_show_option(FillY);
  opt_show_option(RealScaleY);
  opt_show_option(ImagScaleY);
  return vecvec_fill_show_help();
}

const char* vecvec_name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  vecvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &RealScaleX);
  opt_eval_option(argc, argv, &ImagScaleX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &RealScaleY);
  opt_eval_option(argc, argv, &ImagScaleY);
  snprintf(name_buffer, MAX_LINE * sizeof(char), "%s X=%s*(%g+%gi) Y=%s*(%g+%gi)", vecvec_fill_name(argc, argv), util_vec_fill_names[FillX._named.value], RealScaleX._double.value, ImagScaleX._double.value, util_vec_fill_names[FillY._named.value], RealScaleY._double.value, ImagScaleY._double.value);
  return name_buffer;
}

int vecvec_test(int argc, char** argv, int N, int incX, int incY){
  int rc;

  vecvec_fill_options_initialize();

  opt_eval_option(argc, argv, &FillX);
  opt_eval_option(argc, argv, &RealScaleX);
  opt_eval_option(argc, argv, &ImagScaleX);
  opt_eval_option(argc, argv, &FillY);
  opt_eval_option(argc, argv, &RealScaleY);
  opt_eval_option(argc, argv, &ImagScaleY);
  rc = vecvec_fill_test(argc, argv, N, FillX._named.value, RealScaleX._double.value, ImagScaleX._double.value, incX, FillY._named.value, RealScaleY._double.value, ImagScaleY._double.value, incY);
  return rc;
}
