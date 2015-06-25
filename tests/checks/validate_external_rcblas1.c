#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../config.h"
#include "../common/test_util.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rcblas1_wrapper.h"

#include "../common/test_file_header.h"

static opt_option func_type;
static opt_option record;
static void validate_external_rcblas1_options_initialize(void){
  func_type._named.header.type       = opt_named;
  func_type._named.header.short_name = 'w';
  func_type._named.header.long_name  = "w_type";
  func_type._named.header.help       = "wrapped function type";
  func_type._named.required          = 1;
  func_type._named.n_names           = wrap_rcblas1_n_names;
  func_type._named.names             = (char**)wrap_rcblas1_names;
  func_type._named.descs             = (char**)wrap_rcblas1_descs;
  func_type._named.value             = wrap_RCSUM;

  record._flag.header.type       = opt_flag;
  record._flag.header.short_name = 'r';
  record._flag.header.long_name  = "record";
  record._flag.header.help       = "record the run insted of testing";
}

int file_show_help(void){
  validate_external_rcblas1_options_initialize();

  opt_show_option(func_type);
  opt_show_option(record);
  return 0;
}

const char* file_name(int argc, char** argv) {
  static char name_buffer[MAX_LINE];

  validate_external_rcblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE, "Validate %s external", wrap_rcblas1_names[func_type._named.value]);
  return name_buffer;
}

int file_test(int argc, char** argv, char *fname) {
  int N;
  char ref_fname[MAX_NAME];
  char Iref_fname[MAX_NAME];

  validate_external_rcblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  opt_eval_option(argc, argv, &record);

  float complex *X;
  float complex *Y;

  float complex ref;
  float_complex_indexed *Iref = cialloc(SIDEFAULTFOLD);
  cisetzero(SIDEFAULTFOLD, Iref);

  float complex res;
  float_complex_indexed *Ires = cialloc(SIDEFAULTFOLD);
  cisetzero(SIDEFAULTFOLD, Ires);

  file_read_vector(fname, &N, (void**)&X, sizeof(float complex));
  Y = util_cvec_alloc(N, 1);
  //fill Y with -i where necessary
  util_cvec_fill(N, Y, 1, util_Vec_Constant, -_Complex_I, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_NAME, "%s__%s.dat", fname, wrap_rcblas1_names[func_type._named.value]);
  snprintf(Iref_fname, MAX_NAME, "%s__I%s.dat", fname, wrap_rcblas1_names[func_type._named.value]);

  res = (wrap_rcblas1_func(func_type._named.value))(N, X, 1, Y, 1);
  (wrap_ciblas1_func(func_type._named.value))(N, X, 1, Y, 1, Ires);

  if(record._flag.exists){
    free(Iref);
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, Iref, cisize(SIDEFAULTFOLD));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(float complex*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, cisize(SIDEFAULTFOLD));
    free(Iref);
    Iref = (float_complex_indexed*)data;
    if(ref != res){
      printf("%s(%s) = %g + %gi != %g + %gi\n", wrap_rcblas1_names[func_type._named.value], fname, crealf(res), cimagf(res), crealf(ref), cimagf(ref));
      return 1;
    }
    if(memcmp(Iref, Ires, cisize(SIDEFAULTFOLD)) != 0){
      cciconv_sub(SIDEFAULTFOLD, Ires, &res);
      cciconv_sub(SIDEFAULTFOLD, Iref, &ref);
      printf("I%s(%s) = %g + %gi != %g + %gi\n", wrap_rcblas1_names[func_type._named.value], fname, crealf(res), cimagf(res), crealf(ref), cimagf(ref));
      printf("Ref I_float_Complex:\n");
      ciprint(SIDEFAULTFOLD, Iref);
      printf("\nRes I_float_Complex:\n");
      ciprint(SIDEFAULTFOLD, Ires);
      printf("\n");
      return 1;
    }
  }

  free(Iref);
  free(Ires);
  free(X);
  free(Y);
  return 0;
}
