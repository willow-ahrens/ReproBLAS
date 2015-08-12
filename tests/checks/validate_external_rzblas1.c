#include <idxdBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../config.h"
#include "../common/test_util.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rzblas1_wrapper.h"

#include "../common/test_file_header.h"

static opt_option func_type;
static opt_option record;
static void validate_external_rzblas1_options_initialize(void){
  func_type._named.header.type       = opt_named;
  func_type._named.header.short_name = 'w';
  func_type._named.header.long_name  = "w_type";
  func_type._named.header.help       = "wrapped function type";
  func_type._named.required          = 1;
  func_type._named.n_names           = wrap_rzblas1_n_names;
  func_type._named.names             = (char**)wrap_rzblas1_names;
  func_type._named.descs             = (char**)wrap_rzblas1_descs;
  func_type._named.value             = wrap_RZSUM;

  record._flag.header.type       = opt_flag;
  record._flag.header.short_name = 'r';
  record._flag.header.long_name  = "record";
  record._flag.header.help       = "record the run insted of testing";
}

int file_show_help(void){
  validate_external_rzblas1_options_initialize();

  opt_show_option(func_type);
  opt_show_option(record);
  return 0;
}

const char* file_name(int argc, char** argv) {
  static char name_buffer[MAX_LINE];

  validate_external_rzblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE, "Validate %s external", wrap_rzblas1_names[func_type._named.value]);
  return name_buffer;
}

int file_test(int argc, char** argv, char *fname) {
  int N;
  char ref_fname[MAX_NAME];
  char Iref_fname[MAX_NAME];

  validate_external_rzblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  opt_eval_option(argc, argv, &record);

  double complex *X;
  double complex *Y;

  double complex ref;
  double_complex_indexed *Iref = idxd_zialloc(DIDEFAULTFOLD);
  idxd_zisetzero(DIDEFAULTFOLD, Iref);

  double complex res;
  double_complex_indexed *Ires = idxd_zialloc(DIDEFAULTFOLD);
  idxd_zisetzero(DIDEFAULTFOLD, Ires);

  file_read_vector(fname, &N, (void**)&X, sizeof(double complex));
  Y = util_zvec_alloc(N, 1);
  //fill Y with -i where necessary
  util_zvec_fill(N, Y, 1, util_Vec_Constant, -_Complex_I, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_NAME, "%s__%s.dat", fname, wrap_rzblas1_names[func_type._named.value]);
  snprintf(Iref_fname, MAX_NAME, "%s__I%s.dat", fname, wrap_rzblas1_names[func_type._named.value]);

  res = (wrap_rzblas1_func(func_type._named.value))(N, X, 1, Y, 1);
  (wrap_ziblas1_func(func_type._named.value))(N, X, 1, Y, 1, Ires);

  if(record._flag.exists){
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, Iref, idxd_zisize(DIDEFAULTFOLD));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(double complex*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, idxd_zisize(DIDEFAULTFOLD));
    free(Iref);
    Iref = data;
    if(ref != res){
      printf("%s(%s) = %g + %gi != %g + %gi\n", wrap_rzblas1_names[func_type._named.value], fname, creal(res), cimag(res), creal(ref), cimag(ref));
      return 1;
    }
    if(memcmp(&Iref, &Ires, sizeof(Iref)) != 0){
      idxd_zziconv_sub(DIDEFAULTFOLD, Ires, &res);
      idxd_zziconv_sub(DIDEFAULTFOLD, Iref, &ref);
      printf("I%s(%s) = %g + %gi != %g + %gi\n", wrap_rzblas1_names[func_type._named.value], fname, creal(res), cimag(res), creal(ref), cimag(ref));
      printf("Ref I_double_Complex:\n");
      idxd_ziprint(DIDEFAULTFOLD, Iref);
      printf("\nRes I_double_Complex:\n");
      idxd_ziprint(DIDEFAULTFOLD, Ires);
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
