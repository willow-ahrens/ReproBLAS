#include <indexedBLAS.h>
#include <idxd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../config.h"
#include "../common/test_util.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rdblas1_wrapper.h"

#include "../common/test_file_header.h"

static opt_option func_type;
static opt_option record;
static void validate_external_rdblas1_options_initialize(void){
  func_type._named.header.type       = opt_named;
  func_type._named.header.short_name = 'w';
  func_type._named.header.long_name  = "w_type";
  func_type._named.header.help       = "wrapped function type";
  func_type._named.required          = 1;
  func_type._named.n_names           = wrap_rdblas1_n_names;
  func_type._named.names             = (char**)wrap_rdblas1_names;
  func_type._named.descs             = (char**)wrap_rdblas1_descs;
  func_type._named.value             = wrap_RDSUM;

  record._flag.header.type       = opt_flag;
  record._flag.header.short_name = 'r';
  record._flag.header.long_name  = "record";
  record._flag.header.help       = "record the run insted of testing";
}


int file_show_help(void){
  validate_external_rdblas1_options_initialize();

  opt_show_option(func_type);
  opt_show_option(record);
  return 0;
}

const char* file_name(int argc, char** argv) {
  static char name_buffer[MAX_LINE];

  validate_external_rdblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE, "Validate %s external", wrap_rdblas1_names[func_type._named.value]);
  return name_buffer;
}

int file_test(int argc, char** argv, char *fname) {
  int N;
  char ref_fname[MAX_NAME];
  char Iref_fname[MAX_NAME];

  validate_external_rdblas1_options_initialize();

  opt_eval_option(argc, argv, &func_type);
  opt_eval_option(argc, argv, &record);

  double *X;
  double *Y;

  double ref;
  double_indexed *Iref = idxd_dialloc(DIDEFAULTFOLD);
  idxd_disetzero(DIDEFAULTFOLD, Iref);

  double res;
  double_indexed *Ires = idxd_dialloc(DIDEFAULTFOLD);
  idxd_disetzero(DIDEFAULTFOLD, Ires);

  file_read_vector(fname, &N, (void**)&X, sizeof(double));
  Y = util_dvec_alloc(N, 1);
  //fill Y with 1 where necessary
  util_dvec_fill(N, Y, 1, util_Vec_Constant, 1.0, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_NAME, "%s__%s.dat", fname, wrap_rdblas1_names[func_type._named.value]);
  snprintf(Iref_fname, MAX_NAME, "%s__I%s.dat", fname, wrap_rdblas1_names[func_type._named.value]);

  res = (wrap_rdblas1_func(func_type._named.value))(N, X, 1, Y, 1);
  (wrap_diblas1_func(func_type._named.value))(N, X, 1, Y, 1, Ires);

  if(record._flag.exists){
    free(Iref);
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, Iref, idxd_disize(DIDEFAULTFOLD));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(double*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, idxd_disize(DIDEFAULTFOLD));
    free(Iref);
    Iref = data;
    if(ref != res){
      printf("%s(%s) = %g != %g\n", wrap_rdblas1_names[func_type._named.value], fname, res, ref);
      return 1;
    }
    if(memcmp(&Iref, &Ires, idxd_disize(DIDEFAULTFOLD)) != 0){
      printf("I%s(%s) = %g != %g\n", wrap_rdblas1_names[func_type._named.value], fname, idxd_ddiconv(DIDEFAULTFOLD, Ires), idxd_ddiconv(DIDEFAULTFOLD, Iref));
      printf("Ref I_double:\n");
      idxd_diprint(DIDEFAULTFOLD, Iref);
      printf("\nRes I_double:\n");
      idxd_diprint(DIDEFAULTFOLD, Ires);
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
