#include <indexedBLAS.h>
#include <indexed.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/test_util.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rdblas1_wrapper.h"

#include "../common/test_file_header.h"

static opt_option func_type = {._named.header.type       = opt_named,
                               ._named.header.short_name = 'w',
                               ._named.header.long_name  = "w_type",
                               ._named.header.help       = "wrapped function type",
                               ._named.required          = 1,
                               ._named.n_names           = wrap_rdblas1_n_names,
                               ._named.names             = (char**)wrap_rdblas1_names,
                               ._named.descs             = (char**)wrap_rdblas1_descs,
                               ._named.value             = wrap_RDSUM};

static opt_option record    = {._flag.header.type       = opt_flag,
                               ._flag.header.short_name = 'r',
                               ._flag.header.long_name  = "record",
                               ._flag.header.help       = "record the run insted of testing"};

int file_show_help(void){
  opt_show_option(func_type);
  opt_show_option(record);
  return 0;
}

const char* file_name(int argc, char** argv) {
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &func_type);
  snprintf(name_buffer, MAX_LINE, "Validate %s external", wrap_rdblas1_names[func_type._named.value]);
  return name_buffer;
}

int file_test(int argc, char** argv, char *fname) {
  int N;
  char ref_fname[MAX_NAME];
  char Iref_fname[MAX_NAME];

  opt_eval_option(argc, argv, &func_type);
  opt_eval_option(argc, argv, &record);

  double *X;
  double *Y;

  double ref;
  I_double Iref;

  double res;
  I_double Ires;

  file_read_vector(fname, &N, (void**)&X, sizeof(double));
  Y = util_dvec_alloc(N, 1);
  //fill Y with 1 where necessary
  util_dvec_fill(N, Y, 1, util_Vec_Constant, 1.0, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_NAME, "%s__%s.dat", fname, wrap_rdblas1_names[func_type._named.value]);
  snprintf(Iref_fname, MAX_NAME, "%s__I%s.dat", fname, wrap_rdblas1_names[func_type._named.value]);

  res = (wrap_rdblas1_func(func_type._named.value))(N, X, 1, Y, 1);
  Ires = (wrap_Idblas1_func(func_type._named.value))(N, X, 1, Y, 1);

  if(record._flag.exists){
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, &Iref, sizeof(Iref));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(double*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, sizeof(Iref));
    Iref = *(I_double*)data;
    free(data);
    if(ref != res){
      printf("%s(%s) = %g != %g\n", wrap_rdblas1_names[func_type._named.value], fname, res, ref);
      return 1;
    }
    if(memcmp(&Iref, &Ires, sizeof(Iref)) != 0){
      printf("I%s(%s) = %g != %g\n", wrap_rdblas1_names[func_type._named.value], fname, ddiconv(&Ires, DEFAULT_FOLD), ddiconv(&Iref, DEFAULT_FOLD));
      printf("Ref I_double:\n");
      dIprint(Iref);
      printf("\nRes I_double:\n");
      dIprint(Ires);
      printf("\n");
      return 1;
    }
  }

  free(X);
  free(Y);
  return 0;
}
