#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/test_util.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rcblas1_wrapper.h"

#include "../common/test_file_header.h"

static opt_option func_type = {._named.header.type       = opt_named,
                               ._named.header.short_name = 'w',
                               ._named.header.long_name  = "w_type",
                               ._named.header.help       = "wrapped function type",
                               ._named.required          = 1,
                               ._named.n_names           = wrap_rcblas1_n_names,
                               ._named.names             = (char**)wrap_rcblas1_names,
                               ._named.descs             = (char**)wrap_rcblas1_descs,
                               ._named.value             = wrap_RCSUM};

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
  snprintf(name_buffer, MAX_LINE, "Validate %s external", wrap_rcblas1_names[func_type._named.value]);
  return name_buffer;
}

int file_test(int argc, char** argv, char *fname) {
  int i;
  int N;
  char ref_fname[MAX_NAME];
  char Iref_fname[MAX_NAME];

  opt_eval_option(argc, argv, &func_type);
  opt_eval_option(argc, argv, &record);

  float complex *x;
  float complex *y;

  float complex ref;
  I_float_Complex Iref;

  float complex res;
  I_float_Complex Ires;

  file_read_vector(fname, &N, (void**)&x, sizeof(float complex));
  y = cvec_alloc(N, 1);
  //fill y with -i where necessary
  cvec_fill(N, y, 1, vec_fill_CONSTANT, -_Complex_I, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_NAME, "%s__%s.dat", fname, wrap_rcblas1_names[func_type._named.value]);
  snprintf(Iref_fname, MAX_NAME, "%s__I%s.dat", fname, wrap_rcblas1_names[func_type._named.value]);

  res = (wrap_rcblas1_func(func_type._named.value))(N, x, 1, y, 1);
  Ires = (wrap_Icblas1_func(func_type._named.value))(N, x, 1, y, 1);

  if(record._flag.exists){
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, &Iref, sizeof(Iref));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(float complex*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, sizeof(Iref));
    Iref = *(I_float_Complex*)data;
    free(data);
    if(ref != res){
      printf("%s(%s) = %g + %gi != %g + %gi\n", wrap_rcblas1_names[func_type._named.value], fname, CREAL_(res), CIMAG_(res), CREAL_(ref), CIMAG_(ref));
      return 1;
    }
    if(memcmp(&Iref, &Ires, sizeof(Iref)) != 0){
      printf("I%s(%s) = %g + %gi != %g + %gi\n", wrap_rcblas1_names[func_type._named.value], fname, CREAL_(Iconv2c(Ires)), CIMAG_(Iconv2c(Ires)), CREAL_(Iconv2c(Iref)), CIMAG_(Iconv2c(Iref)));
      printf("Ref I_float_Complex:\n");
      cIprint(Iref);
      printf("\nRes I_float_Complex:\n");
      cIprint(Ires);
      printf("\n");
      return 1;
    }
  }

  free(x);
  free(y);
  return 0;
}
