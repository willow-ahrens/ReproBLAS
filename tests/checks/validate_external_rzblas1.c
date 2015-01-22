#include <rblas.h>
#include <IndexedFP.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/test_util.h"
#include "../common/test_opt.h"
#include "../common/test_file.h"
#include "rzblas1_wrapper.h"

#include "../common/test_file_header.h"

static opt_option func_type = {._named.header.type       = opt_named,
                               ._named.header.short_name = 'w',
                               ._named.header.long_name  = "w_type",
                               ._named.header.help       = "wrapped function type",
                               ._named.required          = 1,
                               ._named.n_names           = wrap_rzblas1_n_names,
                               ._named.names             = (char**)wrap_rzblas1_names,
                               ._named.descs             = (char**)wrap_rzblas1_descs,
                               ._named.value             = wrap_RZSUM};

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
  snprintf(name_buffer, MAX_LINE, "Validate %s external", wrap_rzblas1_names[func_type._named.value]);
  return name_buffer;
}

int file_test(int argc, char** argv, char *fname) {
  int i;
  int N;
  char ref_fname[MAX_NAME];
  char Iref_fname[MAX_NAME];

  opt_eval_option(argc, argv, &func_type);
  opt_eval_option(argc, argv, &record);

  double complex *x;
  double complex *y;

  double complex ref;
  I_double_Complex Iref;

  double complex res;
  I_double_Complex Ires;

  file_read_vector(fname, &N, (void**)&x, sizeof(double complex));
  y = util_zvec_alloc(N, 1);
  //fill y with -i where necessary
  util_zvec_fill(N, y, 1, util_Vec_Constant, -_Complex_I, 1.0);

  ((char*)file_ext(fname))[0] = '\0';
  snprintf(ref_fname, MAX_NAME, "%s__%s.dat", fname, wrap_rzblas1_names[func_type._named.value]);
  snprintf(Iref_fname, MAX_NAME, "%s__I%s.dat", fname, wrap_rzblas1_names[func_type._named.value]);

  res = (wrap_rzblas1_func(func_type._named.value))(N, x, 1, y, 1);
  Ires = (wrap_Izblas1_func(func_type._named.value))(N, x, 1, y, 1);

  if(record._flag.exists){
    ref = res;
    Iref = Ires;

    file_write_vector(ref_fname, 1, &ref, sizeof(ref));
    file_write_vector(Iref_fname, 1, &Iref, sizeof(Iref));
  } else {
    void *data;
    int unused0;
    file_read_vector(ref_fname, &unused0, &data, sizeof(ref));
    ref = *(double complex*)data;
    free(data);
    file_read_vector(Iref_fname, &unused0, &data, sizeof(Iref));
    Iref = *(I_double_Complex*)data;
    free(data);
    if(ref != res){
      printf("%s(%s) = %g + %gi != %g + %gi\n", wrap_rzblas1_names[func_type._named.value], fname, ZREAL_(res), ZIMAG_(res), ZREAL_(ref), ZIMAG_(ref));
      return 1;
    }
    if(memcmp(&Iref, &Ires, sizeof(Iref)) != 0){
      printf("I%s(%s) = %g + %gi != %g + %gi\n", wrap_rzblas1_names[func_type._named.value], fname, ZREAL_(Iconv2z(Ires)), ZIMAG_(Iconv2z(Ires)), ZREAL_(Iconv2z(Iref)), ZIMAG_(Iconv2z(Iref)));
      printf("Ref I_double_Complex:\n");
      zIprint(Iref);
      printf("\nRes I_double_Complex:\n");
      zIprint(Ires);
      printf("\n");
      return 1;
    }
  }

  free(x);
  free(y);
  return 0;
}
