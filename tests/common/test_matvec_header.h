#include <stdio.h>
#include "test_opt.h"
#include "test_util.h"

#include "test_header.h"

int matvec_show_help(void);
const char *matvec_name(int argc, char** argv);
int matvec_test(int argc, char** argv, char Order, char TransA, int M, int N, int lda, int incX, int incY);

static const char *Order_names[] = {"RowMajor", "ColMajor"};
static const char *Order_descs[] = {"Row Major", "Column Major"};
static opt_option Order;
static const char *TransA_names[] = {"NoTrans", "Trans", "ConjTrans"};
static const char *TransA_descs[] = {"Don't Transpose", "Transpose", "Conjugate Transpose"};
static opt_option TransA;
static opt_option M;
static opt_option N;
static opt_option lda;
static opt_option incX;
static opt_option incY;

static void matvec_options_initialize(void){
  Order._named.header.type       = opt_named;
  Order._named.header.short_name = 'O';
  Order._named.header.long_name  = "Order";
  Order._named.header.help       = "2D ordering";
  Order._named.required          = 0;
  Order._named.n_names           = 2;
  Order._named.names             = (char**)Order_names;
  Order._named.descs             = (char**)Order_descs;
  Order._named.value             = 0;

  TransA._named.header.type       = opt_named;
  TransA._named.header.short_name = 'T';
  TransA._named.header.long_name  = "TransA";
  TransA._named.header.help       = "transpose A?";
  TransA._named.required          = 0;
  TransA._named.n_names           = 3;
  TransA._named.names             = (char**)TransA_names;
  TransA._named.descs             = (char**)TransA_descs;
  TransA._named.value             = 0;

  M._int.header.type       = opt_int;
  M._int.header.short_name = 'M';
  M._int.header.long_name  = "M_dim";
  M._int.header.help       = "M dimension size";
  M._int.required          = 0;
  M._int.min               = 0;
  M._int.max               = INT_MAX;
  M._int.value             = 2048;

  N._int.header.type       = opt_int;
  N._int.header.short_name = 'N';
  N._int.header.long_name  = "N_dim";
  N._int.header.help       = "N dimension size";
  N._int.required          = 0;
  N._int.min               = 0;
  N._int.max               = INT_MAX;
  N._int.value             = 2048;

  lda._int.header.type       = opt_int;
  lda._int.header.short_name = '\0';
  lda._int.header.long_name  = "lda";
  lda._int.header.help       = "leading A size (0 for auto)";
  lda._int.required          = 0;
  lda._int.min               = INT_MIN;
  lda._int.max               = INT_MAX;
  lda._int.value             = 0;

  incX._int.header.type       = opt_int;
  incX._int.header.short_name = '\0';
  incX._int.header.long_name  = "incX";
  incX._int.header.help       = "X vector increment";
  incX._int.required          = 0;
  incX._int.min               = 1;
  incX._int.max               = INT_MAX;
  incX._int.value             = 1;

  incY._int.header.type       = opt_int;
  incY._int.header.short_name = '\0';
  incY._int.header.long_name  = "incY";
  incY._int.header.help       = "Y vector increment";
  incY._int.required          = 0;
  incY._int.min               = 1;
  incY._int.max               = INT_MAX;
  incY._int.value             = 1;
}

int show_help(void){
  matvec_options_initialize();

  opt_show_option(Order);
  opt_show_option(TransA);
  opt_show_option(M);
  opt_show_option(N);
  opt_show_option(lda);
  opt_show_option(incX);
  opt_show_option(incY);
  return matvec_show_help();
}

const char* name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  matvec_options_initialize();

  opt_eval_option(argc, argv, &Order);
  opt_eval_option(argc, argv, &TransA);
  opt_eval_option(argc, argv, &M);
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &lda);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);
  snprintf(name_buffer, MAX_LINE, "%s Order=%c TransA=%c M=%d N=%d lda=%d incX=%d incY=%d", matvec_name(argc, argv), Order._named.names[Order._named.value][0], TransA._named.names[TransA._named.value][0], M._int.value, N._int.value, lda._int.value, incX._int.value, incY._int.value);
  return name_buffer;
}

int test(int argc, char** argv){
  int rc;

  matvec_options_initialize();

  opt_eval_option(argc, argv, &Order);
  opt_eval_option(argc, argv, &TransA);
  opt_eval_option(argc, argv, &M);
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &lda);
  opt_eval_option(argc, argv, &incX);
  opt_eval_option(argc, argv, &incY);

  switch(Order._named.names[Order._named.value][0]){
    case 'r':
    case 'R':
      if(lda._int.value < 0){
        lda._int.value = N._int.value - lda._int.value;
      }else if(lda._int.value == 0){
        lda._int.value = N._int.value;
      }else if(N._int.value > lda._int.value){
        fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent N=%d, lda=%d\n", N._int.value, lda._int.value);
        return 125;
      }
      break;
    default:
      if(lda._int.value < 0){
        lda._int.value = M._int.value - lda._int.value;
      }else if(lda._int.value == 0){
        lda._int.value = M._int.value;
      }else if(M._int.value > lda._int.value){
        fprintf(stderr, "ReproBLAS error: column major matrix arguments inconsistent M=%d, lda=%d\n", M._int.value, lda._int.value);
        return 125;
      }
      break;
  }
  rc = matvec_test(argc, argv, Order._named.names[Order._named.value][0], TransA._named.names[TransA._named.value][0], M._int.value, N._int.value, lda._int.value, incX._int.value, incY._int.value);
  return rc;
}
