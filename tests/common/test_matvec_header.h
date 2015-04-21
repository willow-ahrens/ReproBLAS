#include <stdio.h>
#include "test_opt.h"
#include "test_util.h"

#include "test_header.h"

int matvec_show_help(void);
const char *matvec_name(int argc, char** argv);
int matvec_test(int argc, char** argv, char Order, char TransA, int M, int N, int lda, int incX);

static const char *Order_names[] = {"RowMajor", "ColMajor"};
static const char *Order_descs[] = {"Row Major", "Column Major"};
static opt_option Order  = {._named.header.type       = opt_named,
                            ._named.header.short_name = 'O',
                            ._named.header.long_name  = "Order",
                            ._named.header.help       = "2D ordering",
                            ._named.required          = 0,
                            ._named.n_names           = 2,
                            ._named.names             = (char**)Order_names,
                            ._named.descs             = (char**)Order_descs,
                            ._named.value             = 0};

static const char *TransA_names[] = {"NoTrans", "Trans"};
static const char *TransA_descs[] = {"Don't Transpose", "Transpose"};
static opt_option TransA = {._named.header.type       = opt_named,
                            ._named.header.short_name = 'T',
                            ._named.header.long_name  = "TransA",
                            ._named.header.help       = "transpose A?",
                            ._named.required          = 0,
                            ._named.n_names           = 2,
                            ._named.names             = (char**)TransA_names,
                            ._named.descs             = (char**)TransA_descs,
                            ._named.value             = 0};

static opt_option M      = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'M',
                            ._int.header.long_name  = "M_dim",
                            ._int.header.help       = "M dimension size",
                            ._int.required          = 0,
                            ._int.min               = 0,
                            ._int.max               = INT_MAX,
                            ._int.value             = 2048};

static opt_option N      = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'N',
                            ._int.header.long_name  = "N_dim",
                            ._int.header.help       = "N dimension size",
                            ._int.required          = 0,
                            ._int.min               = 0,
                            ._int.max               = INT_MAX,
                            ._int.value             = 2048};

static opt_option lda    = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'A',
                            ._int.header.long_name  = "lda",
                            ._int.header.help       = "leading A size (0 for auto)",
                            ._int.required          = 0,
                            ._int.min               = 0,
                            ._int.max               = INT_MAX,
                            ._int.value             = 0};

static opt_option incX   = {._int.header.type       = opt_int,
                            ._int.header.short_name = 'x',
                            ._int.header.long_name  = "incX",
                            ._int.header.help       = "X vector increment",
                            ._int.required          = 0,
                            ._int.min               = 1,
                            ._int.max               = INT_MAX,
                            ._int.value             = 1};

int show_help(void){
  opt_show_option(Order);
  opt_show_option(TransA);
  opt_show_option(M);
  opt_show_option(N);
  opt_show_option(lda);
  opt_show_option(incX);
  return matvec_show_help();
}

const char* name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  opt_eval_option(argc, argv, &Order);
  opt_eval_option(argc, argv, &TransA);
  opt_eval_option(argc, argv, &M);
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &lda);
  opt_eval_option(argc, argv, &incX);
  snprintf(name_buffer, MAX_LINE, "%s Order=%c TransA=%c M=%d N=%d lda=%d incX=%d", matvec_name(argc, argv), Order._named.names[Order._named.value][0], TransA._named.names[TransA._named.value][0], M._int.value, N._int.value, lda._int.value, incX._int.value);
  return name_buffer;
}

int test(int argc, char** argv){
  int rc;

  opt_eval_option(argc, argv, &Order);
  opt_eval_option(argc, argv, &TransA);
  opt_eval_option(argc, argv, &M);
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &lda);
  opt_eval_option(argc, argv, &incX);

  int lda_prime;

  switch(Order._named.names[Order._named.value][0]){
    case 'r':
    case 'R':
      if(lda._int.value == 0){
        lda._int.value = N._int.value;
      }
      if(N._int.value > lda._int.value){
        fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent N=%d, lda=%d\n", N._int.value, lda._int.value);
        return 125;
      }
      break;
    case 'c':
    case 'C':
      if(lda._int.value == 0){
        lda._int.value = M._int.value;
      }
      if(M._int.value > lda._int.value){
        fprintf(stderr, "ReproBLAS error: column major matrix arguments inconsistent M=%d, lda=%d\n", M._int.value, lda._int.value);
        return 125;
      }
      break;
  }
  rc = matvec_test(argc, argv, Order._named.names[Order._named.value][0], TransA._named.names[TransA._named.value][0], M._int.value, N._int.value, lda._int.value, incX._int.value);
  return rc;
}
