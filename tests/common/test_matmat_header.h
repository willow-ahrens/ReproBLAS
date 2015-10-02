#include <stdio.h>
#include "test_opt.h"
#include "test_util.h"

#include "test_header.h"

int matmat_show_help(void);
const char *matmat_name(int argc, char** argv);
int matmat_test(int argc, char** argv, char Order, char TransA, char TransB, int M, int N, int K, int lda, int ldb, int ldc);

static const char *Order_names[] = {"RowMajor", "ColMajor"};
static const char *Order_descs[] = {"Row Major", "Column Major"};
static opt_option Order;
static const char *TransA_names[] = {"NoTrans", "Trans", "ConjTrans"};
static const char *TransA_descs[] = {"Don't Transpose", "Transpose", "Conjugate Transpose"};
static opt_option TransA;
static const char *TransB_names[] = {"NoTrans", "Trans", "ConjTrans"};
static const char *TransB_descs[] = {"Don't Transpose", "Transpose", "Conjugate Transpose"};
static opt_option TransB;
static opt_option M;
static opt_option N;
static opt_option K;
static opt_option lda;
static opt_option ldb;
static opt_option ldc;

static void matmat_options_initialize(void){
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

  TransB._named.header.type       = opt_named;
  TransB._named.header.short_name = 'U';
  TransB._named.header.long_name  = "TransB";
  TransB._named.header.help       = "transpose B?";
  TransB._named.required          = 0;
  TransB._named.n_names           = 3;
  TransB._named.names             = (char**)TransB_names;
  TransB._named.descs             = (char**)TransB_descs;
  TransB._named.value             = 0;

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

  K._int.header.type       = opt_int;
  K._int.header.short_name = 'K';
  K._int.header.long_name  = "K_dim";
  K._int.header.help       = "K dimension size";
  K._int.required          = 0;
  K._int.min               = 0;
  K._int.max               = INT_MAX;
  K._int.value             = 2048;

  lda._int.header.type       = opt_int;
  lda._int.header.short_name = '\0';
  lda._int.header.long_name  = "lda";
  lda._int.header.help       = "leading A size (0 for auto)";
  lda._int.required          = 0;
  lda._int.min               = INT_MIN;
  lda._int.max               = INT_MAX;
  lda._int.value             = 0;

  ldb._int.header.type       = opt_int;
  ldb._int.header.short_name = '\0';
  ldb._int.header.long_name  = "ldb";
  ldb._int.header.help       = "leading B size (0 for auto)";
  ldb._int.required          = 0;
  ldb._int.min               = INT_MIN;
  ldb._int.max               = INT_MAX;
  ldb._int.value             = 0;

  ldc._int.header.type       = opt_int;
  ldc._int.header.short_name = '\0';
  ldc._int.header.long_name  = "ldc";
  ldc._int.header.help       = "leading c size (0 for auto)";
  ldc._int.required          = 0;
  ldc._int.min               = INT_MIN;
  ldc._int.max               = INT_MAX;
  ldc._int.value             = 0;
}

int show_help(void){
  matmat_options_initialize();

  opt_show_option(Order);
  opt_show_option(TransA);
  opt_show_option(TransB);
  opt_show_option(M);
  opt_show_option(N);
  opt_show_option(K);
  opt_show_option(lda);
  opt_show_option(ldb);
  opt_show_option(ldc);
  return matmat_show_help();
}

const char* name(int argc, char** argv){
  static char name_buffer[MAX_LINE];

  matmat_options_initialize();

  opt_eval_option(argc, argv, &Order);
  opt_eval_option(argc, argv, &TransA);
  opt_eval_option(argc, argv, &TransB);
  opt_eval_option(argc, argv, &M);
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &K);
  opt_eval_option(argc, argv, &lda);
  opt_eval_option(argc, argv, &ldb);
  opt_eval_option(argc, argv, &ldc);

  snprintf(name_buffer, MAX_LINE, "%s Order=%c TransA=%c TransB=%c M=%d N=%d K=%d lda=%d ldb=%d ldc=%d", matmat_name(argc, argv), Order._named.names[Order._named.value][0], TransA._named.names[TransA._named.value][0], TransB._named.names[TransB._named.value][0], M._int.value, N._int.value, K._int.value, lda._int.value, ldb._int.value, ldc._int.value);
  return name_buffer;
}

int test(int argc, char** argv){
  int rc;

  matmat_options_initialize();

  opt_eval_option(argc, argv, &Order);
  opt_eval_option(argc, argv, &TransA);
  opt_eval_option(argc, argv, &TransB);
  opt_eval_option(argc, argv, &M);
  opt_eval_option(argc, argv, &N);
  opt_eval_option(argc, argv, &K);
  opt_eval_option(argc, argv, &lda);
  opt_eval_option(argc, argv, &ldb);
  opt_eval_option(argc, argv, &ldc);

  switch(Order._named.names[Order._named.value][0]){
    case 'r':
    case 'R':
      switch(TransA._named.names[TransA._named.value][0]){
        case 'n':
        case 'N':
          if(lda._int.value < 0){
            lda._int.value = K._int.value - lda._int.value;
          }else if(lda._int.value == 0){
            lda._int.value = K._int.value;
          }else if(K._int.value > lda._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent K=%d, lda=%d\n", K._int.value, lda._int.value);
            return 125;
          }
          break;
        default:
          if(lda._int.value < 0){
            lda._int.value = M._int.value - lda._int.value;
          }else if(lda._int.value == 0){
            lda._int.value = M._int.value;
          }else if(M._int.value > lda._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent M=%d, lda=%d\n", M._int.value, lda._int.value);
            return 125;
          }
          break;
      }
      switch(TransB._named.names[TransB._named.value][0]){
        case 'n':
        case 'N':
          if(ldb._int.value < 0){
            ldb._int.value = N._int.value - ldb._int.value;
          }else if(ldb._int.value == 0){
            ldb._int.value = N._int.value;
          }else if(N._int.value > ldb._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent N=%d, ldb=%d\n", N._int.value, ldb._int.value);
            return 125;
          }
          break;
        default:
          if(ldb._int.value < 0){
            ldb._int.value = K._int.value - ldb._int.value;
          }else if(ldb._int.value == 0){
            ldb._int.value = K._int.value;
          }else if(K._int.value > ldb._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent K=%d, ldb=%d\n", K._int.value, ldb._int.value);
            return 125;
          }
          break;
      }
      if(ldc._int.value < 0){
        ldc._int.value = N._int.value - ldc._int.value;
      }else if(ldc._int.value == 0){
        ldc._int.value = N._int.value;
      }else if(N._int.value > ldc._int.value){
        fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent N=%d, ldc=%d\n", N._int.value, ldc._int.value);
        return 125;
      }
      break;
    default:
      switch(TransA._named.names[TransA._named.value][0]){
        case 'n':
        case 'N':
          if(lda._int.value < 0){
            lda._int.value = M._int.value - lda._int.value;
          }else if(lda._int.value == 0){
            lda._int.value = M._int.value;
          }else if(M._int.value > lda._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent M=%d, lda=%d\n", M._int.value, lda._int.value);
            return 125;
          }
          break;
        default:
          if(lda._int.value < 0){
            lda._int.value = K._int.value - lda._int.value;
          }else if(lda._int.value == 0){
            lda._int.value = K._int.value;
          }else if(K._int.value > lda._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent K=%d, lda=%d\n", K._int.value, lda._int.value);
            return 125;
          }
          break;
      }
      switch(TransB._named.names[TransB._named.value][0]){
        case 'n':
        case 'N':
          if(ldb._int.value < 0){
            ldb._int.value = K._int.value - ldb._int.value;
          }else if(ldb._int.value == 0){
            ldb._int.value = K._int.value;
          }else if(K._int.value > ldb._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent K=%d, ldb=%d\n", K._int.value, ldb._int.value);
            return 125;
          }
          break;
        default:
          if(ldb._int.value < 0){
            ldb._int.value = N._int.value - ldb._int.value;
          }else if(ldb._int.value == 0){
            ldb._int.value = N._int.value;
          }else if(N._int.value > ldb._int.value){
            fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent N=%d, ldb=%d\n", N._int.value, ldb._int.value);
            return 125;
          }
          break;
      }
      if(ldc._int.value < 0){
        ldc._int.value = M._int.value - ldc._int.value;
      }else if(ldc._int.value == 0){
        ldc._int.value = M._int.value;
      }else if(M._int.value > ldc._int.value){
        fprintf(stderr, "ReproBLAS error: row major matrix arguments inconsistent M=%d, ldc=%d\n", M._int.value, ldc._int.value);
        return 125;
      }
      break;
  }
  rc = matmat_test(argc, argv, Order._named.names[Order._named.value][0], TransA._named.names[TransA._named.value][0], TransB._named.names[TransB._named.value][0], M._int.value, N._int.value, K._int.value, lda._int.value, ldb._int.value, ldc._int.value);
  return rc;
}
