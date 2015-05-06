#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "test_util.h"

const int util_vec_fill_n_names  = 15;
const char *util_vec_fill_names[] = {"constant",
                                     "rand",
                                     "2*rand-1",
                                     "rand+(rand-1)",
                                     "normal",
                                     "sine",
                                     "small+grow*big",
                                     "small+rand*big",
                                     "rand_cond",
                                     "constant[drop]",
                                     "rand[drop]",
                                     "2*rand-1[drop]",
                                     "rand+(rand-1)[drop]",
                                     "normal[drop]",
                                     "sine[drop]"};
const char *util_vec_fill_descs[] = {"Constant",
                                     "Random",
                                     "2*Random-1",
                                     "Random+(Random-1)",
                                     "Normal",
                                     "Sine(2pi*(i/n))",
                                     "Small+(i/n)*Big",
                                     "Small+Rand*Big",
                                     "RandomConditioned",
                                     "Constant[drop]",
                                     "Random[drop]",
                                     "2*Random-1[drop]",
                                     "Random+(Random-1)[drop]",
                                     "Normal[drop]",
                                     "Sine(2pi*(i/n))[drop]"};

const int  util_mat_fill_n_names  = 16;
const char *util_mat_fill_names[] = {"constant",
                                     "rand",
                                     "2*rand-1",
                                     "rand+(rand-1)",
                                     "normal",
                                     "sine",
                                     "small+grow*big",
                                     "small+rand*big",
                                     "rand_cond",
                                     "constant[drop]",
                                     "rand[drop]",
                                     "2*rand-1[drop]",
                                     "rand+(rand-1)[drop]",
                                     "normal[drop]",
                                     "sine[drop]",
                                     "identity"};
const char *util_mat_fill_descs[] = {"Constant",
                                     "Random",
                                     "2*Random-1",
                                     "Random+(Random-1)",
                                     "Normal",
                                     "Sine(2pi*(i/n))",
                                     "Small+(i/n)*Big",
                                     "Small+Rand*Big",
                                     "RandomConditioned",
                                     "Constant[drop]",
                                     "Random[drop]",
                                     "2*Random-1[drop]",
                                     "Random+(Random-1)[drop]",
                                     "Normal[drop]",
                                     "Sine(2pi*(i/n))[drop]",
                                     "Identity"};

#define PI 3.14159265358979323846

void util_random_seed(void) {
  struct timeval st;
  gettimeofday( &st, NULL );
  srand((long)(st.tv_usec + 1e6*st.tv_sec));
}

//TODO add random in range support (or find an external library to generate input data)
static double util_drand48(){
  unsigned long l = 0;
  int i;
  int r;
  for(i = 0; i < 4; i++){
    do{
      r = rand();
    }while(r >= (RAND_MAX/256)*256);
    l <<= 8;
    l += r % 256;
  }
  double ret =  ((double)l/1.0l);
  printf("%g\n", ret);
  return ret;
}

typedef int (*compare_func)(int i, int j, void *data);
typedef void (*swap_func)(int i, int j, void *data);
typedef int (*elem_compare_func)(void *a, void *b, util_comp_t comp);

int util_dcompare(void *a, void *b, util_comp_t comp){
  double a_prime;
  double b_prime;
  switch(comp){
    case util_Increasing:
      if(*((double*)a) < *((double*)b)){
        return -1;
      }else if(*((double*)a) > *((double*)b)){
        return 1;
      }else{
        return 0;
      }
    case util_Decreasing:
      return -1 * util_dcompare(a, b, util_Increasing);
    case util_Increasing_Magnitude:
      a_prime = fabs(*((double*)a));
      b_prime = fabs(*((double*)b));
      return util_dcompare(&a_prime, &b_prime, util_Increasing);
    case util_Decreasing_Magnitude:
      return -1 * util_dcompare(a, b, util_Increasing_Magnitude);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

int util_scompare(void *a, void *b, util_comp_t comp){
  float a_prime;
  float b_prime;
  switch(comp){
    case util_Increasing:
      if(*((float*)a) < *((float*)b)){
        return -1;
      }else if(*((float*)a) > *((float*)b)){
        return 1;
      }else{
        return 0;
      }
    case util_Decreasing:
      return -1 * util_scompare(a, b, util_Increasing);
    case util_Increasing_Magnitude:
      a_prime = fabs(*((float*)a));
      b_prime = fabs(*((float*)b));
      return util_scompare(&a_prime, &b_prime, util_Increasing);
    case util_Decreasing_Magnitude:
      return -1 * util_scompare(a, b, util_Increasing_Magnitude);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

int util_zcompare(void *a, void *b, util_comp_t comp){
  double a_prime;
  double b_prime;
  switch(comp){
    case util_Increasing:
      a_prime = ZREAL_(*((double complex*)a));
      b_prime = ZREAL_(*((double complex*)b));
      return util_dcompare(&a_prime, &b_prime, util_Increasing);
    case util_Decreasing:
      return -1 * util_zcompare(a, b, util_Increasing);
    case util_Increasing_Magnitude:
      a_prime = cabs(*((double complex*)a));
      b_prime = cabs(*((double complex*)b));
      return util_dcompare(&a_prime, &b_prime, util_Increasing);
    case util_Decreasing_Magnitude:
      return -1 * util_zcompare(a, b, util_Increasing_Magnitude);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

int util_ccompare(void *a, void *b, util_comp_t comp){
  float a_prime;
  float b_prime;
  switch(comp){
    case util_Increasing:
      a_prime = CREAL_(*((float complex*)a));
      b_prime = CREAL_(*((float complex*)b));
      return util_scompare(&a_prime, &b_prime, util_Increasing);
    case util_Decreasing:
      return -1 * util_ccompare(a, b, util_Increasing);
    case util_Increasing_Magnitude:
      a_prime = cabsf(*((float complex*)a));
      b_prime = cabsf(*((float complex*)b));
      return util_scompare(&a_prime, &b_prime, util_Increasing);
    case util_Decreasing_Magnitude:
      return -1 * util_ccompare(a, b, util_Increasing_Magnitude);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

static void sort(int start, int N, compare_func compare, void *compare_data, swap_func swap, void *swap_data) {
  int i, j;

  if (N <= 1) {
    return;
  }

  i = (rand() % N);

  swap(i, N - 1, swap_data); //pivot rests at position N - 1

  i = start - 1;
  j = N - 1;
  while (i < j) {
    do {
      i++;
    } while (compare(i, N - 1, compare_data) > 0);
    do {
      j--;
    } while (compare(j, N - 1, compare_data) < 0 && i < j);

    if (i < j) {
      swap(i, j, swap_data);
    }
  }
  sort(0, i, compare, compare_data, swap, swap_data);
  swap(i, N - 1, swap_data);
  sort(i + 1, N - i - 1, compare, compare_data, swap, swap_data);
}

static void elem_swap(void *a, void *b, size_t elem_size){
  void *t = malloc(elem_size);
  memcpy(t, b, elem_size);
  memcpy(b, a, elem_size);
  memcpy(a, t, elem_size);
  free(t);
}

typedef struct vec_swap_data{
  void   *V;
  int    incV;
  size_t elem_size;
  int    *P;
  int    incP;
} vec_swap_data_t;

typedef struct mat_row_swap_data{
  char     Order;
  char TransA;
  int               M;
  int               N;
  void              *A;
  int               lda;
  size_t            elem_size;
  int               *P;
  int               incP;
} mat_row_swap_data_t;

static void vec_swap(int a, int b, void *data){
  vec_swap_data_t *d = (vec_swap_data_t*)data;
  elem_swap(d->V + a * d->incV * d->elem_size, d->V + b * d->incV * d->elem_size, d->elem_size);
  if(d->P){
    elem_swap(d->P + a * d->incP, d->P + b * d->incP, sizeof(int));
  }
}

static void mat_row_swap(int a, int b, void *data){
  int i;
  mat_row_swap_data_t *d = (mat_row_swap_data_t*)data;
  switch(d->Order){
    case 'r':
    case 'R':
      switch(d->TransA){
        case 'n':
        case 'N':
          elem_swap(d->A + a * d->lda * d->elem_size, d->A + b * d->lda * d->elem_size, d->N * d->elem_size);
          break;
        default:
          for(i = 0; i < d->M; i++){
            elem_swap(d->A + (a + (i * d->lda)) * d->elem_size, d->A + (b + (i * d->lda)) * d->elem_size, d->elem_size);
          }
          break;
      }
      break;
    default:
      switch(d->TransA){
        case 'n':
        case 'N':
          for(i = 0; i < d->N; i++){
            elem_swap(d->A + (a + (i * d->lda)) * d->elem_size, d->A + (b + (i * d->lda)) * d->elem_size, d->elem_size);
          }
          break;
        default:
          elem_swap(d->A + a * d->lda * d->elem_size, d->A + b * d->lda * d->elem_size, d->M * d->elem_size);
          break;
      }
      break;
  }
  if(d->P){
    elem_swap(d->P + a * d->incP, d->P + b * d->incP, sizeof(int));
  }
}

typedef struct vec_compare_data{
  void              *V;
  int               incV;
  size_t            elem_size;
  elem_compare_func elem_compare;
  int               comp;
} vec_compare_data_t;

typedef struct mat_row_compare_data{
  char              Order;
  char              TransA;
  int               M;
  int               N;
  void              *A;
  int               lda;
  size_t            elem_size;
  elem_compare_func elem_compare;
  int               comp;
  int               col;
} mat_row_compare_data_t;

static int vec_compare(int a, int b, void *data){
  vec_compare_data_t *d = (vec_compare_data_t*)data;
  return d->elem_compare(d->V + a * d->incV * d->elem_size, d->V + b * d->incV * d->elem_size, d->comp);
}

static int mat_row_compare(int a, int b, void *data){
  mat_row_compare_data_t *d = (mat_row_compare_data_t*)data;
  switch(d->Order){
    case 'r':
    case 'R':
      switch(d->TransA){
        case 'n':
        case 'N':
          return d->elem_compare(d->A + (a * d->lda + d->col) * d->elem_size, d->A + (b * d->lda + d->col) * d->elem_size, d->elem_size);
        default:
          return d->elem_compare(d->A + (a + d->lda * d->col) * d->elem_size, d->A + (b + d->lda * d->col) * d->elem_size, d->elem_size);
      }
      break;
    default:
      switch(d->TransA){
        case 'n':
        case 'N':
          return d->elem_compare(d->A + (a + d->lda * d->col) * d->elem_size, d->A + (b + d->lda * d->col) * d->elem_size, d->elem_size);
        default:
          return d->elem_compare(d->A + (a * d->lda + d->col) * d->elem_size, d->A + (b * d->lda + d->col) * d->elem_size, d->elem_size);
      }
      break;
  }
}

static int mat_row_permute_size(char Order, char TransA, int M, int N){
  switch(Order){
    case 'r':
    case 'R':
      switch(TransA){
        case 'n':
        case 'N':
          return M;
        default:
          return N;
      }
    default:
      switch(TransA){
        case 'n':
        case 'N':
          return M;
        default:
          return N;
      }
  }
}

void util_dvec_sort(int N, double *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double),
                               .P         = P,
                               .incP      = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .elem_size    = sizeof(double),
                                     .elem_compare = util_dcompare,
                                     .comp         = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_svec_sort(int N, float *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float),
                               .P         = P,
                               .incP      = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .elem_size    = sizeof(float),
                                     .elem_compare = util_scompare,
                                     .comp         = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_zvec_sort(int N, double complex *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double complex),
                               .P         = P,
                               .incP      = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .elem_size    = sizeof(double complex),
                                     .elem_compare = util_zcompare,
                                     .comp         = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_cvec_sort(int N, float complex *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float complex),
                               .P         = P,
                               .incP      = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .elem_size    = sizeof(float complex),
                                     .elem_compare = util_ccompare,
                                     .comp         = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_dmat_row_sort(char Order, char TransA, int M, int N, double *A, int lda, int *P, int incP, util_comp_t comp, int col){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double),
                                   .P         = P,
                                   .incP      = incP};
  mat_row_compare_data_t compare_data = {.Order        = Order,
                                         .TransA       = TransA,
                                         .M            = M,
                                         .N            = N,
                                         .A            = A,
                                         .lda          = lda,
                                         .elem_size    = sizeof(double),
                                         .elem_compare = util_dcompare,
                                         .comp         = comp,
                                         .col          = col};
  sort(0, mat_row_permute_size(Order, TransA, M, N), &mat_row_compare, &compare_data, &mat_row_swap, &swap_data);
}

void util_smat_row_sort(char Order, char TransA, int M, int N, float *A, int lda, int *P, int incP, util_comp_t comp, int col){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float),
                                   .P         = P,
                                   .incP      = incP};
  mat_row_compare_data_t compare_data = {.Order        = Order,
                                         .TransA       = TransA,
                                         .M            = M,
                                         .N            = N,
                                         .A            = A,
                                         .lda          = lda,
                                         .elem_size    = sizeof(float),
                                         .elem_compare = util_scompare,
                                         .comp         = comp,
                                         .col          = col};
  sort(0, mat_row_permute_size(Order, TransA, M, N), &mat_row_compare, &compare_data, &mat_row_swap, &swap_data);
}

void util_zmat_row_sort(char Order, char TransA, int M, int N, double complex *A, int lda, int *P, int incP, util_comp_t comp, int col){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double complex),
                                   .P         = P,
                                   .incP      = incP};
  mat_row_compare_data_t compare_data = {.Order        = Order,
                                         .TransA       = TransA,
                                         .M            = M,
                                         .N            = N,
                                         .A            = A,
                                         .lda          = lda,
                                         .elem_size    = sizeof(double complex),
                                         .elem_compare = util_zcompare,
                                         .comp         = comp,
                                         .col          = col};
  sort(0, mat_row_permute_size(Order, TransA, M, N), &mat_row_compare, &compare_data, &mat_row_swap, &swap_data);
}

void util_cmat_row_sort(char Order, char TransA, int M, int N, float complex *A, int lda, int *P, int incP, util_comp_t comp, int col){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float complex),
                                   .P         = P,
                                   .incP      = incP};
  mat_row_compare_data_t compare_data = {.Order        = Order,
                                         .TransA       = TransA,
                                         .M            = M,
                                         .N            = N,
                                         .A            = A,
                                         .lda          = lda,
                                         .elem_size    = sizeof(float complex),
                                         .elem_compare = util_ccompare,
                                         .comp         = comp,
                                         .col          = col};
  sort(0, mat_row_permute_size(Order, TransA, M, N), &mat_row_compare, &compare_data, &mat_row_swap, &swap_data);
}

static void reverse(int N, swap_func swap, void *swap_data) {
  int i;
  for (i = 0 ; i < N / 2; i++)
  {
    swap(i, N -i-1, swap_data);
  }
}

void util_dvec_reverse(int N, double* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double),
                               .P         = P,
                               .incP      = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_svec_reverse(int N, float* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float),
                               .P         = P,
                               .incP      = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_zvec_reverse(int N, double complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double complex),
                               .P         = P,
                               .incP      = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_cvec_reverse(int N, float complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float complex),
                               .P         = P,
                               .incP      = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_dmat_row_reverse(char Order, char TransA, int M, int N, double *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double),
                                   .P         = P,
                                   .incP      = incP};
  reverse(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

void util_smat_row_reverse(char Order, char TransA, int M, int N, float *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float),
                                   .P         = P,
                                   .incP      = incP};
  reverse(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

void util_zmat_row_reverse(char Order, char TransA, int M, int N, double complex *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double complex),
                                   .P         = P,
                                   .incP      = incP};
  reverse(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

void util_cmat_row_reverse(char Order, char TransA, int M, int N, float complex *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float complex),
                                   .P         = P,
                                   .incP      = incP};
  reverse(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

static void shuffle(int N, swap_func swap, void *swap_data) {
  int i;
  int gap;
  for (i = 0; i < N - 2; i++)
  {
    gap = rand() % (N - i);
    swap(i, i + gap, swap_data);
  }
}

void util_dvec_shuffle(int N, double* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double),
                               .P         = P,
                               .incP      = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_svec_shuffle(int N, float* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float),
                               .P         = P,
                               .incP      = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_zvec_shuffle(int N, double complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double complex),
                               .P         = P,
                               .incP      = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_cvec_shuffle(int N, float complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float complex),
                               .P         = P,
                               .incP      = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_dmat_row_shuffle(char Order, char TransA, int M, int N, double *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double),
                                   .P         = P,
                                   .incP      = incP};
  shuffle(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

void util_smat_row_shuffle(char Order, char TransA, int M, int N, float *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float),
                                   .P         = P,
                                   .incP      = incP};
  shuffle(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

void util_zmat_row_shuffle(char Order, char TransA, int M, int N, double complex *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double complex),
                                   .P         = P,
                                   .incP      = incP};
  shuffle(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

void util_cmat_row_shuffle(char Order, char TransA, int M, int N, float complex *A, int lda, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float complex),
                                   .P         = P,
                                   .incP      = incP};
  shuffle(mat_row_permute_size(Order, TransA, M, N), &mat_row_swap, &swap_data);
}

int* util_identity_permutation(int N){
  int i;
  int *identity = (int*)malloc(N * sizeof(int));
  for (i = 0; i < N; i++)
  {
    identity[i] = i;
  }
  return identity;
}

int* util_inverse_permutation(int N, int *P, int incP){
  int i;
  int *P_inverse = (int*)malloc(N * sizeof(int));
  for (i = 0; i < N; i++)
  {
    P_inverse[P[i * incP]] = i;
  }
  return P_inverse;
}

static void permute(int N, int *P, int incP, swap_func swap, void *swap_data) {
  int i, j, k;
  int *Q = util_inverse_permutation(N, P, incP);
  for (i = 0; i < N; i++)
  {
    if(Q[i] != -1){
      j = Q[i];
      while(j != i){
        swap(i, j, swap_data);
        k = j;
        j = Q[j];
        Q[k] = -1;
      }
      Q[i] = -1;
    }
  }
  free(Q);
}

void util_dvec_permute(int N, double* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double),
                               .P         = P,
                               .incP      = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_svec_permute(int N, float* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float),
                               .P         = P,
                               .incP      = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_zvec_permute(int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(double complex),
                               .P         = P,
                               .incP      = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_cvec_permute(int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V         = V,
                               .incV      = incV,
                               .elem_size = sizeof(float complex),
                               .P         = P,
                               .incP      = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_dmat_row_permute(char Order, char TransA, int M, int N, double *A, int lda, int *Q, int incQ, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double),
                                   .P         = P,
                                   .incP      = incP};
  permute(mat_row_permute_size(Order, TransA, M, N), Q, incQ, &mat_row_swap, &swap_data);
}

void util_smat_row_permute(char Order, char TransA, int M, int N, float *A, int lda, int *Q, int incQ, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float),
                                   .P         = P,
                                   .incP      = incP};
  permute(mat_row_permute_size(Order, TransA, M, N), Q, incQ, &mat_row_swap, &swap_data);
}

void util_zmat_row_permute(char Order, char TransA, int M, int N, double complex *A, int lda, int *Q, int incQ, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(double complex),
                                   .P         = P,
                                   .incP      = incP};
  permute(mat_row_permute_size(Order, TransA, M, N), Q, incQ, &mat_row_swap, &swap_data);
}

void util_cmat_row_permute(char Order, char TransA, int M, int N, float complex *A, int lda, int *Q, int incQ, int *P, int incP){
  mat_row_swap_data_t swap_data = {.Order     = Order,
                                   .TransA    = TransA,
                                   .M         = M,
                                   .N         = N,
                                   .A         = A,
                                   .lda       = lda,
                                   .elem_size = sizeof(float complex),
                                   .P         = P,
                                   .incP      = incP};
  permute(mat_row_permute_size(Order, TransA, M, N), Q, incQ, &mat_row_swap, &swap_data);
}

double* util_dvec_alloc(int N, int incV) {
  double *V = (double*)malloc(N * incV * sizeof(double));
  if(incV != 1){
    //fill empty space with random data to check increments
    util_dvec_fill(N * incV, V, 1, util_Vec_Rand, 1.0, 1.0);
    util_dvec_fill(N, V, incV, util_Vec_Constant, 0.0, 1.0);
  }
  return V;
}

float* util_svec_alloc(int N, int incV) {
  float *V = (float*)malloc(N * incV * sizeof(float));
  if(incV != 1){
    //fill empty space with random data to check increments
    util_svec_fill(N * incV, V, 1, util_Vec_Rand, 1.0, 1.0);
    util_svec_fill(N, V, incV, util_Vec_Constant, 0.0, 1.0);
  }
  return V;
}

double complex* util_zvec_alloc(int N, int incV) {
  double complex *V = (double complex*)malloc(N * incV * sizeof(double complex));
  if(incV != 1){
    //fill empty space with random data to check increments
    util_zvec_fill(N * incV, V, 1, util_Vec_Rand, 1.0, 1.0);
    util_zvec_fill(N, V, incV, util_Vec_Constant, 0.0, 1.0);
  }
  return V;
}

float complex* util_cvec_alloc(int N, int incV) {
  float complex *V = (float complex*)malloc(N * incV * sizeof(float complex));
  if(incV != 1){
    //fill empty space with random data to check increments
    util_cvec_fill(N * incV, V, 1, util_Vec_Rand, 1.0, 1.0);
    util_cvec_fill(N, V, incV, util_Vec_Constant, 0.0, 1.0);
  }
  return V;
}

double* util_dmat_alloc(char Order, int M, int N, int lda) {
  double *A;
  switch(Order){
    case 'r':
    case 'R':
      A = (double*)malloc(lda * M * sizeof(double));
      //fill empty space with random data to check lda
      util_dmat_fill(Order, 'n', M, lda, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_dmat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
    default:
      A = (double*)malloc(lda * N * sizeof(double));
      //fill empty space with random data to check lda
      util_dmat_fill(Order, 'n', lda, N, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_dmat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
  }
  return A;
}

float* util_smat_alloc(char Order, int M, int N, int lda) {
  float *A;
  switch(Order){
    case 'r':
    case 'R':
      A = (float*)malloc(lda * M * sizeof(float));
      //fill empty space with random data to check lda
      util_smat_fill(Order, 'n', M, lda, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_smat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
    default:
      A = (float*)malloc(lda * N * sizeof(float));
      //fill empty space with random data to check lda
      util_smat_fill(Order, 'n', lda, N, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_smat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
  }
  return A;
}

double complex* util_zmat_alloc(char Order, int M, int N, int lda) {
  double complex *A;
  switch(Order){
    case 'r':
    case 'R':
      A = (double complex*)malloc(lda * M * sizeof(double complex));
      //fill empty space with random data to check lda
      util_zmat_fill(Order, 'n', M, lda, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_zmat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
    default:
      A = (double complex*)malloc(lda * N * sizeof(double complex));
      //fill empty space with random data to check lda
      util_zmat_fill(Order, 'n', lda, N, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_zmat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
  }
  return A;
}

float complex* util_cmat_alloc(char Order, int M, int N, int lda) {
  float complex *A;
  switch(Order){
    case 'r':
    case 'R':
      A = (float complex*)malloc(lda * M * sizeof(float complex));
      //fill empty space with random data to check lda
      util_cmat_fill(Order, 'n', M, lda, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_cmat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
    default:
      A = (float complex*)malloc(lda * N * sizeof(float complex));
      //fill empty space with random data to check lda
      util_cmat_fill(Order, 'n', lda, N, A, lda, util_Mat_Row_Rand, 1.0, 1.0);
      util_cmat_fill(Order, 'n', M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      break;
  }
  return A;
}

void util_dvec_fill(int N, double* V, int incV, util_vec_fill_t fill, double a, double b) {
  int i;
  double small = 1.0 / (1024.0 * 1024.0);			// 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
  switch(fill){
    case util_Vec_Constant_Drop:
    case util_Vec_Constant:
      for (i = 0; i < N; i++) {
        V[i*incV] = 1.0;
      }
      break;
    case util_Vec_Rand_Drop:
    case util_Vec_Rand:
      for (i = 0; i < N; i++) {
        V[i*incV] = util_drand48() * (1+1e-9);
      }
      break;
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_2_Times_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (2 * util_drand48() * (1+1e-9) - 1);
      }
      break;
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = util_drand48() * (1+1e-9) + (util_drand48() * (1+1e-9) - 1);
      }
      break;
    case util_Vec_Normal_Drop:
    case util_Vec_Normal:
      for (i = 0; i < N; i++) {
        double t1 = util_drand48();
        double t2 = util_drand48();
        V[i * incV] = sqrt(-2.0 * log(t1)) * cos(2.0 * PI * t2);
      }
      break;
    case util_Vec_Sine_Drop:
    case util_Vec_Sine:
      for (i = 0; i < N; i++) {
        V[i*incV] = sin(2.0 * PI * ((double)i / (double)N));
      }
      break;
    case util_Vec_Rand_Cond:
      {
        int quart = (N / 8) & ~1;
        int mid = N - quart * 2;
        double c1, c2, c, f;
        int i;

        util_dvec_fill(quart / 2, V, incV, util_Vec_Rand, 1.0e-10, 1.0);
        util_dvec_fill(quart / 2, V + incV*(quart / 2), incV, util_Vec_Rand, 1, 1.0);
        util_dvec_fill(quart / 8, V, incV, util_Vec_Rand, 1e-20, 1.0);
        c1 = 0.0;
        for (i = 0; i < quart; i++) {
          c1 += V[i*incV];
          V[(i + quart)*incV] = -V[i*incV];
        }
        util_dvec_fill(quart, V + 2 * quart*incV, incV, util_Vec_Rand, 1.0, 1.0);
        util_dvec_fill(mid - quart, V + 3 * quart *incV, incV, util_Vec_Rand, 1e-8, 1.0);
        c2 = 0.0; c = 0.0;
        for (i = 2 * quart; i < N; i++) {
          c2 += fabs(V[i*incV]);
          c  += V[i*incV];
        }

        f = 2 * c1 / (b * fabs(c) - c2);
        if (f < 1e-16) f = 1e-16;
        if (f > 1e16) f = 1e16;
        for (i = 2 * quart; i < N; i++) {
          V[i*incV] *= f;
        }

        util_dvec_shuffle(N, V, incV, NULL, 1);
      }
      break;
    case util_Vec_Small_Plus_Increasing_Big:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * i / N;
      }
      break;
    case util_Vec_Small_Plus_Rand_Big:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * util_drand48() * (1+1e-9);
      }
      break;
  }

  if (a != 1.0) {
    for (i = 0; i < N; i++) {
      V[i*incV] *= a;
    }
  }

  switch(fill){
    case util_Vec_Constant_Drop:
    case util_Vec_Rand_Drop:
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Normal_Drop:
    case util_Vec_Sine_Drop:
      for (i = N/2; i < N; i++) {
        V[i*incV] *= 1e-12;
      }
      break;
    default:
      break;
  }
}


void util_svec_fill(int N, float* V, int incV, util_vec_fill_t fill, float a, float b) {
  int i;
  float small = 1.0 / 1024.0;			// 2^-10
  float big   = 1024.0 * 32;		// 2^15
  switch(fill){
    case util_Vec_Constant_Drop:
    case util_Vec_Constant:
      for (i = 0; i < N; i++) {
        V[i*incV] = 1.0;
      }
      break;
    case util_Vec_Rand_Drop:
    case util_Vec_Rand:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)util_drand48() * (1+1e-4);
      }
      break;
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_2_Times_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (2 * (float)util_drand48() * (1+1e-4) - 1);
      }
      break;
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)util_drand48() * (1+1e-4) + ((float)util_drand48() * (1+1e-4) - 1);
      }
      break;
    case util_Vec_Normal_Drop:
    case util_Vec_Normal:
      for (i = 0; i < N; i++) {
        double t1 = util_drand48();
        double t2 = util_drand48();
        V[i * incV] = (float)(sqrt(-2.0 * log(t1)) * cos(2.0 * PI * t2));
      }
      break;
    case util_Vec_Sine_Drop:
    case util_Vec_Sine:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)sin(2.0 * PI * ((float)i / (float)N));
      }
      break;
    case util_Vec_Rand_Cond:
      {
        int quart = (N / 8) & ~1;
        int mid = N - quart * 2;
        double c1, c2, c, f;
        int i;

        util_svec_fill(quart / 2, V, incV, util_Vec_Rand, 1.0e-10, 1.0);
        util_svec_fill(quart / 2, V + (quart / 2) * incV, incV, util_Vec_Rand, 1.0, 1.0);
        util_svec_fill(quart / 8, V, incV, util_Vec_Rand, 1e-20, 1.0);
        c1 = 0.0;
        for (i = 0; i < quart; i++) {
          c1 += V[i*incV];
          V[i + quart] = -V[i*incV];
        }
        util_svec_fill(quart, V + 2 * quart * incV, incV, util_Vec_Rand, 1.0, 1.0);
        util_svec_fill(mid - quart, V + 3 * quart * incV, incV, util_Vec_Rand, 1e-8, 1.0);
        c2 = 0.0; c = 0.0;
        for (i = 2 * quart; i < N; i++) {
          c2 += fabs(V[i*incV]);
          c  += V[i*incV];
        }

        f = 2 * c1 / (b * fabs(c) - c2);
        if (f < 1e-16) f = 1e-16;
        if (f > 1e16) f = 1e16;
        for (i = 2 * quart; i < N; i++) {
          V[i*incV] *= f;
        }

        util_svec_shuffle(N, V, incV, NULL, 1);
      }
      break;
    case util_Vec_Small_Plus_Increasing_Big:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * i / N;
      }
      break;
    case util_Vec_Small_Plus_Rand_Big:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * (float)util_drand48() * (1+1e-4);
      }
      break;
  }

  if (a != 1.0) {
    for (i = 0; i < N; i++) {
      V[i*incV] *= a;
    }
  }

  switch(fill){
    case util_Vec_Constant_Drop:
    case util_Vec_Rand_Drop:
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Normal_Drop:
    case util_Vec_Sine_Drop:
      for (i = N/2; i < N; i++) {
        V[i*incV] *= 1e-12;
      }
      break;
    default:
      break;
  }
}

void util_zvec_fill(int N, double complex* V, int incV, util_vec_fill_t fill, double complex a, double complex b) {
  int i;
  util_dvec_fill(N, (double*)V, incV * 2, fill, 1.0, CREAL_(b));
  switch(fill){
    case util_Vec_Rand_Drop:
    case util_Vec_Rand:
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_2_Times_Rand_Minus_1:
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1:
    case util_Vec_Small_Plus_Rand_Big:
      util_dvec_fill(N, (double*)V + 1, incV * 2, fill, 1.0, 1.0);
      break;
    case util_Vec_Constant_Drop:
    case util_Vec_Constant:
    case util_Vec_Normal_Drop:
    case util_Vec_Normal:
    case util_Vec_Sine_Drop:
    case util_Vec_Sine:
    case util_Vec_Rand_Cond:
    case util_Vec_Small_Plus_Increasing_Big:
      util_dvec_fill(N, (double*)V + 1, incV * 2, util_Vec_Constant, 0.0, 1.0);
      break;
  }
  if (a != 1.0) {
    for (i = 0; i < N; i++) {
      V[i*incV] *= a;
    }
  }
}

void util_cvec_fill(int N, float complex* V, int incV, util_vec_fill_t fill, float complex a, float complex b) {
  int i;
  util_svec_fill(N, (float*)V, incV * 2, fill, 1.0, CREAL_(b));
  switch(fill){
    case util_Vec_Rand_Drop:
    case util_Vec_Rand:
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_2_Times_Rand_Minus_1:
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1:
    case util_Vec_Small_Plus_Rand_Big:
      util_svec_fill(N, (float*)V + 1, incV * 2, fill, 1.0, 1.0);
      break;
    case util_Vec_Constant_Drop:
    case util_Vec_Constant:
    case util_Vec_Normal_Drop:
    case util_Vec_Normal:
    case util_Vec_Sine_Drop:
    case util_Vec_Sine:
    case util_Vec_Rand_Cond:
    case util_Vec_Small_Plus_Increasing_Big:
      util_svec_fill(N, (float*)V + 1, incV * 2, util_Vec_Constant, 0.0, 1.0);
      break;
  }
  if (a != 1.0) {
    for (i = 0; i < N; i++) {
      V[i*incV] *= a;
    }
  }
}

void util_dmat_fill(char Order, char TransA, int M, int N, double* A, int lda, util_mat_fill_t fill, double a, double b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      util_dmat_fill(Order, TransA, M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * lda + i] = 1.0;
      }
      return;
    case util_Mat_Row_Constant_Drop:
      row_fill = util_Vec_Constant_Drop;
      break;
    case util_Mat_Row_Constant:
      row_fill = util_Vec_Constant;
      break;
    case util_Mat_Row_Rand_Drop:
      row_fill = util_Vec_Rand_Drop;
      break;
    case util_Mat_Row_Rand:
      row_fill = util_Vec_Rand;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1_Drop:
      row_fill = util_Vec_2_Times_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1:
      row_fill = util_Vec_2_Times_Rand_Minus_1;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1_Drop:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1;
      break;
    case util_Mat_Row_Normal_Drop:
      row_fill = util_Vec_Normal_Drop;
      break;
    case util_Mat_Row_Normal:
      row_fill = util_Vec_Normal;
      break;
    case util_Mat_Row_Sine_Drop:
      row_fill = util_Vec_Sine_Drop;
      break;
    case util_Mat_Row_Sine:
      row_fill = util_Vec_Sine;
      break;
    case util_Mat_Row_Rand_Cond:
      row_fill = util_Vec_Rand_Cond;
      break;
    case util_Mat_Row_Small_Plus_Increasing_Big:
      row_fill = util_Vec_Small_Plus_Increasing_Big;
      break;
    case util_Mat_Row_Small_Plus_Rand_Big:
      row_fill = util_Vec_Small_Plus_Rand_Big;
      break;
    default:
      exit(125);//TODO better error here
  }
  switch(Order){
    case 'r':
    case 'R':
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_dvec_fill(N, A + i * lda, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_dvec_fill(M, A + i, lda, row_fill, a, b);
          }
          break;
      }
      break;
    default:
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_dvec_fill(N, A + i, lda, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_dvec_fill(M, A + i * lda, 1, row_fill, a, b);
          }
          break;
      }
      break;
  }
}

void util_smat_fill(char Order, char TransA, int M, int N, float* A, int lda, util_mat_fill_t fill, float a, float b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      util_smat_fill(Order, TransA, M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * lda + i] = 1.0;
      }
      return;
    case util_Mat_Row_Constant_Drop:
      row_fill = util_Vec_Constant_Drop;
      break;
    case util_Mat_Row_Constant:
      row_fill = util_Vec_Constant;
      break;
    case util_Mat_Row_Rand_Drop:
      row_fill = util_Vec_Rand_Drop;
      break;
    case util_Mat_Row_Rand:
      row_fill = util_Vec_Rand;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1_Drop:
      row_fill = util_Vec_2_Times_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1:
      row_fill = util_Vec_2_Times_Rand_Minus_1;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1_Drop:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1;
      break;
    case util_Mat_Row_Normal_Drop:
      row_fill = util_Vec_Normal_Drop;
      break;
    case util_Mat_Row_Normal:
      row_fill = util_Vec_Normal;
      break;
    case util_Mat_Row_Sine_Drop:
      row_fill = util_Vec_Sine_Drop;
      break;
    case util_Mat_Row_Sine:
      row_fill = util_Vec_Sine;
      break;
    case util_Mat_Row_Rand_Cond:
      row_fill = util_Vec_Rand_Cond;
      break;
    case util_Mat_Row_Small_Plus_Increasing_Big:
      row_fill = util_Vec_Small_Plus_Increasing_Big;
      break;
    case util_Mat_Row_Small_Plus_Rand_Big:
      row_fill = util_Vec_Small_Plus_Rand_Big;
      break;
    default:
      exit(125); //TODO better error here
  }
  switch(Order){
    case 'r':
    case 'R':
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_svec_fill(N, A + i * lda, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_svec_fill(M, A + i, lda, row_fill, a, b);
          }
          break;
      }
      break;
    default:
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_svec_fill(N, A + i, lda, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_svec_fill(M, A + i * lda, 1, row_fill, a, b);
          }
          break;
      }
    break;
  }
}

void util_zmat_fill(char Order, char TransA, int M, int N, double complex* A, int lda, util_mat_fill_t fill, double complex a, double complex b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      util_zmat_fill(Order, TransA, M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * lda + i] = 1.0;
      }
      return;
    case util_Mat_Row_Constant_Drop:
      row_fill = util_Vec_Constant_Drop;
      break;
    case util_Mat_Row_Constant:
      row_fill = util_Vec_Constant;
      break;
    case util_Mat_Row_Rand_Drop:
      row_fill = util_Vec_Rand_Drop;
      break;
    case util_Mat_Row_Rand:
      row_fill = util_Vec_Rand;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1_Drop:
      row_fill = util_Vec_2_Times_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1:
      row_fill = util_Vec_2_Times_Rand_Minus_1;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1_Drop:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1;
      break;
    case util_Mat_Row_Normal_Drop:
      row_fill = util_Vec_Normal_Drop;
      break;
    case util_Mat_Row_Normal:
      row_fill = util_Vec_Normal;
      break;
    case util_Mat_Row_Sine_Drop:
      row_fill = util_Vec_Sine_Drop;
      break;
    case util_Mat_Row_Sine:
      row_fill = util_Vec_Sine;
      break;
    case util_Mat_Row_Rand_Cond:
      row_fill = util_Vec_Rand_Cond;
      break;
    case util_Mat_Row_Small_Plus_Increasing_Big:
      row_fill = util_Vec_Small_Plus_Increasing_Big;
      break;
    case util_Mat_Row_Small_Plus_Rand_Big:
      row_fill = util_Vec_Small_Plus_Rand_Big;
      break;
    default:
      exit(125);//TODO better error
  }
  switch(Order){
    case 'r':
    case 'R':
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_zvec_fill(N, A + i * lda, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_zvec_fill(M, A + i, lda, row_fill, a, b);
          }
          break;
      }
      break;
    default:
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_zvec_fill(N, A + i, lda, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_zvec_fill(M, A + i * lda, 1, row_fill, a, b);
          }
          break;
      }
    break;
  }
}

void util_cmat_fill(char Order, char TransA, int M, int N, float complex* A, int lda, util_mat_fill_t fill, float complex a, float complex b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      util_cmat_fill(Order, TransA, M, N, A, lda, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * lda + i] = 1.0;
      }
      return;
    case util_Mat_Row_Constant_Drop:
      row_fill = util_Vec_Constant_Drop;
      break;
    case util_Mat_Row_Constant:
      row_fill = util_Vec_Constant;
      break;
    case util_Mat_Row_Rand_Drop:
      row_fill = util_Vec_Rand_Drop;
      break;
    case util_Mat_Row_Rand:
      row_fill = util_Vec_Rand;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1_Drop:
      row_fill = util_Vec_2_Times_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_2_Times_Rand_Minus_1:
      row_fill = util_Vec_2_Times_Rand_Minus_1;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1_Drop:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1_Drop;
      break;
    case util_Mat_Row_Rand_Plus_Rand_Minus_1:
      row_fill = util_Vec_Rand_Plus_Rand_Minus_1;
      break;
    case util_Mat_Row_Normal_Drop:
      row_fill = util_Vec_Normal_Drop;
      break;
    case util_Mat_Row_Normal:
      row_fill = util_Vec_Normal;
      break;
    case util_Mat_Row_Sine_Drop:
      row_fill = util_Vec_Sine_Drop;
      break;
    case util_Mat_Row_Sine:
      row_fill = util_Vec_Sine;
      break;
    case util_Mat_Row_Rand_Cond:
      row_fill = util_Vec_Rand_Cond;
      break;
    case util_Mat_Row_Small_Plus_Increasing_Big:
      row_fill = util_Vec_Small_Plus_Increasing_Big;
      break;
    case util_Mat_Row_Small_Plus_Rand_Big:
      row_fill = util_Vec_Small_Plus_Rand_Big;
      break;
    default:
      exit(125); //TODO better error here
  }
  switch(Order){
    case 'r':
    case 'R':
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_cvec_fill(N, A + i * lda, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_cvec_fill(M, A + i, lda, row_fill, a, b);
          }
          break;
      }
      break;
    default:
      switch(TransA){
        case 'n':
        case 'N':
          for(i = 0; i < M; i++){
            util_cvec_fill(N, A + i, lda, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_cvec_fill(M, A + i * lda, 1, row_fill, a, b);
          }
          break;
      }
    break;
  }
}


const char* vec_fill_name(int type) {
  switch(type){
    case util_Vec_Constant_Drop:
      return "Constant[Drop]";
    case util_Vec_Constant:
      return "Constant";
    case util_Vec_Rand_Drop:
      return "Random[Drop]";
    case util_Vec_Rand:
      return "Random";
    case util_Vec_2_Times_Rand_Minus_1_Drop:
      return "2*Random-1[Drop]";
    case util_Vec_2_Times_Rand_Minus_1:
      return "2*Random-1";
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
      return "Random+(Random-1)[Drop]";
    case util_Vec_Rand_Plus_Rand_Minus_1:
      return "Random+(Random-1)";
    case util_Vec_Normal_Drop:
      return "Normal[Drop]";
    case util_Vec_Normal:
      return "Normal";
    case util_Vec_Sine_Drop:
      return "Sine[Drop]";
    case util_Vec_Sine:
      return "Sine";
    case util_Vec_Rand_Cond:
      return "RandomConditioned";
    case util_Vec_Small_Plus_Increasing_Big:
      return "Small+(i/N*Big)";
    case util_Vec_Small_Plus_Rand_Big:
      return "Small+(Random*Big)";
  }
  return "";
}
