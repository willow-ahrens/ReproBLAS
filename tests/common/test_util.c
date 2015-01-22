#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "test_util.h"

#ifdef _MATH__COMPLEX_
#define fabsz(X) cabs(X)
#else
#define fabsz(X) (fabs(X.real)+fabs(X.imag))
#endif

#ifdef _MATH__COMPLEX_
#define fabsc(X) cabsf(X)
#else
#define fabsc(X) (fabs(X.real)+fabs(X.imag))
#endif

typedef int (*compare_func)(int i, int j, void *data);
typedef void (*swap_func)(int i, int j, void *data);
typedef int (*compare_elem_func)(void *a, void *b, util_comp_t comp);

void util_random_seed(void) {
  struct timeval st;
  gettimeofday( &st, NULL );
  srand48((long)(st.tv_usec + 1e6*st.tv_sec));
}

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
      a_prime = fabsz(*((double complex*)a));
      b_prime = fabsz(*((double complex*)b));
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
      a_prime = fabsc(*((float complex*)a));
      b_prime = fabsc(*((float complex*)b));
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

typedef struct vec_swap_data_t{
  void   *V;
  int    incV;
  size_t sizeV;
  int    *P;
  int    incP;
} vec_swap_data_t;

static void vec_swap(int a, int b, void *data){
  int tmpP;
  vec_swap_data_t *d = (vec_swap_data_t*)data;
  void *tmpV = malloc(d->sizeV);
  memcpy(tmpV, d->V + b * d->incV * d->sizeV, d->sizeV);
  memcpy(d->V + b * d->incV * d->sizeV, d->V + a * d->incV * d->sizeV, d->sizeV);
  memcpy(d->V + a * d->incV * d->sizeV, tmpV, d->sizeV);
  free(tmpV);
  if(d->P){
    tmpP = d->P[b * d->incP];
    d->P[b * d->incP] = d->P[a * d->incP];
    d->P[a * d->incP] = tmpP;
  }
}

typedef struct vec_compare_data_t{
  void              *V;
  int               incV;
  size_t            sizeV;
  compare_elem_func compare_elem;
  int               comp;
} vec_compare_data_t;

static int vec_compare(int a, int b, void *data){
  vec_compare_data_t *d = (vec_compare_data_t*)data;
  return d->compare_elem(d->V + a * d->incV * d->sizeV, d->V + a * d->incV * d->sizeV, d->comp);
}

void util_dvec_sort(int N, double *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(double),
                                     .compare_elem = util_dcompare,
                                     .comp        = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_svec_sort(int N, float *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(float),
                                     .compare_elem = util_scompare,
                                     .comp        = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_zvec_sort(int N, double complex *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(double complex),
                                     .compare_elem = util_zcompare,
                                     .comp        = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void util_cvec_sort(int N, float complex *V, int incV, int *P, int incP, util_comp_t comp){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float complex),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(float complex),
                                     .compare_elem = util_ccompare,
                                     .comp        = comp};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

static void reverse(int N, swap_func swap, void *swap_data) {
  int i;
  for (i = 0 ; i < N / 2; i++)
  {
    swap(i, N -i-1, swap_data);
  }
}

void util_dvec_reverse(int N, double* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_svec_reverse(int N, float* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_zvec_reverse(int N, double complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
}

void util_cvec_reverse(int N, float complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float complex),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
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
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_svec_shuffle(int N, float* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_zvec_shuffle(int N, double complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void util_cvec_shuffle(int N, float complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float complex),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

static void permute(int N, int *P, int incP, swap_func swap, void *swap_data) {
  int i, j;
  int *permuted = (int*)malloc(N * sizeof(int));
  for (i = 0; i < N; i++)
  {
    permuted[i] = 0;
  }
  for (i = 0; i < N; i++)
  {
    if(!permuted[i]){
      j = P[i * incP];
      while(j != i){
        swap(i, j, swap_data);
        j = P[j * incP];
        permuted[j] = 1;
      }
      permuted[i] = 1;
    }
  }
  free(permuted);
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

void util_dvec_permute(int N, double* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_svec_permute(int N, float* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_zvec_permute(int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void util_cvec_permute(int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float complex),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
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

double* util_dmat_alloc(rblas_order_t order, int M, int N, int ldA) {
  double *A;
  switch(order){
    case rblas_Row_Major:
      A = (double*)malloc(ldA * M * sizeof(double));
      //fill empty space with random data to check ldA
      util_dmat_fill(order, rblas_No_Trans, M, ldA, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_dmat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
    case rblas_Col_Major:
      A = (double*)malloc(ldA * N * sizeof(double));
      //fill empty space with random data to check ldA
      util_dmat_fill(order, rblas_No_Trans, ldA, N, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_dmat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
  }
  return A;
}

float* util_smat_alloc(rblas_order_t order, int M, int N, int ldA) {
  float *A;
  switch(order){
    case rblas_Row_Major:
      A = (float*)malloc(ldA * M * sizeof(float));
      //fill empty space with random data to check ldA
      util_smat_fill(order, rblas_No_Trans, M, ldA, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_smat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
    case rblas_Col_Major:
      A = (float*)malloc(ldA * N * sizeof(float));
      //fill empty space with random data to check ldA
      util_smat_fill(order, rblas_No_Trans, ldA, N, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_smat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
  }
  return A;
}

double complex* util_zmat_alloc(rblas_order_t order, int M, int N, int ldA) {
  double complex *A;
  switch(order){
    case rblas_Row_Major:
      A = (double complex*)malloc(ldA * M * sizeof(double complex));
      //fill empty space with random data to check ldA
      util_zmat_fill(order, rblas_No_Trans, M, ldA, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_zmat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
    case rblas_Col_Major:
      A = (double complex*)malloc(ldA * N * sizeof(double complex));
      //fill empty space with random data to check ldA
      util_zmat_fill(order, rblas_No_Trans, ldA, N, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_zmat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
  }
  return A;
}

float complex* util_cmat_alloc(rblas_order_t order, int M, int N, int ldA) {
  float complex *A;
  switch(order){
    case rblas_Row_Major:
      A = (float complex*)malloc(ldA * M * sizeof(float complex));
      //fill empty space with random data to check ldA
      util_cmat_fill(order, rblas_No_Trans, M, ldA, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_cmat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
    case rblas_Col_Major:
      A = (float complex*)malloc(ldA * N * sizeof(float complex));
      //fill empty space with random data to check ldA
      util_cmat_fill(order, rblas_No_Trans, ldA, N, A, ldA, util_Mat_Row_Rand, 1.0, 1.0);
      util_cmat_fill(order, rblas_No_Trans, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
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
        V[i*incV] = drand48() * (1+1e-9);
      }
      break;
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_2_Times_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (2 * drand48() * (1+1e-9) - 1);
      }
      break;
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = drand48() * (1+1e-9) + (drand48() * (1+1e-9) - 1);
      }
      break;
    case util_Vec_Normal_Drop:
    case util_Vec_Normal:
      for (i = 0; i < N; i++) {
        double t1 = drand48();
        double t2 = drand48();
        V[i * incV] = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
      }
      break;
    case util_Vec_Sine_Drop:
    case util_Vec_Sine:
      for (i = 0; i < N; i++) {
        V[i*incV] = sin(2.0 * M_PI * ((double)i / (double)N));
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
        V[i*incV] = small + (big - small) * drand48() * (1+1e-9);
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
        V[i*incV] = (float)drand48() * (1+1e-4);
      }
      break;
    case util_Vec_2_Times_Rand_Minus_1_Drop:
    case util_Vec_2_Times_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (2 * (float)drand48() * (1+1e-4) - 1);
      }
      break;
    case util_Vec_Rand_Plus_Rand_Minus_1_Drop:
    case util_Vec_Rand_Plus_Rand_Minus_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)drand48() * (1+1e-4) + ((float)drand48() * (1+1e-4) - 1);
      }
      break;
    case util_Vec_Normal_Drop:
    case util_Vec_Normal:
      for (i = 0; i < N; i++) {
        double t1 = drand48();
        double t2 = drand48();
        V[i * incV] = (float)(sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2));
      }
      break;
    case util_Vec_Sine_Drop:
    case util_Vec_Sine:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)sin(2.0 * M_PI * ((float)i / (float)N));
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
        V[i*incV] = small + (big - small) * (float)drand48() * (1+1e-4);
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

void dmat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, double* A, int ldA, util_mat_fill_t fill, double a, double b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      dmat_fill(order, transA, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * ldA + i] = 1.0;
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
  }
  switch(order){
    case rblas_Row_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_dvec_fill(N, A + i * ldA, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_dvec_fill(M, A + i, ldA, row_fill, a, b);
          }
          break;
      }
      break;
    case rblas_Col_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_dvec_fill(N, A + i, ldA, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_dvec_fill(M, A + i * ldA, 1, row_fill, a, b);
          }
          break;
      }
    break;
  }
}

void smat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, float* A, int ldA, util_mat_fill_t fill, float a, float b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      smat_fill(order, transA, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * ldA + i] = 1.0;
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
  }
  switch(order){
    case rblas_Row_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_svec_fill(N, A + i * ldA, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_svec_fill(M, A + i, ldA, row_fill, a, b);
          }
          break;
      }
      break;
    case rblas_Col_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_svec_fill(N, A + i, ldA, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_svec_fill(M, A + i * ldA, 1, row_fill, a, b);
          }
          break;
      }
    break;
  }
}

void zmat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, double complex* A, int ldA, util_mat_fill_t fill, double complex a, double complex b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      zmat_fill(order, transA, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * ldA + i] = 1.0;
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
  }
  switch(order){
    case rblas_Row_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_zvec_fill(N, A + i * ldA, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_zvec_fill(M, A + i, ldA, row_fill, a, b);
          }
          break;
      }
      break;
    case rblas_Col_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_zvec_fill(N, A + i, ldA, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_zvec_fill(M, A + i * ldA, 1, row_fill, a, b);
          }
          break;
      }
    break;
  }
}

void cmat_fill(rblas_order_t order, rblas_transpose_t transA, int M, int N, float complex* A, int ldA, util_mat_fill_t fill, float complex a, float complex b) {
  int i;
  util_vec_fill_t row_fill;
  switch(fill){
    case util_Mat_Identity:
      cmat_fill(order, transA, M, N, A, ldA, util_Mat_Row_Constant, 0.0, 1.0);
      for(i = 0; i < M && i < N; i++){
        A[i * ldA + i] = 1.0;
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
  }
  switch(order){
    case rblas_Row_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_cvec_fill(N, A + i * ldA, 1, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_cvec_fill(M, A + i, ldA, row_fill, a, b);
          }
          break;
      }
      break;
    case rblas_Col_Major:
      switch(transA){
        case rblas_No_Trans:
          for(i = 0; i < M; i++){
            util_cvec_fill(N, A + i, ldA, row_fill, a, b);
          }
          break;
        default:
          for(i = 0; i < N; i++){
            util_cvec_fill(M, A + i * ldA, 1, row_fill, a, b);
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
