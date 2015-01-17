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
typedef int (*compare_elem_func)(void *a, void *b, int order);

int dutil_compare(void *a, void *b, int order){
  double a_prime;
  double b_prime;
  switch(order){
    case vec_order_INCREASING:
      if(*((double*)a) < *((double*)b)){
        return -1;
      }else if(*((double*)a) > *((double*)b)){
        return 1;
      }else{
        return 0;
      }
    case vec_order_DECREASING:
      return -1 * dutil_compare(a, b, vec_order_INCREASING);
    case vec_order_INCREASING_MAGNITUDE:
      a_prime = fabs(*((double*)a));
      b_prime = fabs(*((double*)b));
      return dutil_compare(&a_prime, &b_prime, vec_order_INCREASING);
    case vec_order_DECREASING_MAGNITUDE:
      return -1 * dutil_compare(a, b, vec_order_INCREASING_MAGNITUDE);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

int sutil_compare(void *a, void *b, int order){
  float a_prime;
  float b_prime;
  switch(order){
    case vec_order_INCREASING:
      if(*((float*)a) < *((float*)b)){
        return -1;
      }else if(*((float*)a) > *((float*)b)){
        return 1;
      }else{
        return 0;
      }
    case vec_order_DECREASING:
      return -1 * sutil_compare(a, b, vec_order_INCREASING);
    case vec_order_INCREASING_MAGNITUDE:
      a_prime = fabs(*((float*)a));
      b_prime = fabs(*((float*)b));
      return sutil_compare(&a_prime, &b_prime, vec_order_INCREASING);
    case vec_order_DECREASING_MAGNITUDE:
      return -1 * sutil_compare(a, b, vec_order_INCREASING_MAGNITUDE);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

int zutil_compare(void *a, void *b, int order){
  double a_prime;
  double b_prime;
  switch(order){
    case vec_order_INCREASING:
      a_prime = ZREAL_(*((double complex*)a));
      b_prime = ZREAL_(*((double complex*)b));
      return dutil_compare(&a_prime, &b_prime, vec_order_INCREASING);
    case vec_order_DECREASING:
      return -1 * zutil_compare(a, b, vec_order_INCREASING);
    case vec_order_INCREASING_MAGNITUDE:
      a_prime = fabsz(*((double complex*)a));
      b_prime = fabsz(*((double complex*)b));
      return dutil_compare(&a_prime, &b_prime, vec_order_INCREASING);
    case vec_order_DECREASING_MAGNITUDE:
      return -1 * zutil_compare(a, b, vec_order_INCREASING_MAGNITUDE);
    default:
      fprintf(stderr, "error: default case (%s:%d)\n", __FILE__, __LINE__);
      exit(125);
  }
}

int cutil_compare(void *a, void *b, int order){
  float a_prime;
  float b_prime;
  switch(order){
    case vec_order_INCREASING:
      a_prime = CREAL_(*((float complex*)a));
      b_prime = CREAL_(*((float complex*)b));
      return sutil_compare(&a_prime, &b_prime, vec_order_INCREASING);
    case vec_order_DECREASING:
      return -1 * cutil_compare(a, b, vec_order_INCREASING);
    case vec_order_INCREASING_MAGNITUDE:
      a_prime = fabsc(*((float complex*)a));
      b_prime = fabsc(*((float complex*)b));
      return sutil_compare(&a_prime, &b_prime, vec_order_INCREASING);
    case vec_order_DECREASING_MAGNITUDE:
      return -1 * cutil_compare(a, b, vec_order_INCREASING_MAGNITUDE);
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
  int               order;
} vec_compare_data_t;

static int vec_compare(int a, int b, void *data){
  vec_compare_data_t *d = (vec_compare_data_t*)data;
  return d->compare_elem(d->V + a * d->incV * d->sizeV, d->V + a * d->incV * d->sizeV, d->order);
}

void dvec_sort(int N, double *V, int incV, int *P, int incP, int order){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(double),
                                     .compare_elem = dutil_compare,
                                     .order        = order};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void svec_sort(int N, float *V, int incV, int *P, int incP, int order){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(float),
                                     .compare_elem = sutil_compare,
                                     .order        = order};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void zvec_sort(int N, double complex *V, int incV, int *P, int incP, int order){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(double complex),
                                     .compare_elem = zutil_compare,
                                     .order        = order};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

void cvec_sort(int N, float complex *V, int incV, int *P, int incP, int order){
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float complex),
                               .P     = P,
                               .incP  = incP};
  vec_compare_data_t compare_data = {.V            = V,
                                     .incV         = incV,
                                     .sizeV        = sizeof(float complex),
                                     .compare_elem = cutil_compare,
                                     .order        = order};
  sort(0, N, &vec_compare, &compare_data, &vec_swap, &swap_data);
}

static void reverse(int N, swap_func swap, void *swap_data) {
  int i;
  for (i = 0 ; i < N / 2; i++)
  {
    swap(i, N -i-1, swap_data);
  }
}

void dvec_reverse(int N, double* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
}

void svec_reverse(int N, float* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
}

void zvec_reverse(int N, double complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  reverse(N, &vec_swap, &swap_data);
}

void cvec_reverse(int N, float complex* V, int incV, int *P, int incP) {
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

void dvec_shuffle(int N, double* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void svec_shuffle(int N, float* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void zvec_shuffle(int N, double complex* V, int incV, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  shuffle(N, &vec_swap, &swap_data);
}

void cvec_shuffle(int N, float complex* V, int incV, int *P, int incP) {
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

void dvec_permute(int N, double* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void svec_permute(int N, float* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void zvec_permute(int N, double complex* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(double complex),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

void cvec_permute(int N, float complex* V, int incV, int *Q, int incQ, int *P, int incP) {
  vec_swap_data_t swap_data = {.V     = V,
                               .incV  = incV,
                               .sizeV = sizeof(float complex),
                               .P     = P,
                               .incP  = incP};
  permute(N, Q, incQ, &vec_swap, &swap_data);
}

double* dvec_alloc(int N, int incV) {
  return (double*)calloc(N * incV, sizeof(double));
}

float* svec_alloc(int N, int incV) {
  return (float*)calloc(N * incV, sizeof(float));
}

double complex* zvec_alloc(int N, int incV) {
  return (double complex*)calloc(N * incV, sizeof(double complex));
}

float complex* cvec_alloc(int N, int incV) {
  return (float complex*)calloc(N * incV, sizeof(float complex));
}

double* dmat_alloc(int M, int N, int ldA) {
  return (double*)calloc(ldA * N, sizeof(double));
}

float* smat_alloc(int M, int N, int ldA) {
  return (float*)calloc(ldA * N, sizeof(float));
}

double complex* zmat_alloc(int M, int N, int ldA) {
  return (double complex*)calloc(ldA * N, sizeof(double complex));
}

float complex* cmat_alloc(int M, int N, int ldA) {
  return (float complex*)calloc(ldA * N, sizeof(float complex));
}

void util_random_seed(void) {
  struct timeval st;
  gettimeofday( &st, NULL );
  srand48((long)(st.tv_usec + 1e6*st.tv_sec));
}

void dvec_fill(int N, double* V, int incV, int type, double a, double b) {
  int i;
  double small = 1.0 / (1024.0 * 1024.0);			// 2^-20
  double big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
      for (i = 0; i < N; i++) {
        V[i*incV] = 1.0;
      }
      break;
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
      for (i = 0; i < N; i++) {
        V[i*incV] = drand48() * (1+1e-9);
      }
      break;
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (2 * drand48() * (1+1e-9) - 1);
      }
      break;
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = drand48() * (1+1e-9) + (drand48() * (1+1e-9) - 1);
      }
      break;
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
      for (i = 0; i < N; i++) {
        double t1 = drand48();
        double t2 = drand48();
        V[i * incV] = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
      }
      break;
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      for (i = 0; i < N; i++) {
        V[i*incV] = sin(2.0 * M_PI * ((double)i / (double)N));
      }
      break;
    case vec_fill_RAND_COND:
      {
        int quart = (N / 8) & ~1;
        int mid = N - quart * 2;
        double c1, c2, c, f;
        int i;

        dvec_fill(quart / 2, V, incV, vec_fill_RAND, 1.0e-10, 1.0);
        dvec_fill(quart / 2, V + incV*(quart / 2), incV, vec_fill_RAND, 1, 1.0);
        dvec_fill(quart / 8, V, incV, vec_fill_RAND, 1e-20, 1.0);
        c1 = 0.0;
        for (i = 0; i < quart; i++) {
          c1 += V[i*incV];
          V[(i + quart)*incV] = -V[i*incV];
        }
        dvec_fill(quart, V + 2 * quart*incV, incV, vec_fill_RAND, 1.0, 1.0);
        dvec_fill(mid - quart, V + 3 * quart *incV, incV, vec_fill_RAND, 1e-8, 1.0);
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

        dvec_shuffle(N, V, incV, NULL, 1);
      }
      break;
    case vec_fill_SMALL_PLUS_INCREASING_BIG:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * i / N;
      }
      break;
    case vec_fill_SMALL_PLUS_RAND_BIG:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * drand48() * (1+1e-9);
      }
      break;
  }

  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      if (a != 1.0) {
        for (i = 0; i < N; i++) {
          V[i*incV] *= a;
        }
      }
      break;
  }
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_RAND_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_NORMAL_DROP:
    case vec_fill_SINE_DROP:
      for (i = N/2; i < N; i++) {
        V[i*incV] *= 1e-12;
      }
      break;
  }
}


void svec_fill(int N, float* V, int incV, int type, float a, float b) {
  int i;
  float small = 1.0 / 1024.0;			// 2^-10
  float big   = 1024.0 * 32;		// 2^15
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
      for (i = 0; i < N; i++) {
        V[i*incV] = 1.0;
      }
      break;
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)drand48() * (1+1e-4);
      }
      break;
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (2 * (float)drand48() * (1+1e-4) - 1);
      }
      break;
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)drand48() * (1+1e-4) + ((float)drand48() * (1+1e-4) - 1);
      }
      break;
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
      for (i = 0; i < N; i++) {
        double t1 = drand48();
        double t2 = drand48();
        V[i * incV] = (float)(sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2));
      }
      break;
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      for (i = 0; i < N; i++) {
        V[i*incV] = (float)sin(2.0 * M_PI * ((float)i / (float)N));
      }
      break;
    case vec_fill_RAND_COND:
      {
        int quart = (N / 8) & ~1;
        int mid = N - quart * 2;
        double c1, c2, c, f;
        int i;

        svec_fill(quart / 2, V, incV, vec_fill_RAND, 1.0e-10, 1.0);
        svec_fill(quart / 2, V + (quart / 2) * incV, incV, vec_fill_RAND, 1.0, 1.0);
        svec_fill(quart / 8, V, incV, vec_fill_RAND, 1e-20, 1.0);
        c1 = 0.0;
        for (i = 0; i < quart; i++) {
          c1 += V[i*incV];
          V[i + quart] = -V[i*incV];
        }
        svec_fill(quart, V + 2 * quart * incV, incV, vec_fill_RAND, 1.0, 1.0);
        svec_fill(mid - quart, V + 3 * quart * incV, incV, vec_fill_RAND, 1e-8, 1.0);
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

        svec_shuffle(N, V, incV, NULL, 1);
      }
      break;
    case vec_fill_SMALL_PLUS_INCREASING_BIG:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * i / N;
      }
      break;
    case vec_fill_SMALL_PLUS_RAND_BIG:
      for (i = 0; i < N; i++) {
        V[i*incV] = small + (big - small) * (float)drand48() * (1+1e-4);
      }
      break;
  }

  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      if (a != 1.0) {
        for (i = 0; i < N; i++) {
          V[i*incV] *= a;
        }
      }
      break;
  }
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_RAND_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_NORMAL_DROP:
    case vec_fill_SINE_DROP:
      for (i = N/2; i < N; i++) {
        V[i*incV] *= 1e-12;
      }
      break;
  }
}

void zvec_fill(int N, double complex* V, int incV, int type, double complex a, double complex b) {
  int i;
  dvec_fill(N, (double*)V, incV * 2, type, 1.0, CREAL_(b));
  switch(type){
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      dvec_fill(N, (double*)V + 1, incV * 2, type, 1.0, 1.0);
      break;
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
    case vec_fill_RAND_COND:
      dvec_fill(N, (double*)V + 1, incV * 2, vec_fill_CONSTANT, 0.0, 1.0);
      break;
  }
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      if (a != 1.0) {
        for (i = 0; i < N; i++) {
          V[i*incV] *= a;
        }
      }
      break;
  }
}

void cvec_fill(int N, float complex* V, int incV, int type, float complex a, float complex b) {
  int i;
  svec_fill(N, (float*)V, incV * 2, type, 1.0, CREAL_(b));
  switch(type){
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      svec_fill(N, (float*)V + 1, incV * 2, type, 1.0, 1.0);
      break;
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
    case vec_fill_RAND_COND:
      svec_fill(N, (float*)V + 1, incV * 2, vec_fill_CONSTANT, 0.0, 1.0);
      break;
  }
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      if (a != 1.0) {
        for (i = 0; i < N; i++) {
          V[i*incV] *= a;
        }
      }
      break;
  }
}

const char* vec_fill_name(int type) {
  switch(type){
    case vec_fill_CONSTANT_DROP:
      return "Constant[Drop]";
    case vec_fill_CONSTANT:
      return "Constant";
    case vec_fill_RAND_DROP:
      return "Random[Drop]";
    case vec_fill_RAND:
      return "Random";
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
      return "2*Random-1[Drop]";
    case vec_fill_2_TIMES_RAND_MINUS_1:
      return "2*Random-1";
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
      return "Random+(Random-1)[Drop]";
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      return "Random+(Random-1)";
    case vec_fill_NORMAL_DROP:
      return "Normal[Drop]";
    case vec_fill_NORMAL:
      return "Normal";
    case vec_fill_SINE_DROP:
      return "Sine[Drop]";
    case vec_fill_SINE:
      return "Sine";
    case vec_fill_RAND_COND:
      return "RandomConditioned";
    case vec_fill_SMALL_PLUS_INCREASING_BIG:
      return "Small+(i/N*Big)";
    case vec_fill_SMALL_PLUS_RAND_BIG:
      return "Small+(Random*Big)";
  }
  return "";
}
