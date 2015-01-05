#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "test_vec.h"


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

#define SWAP_(A,B,T) {(T) = (A); (A) = (B); (B) = (T);}

void svec_sort(int N, float* v, int inc, int type) {
  float pivot;
  float t;
  int i, j;
  if (N <= 1) {
    return;
  }
  i = (rand() % N);

  pivot = v[i * inc];
  v[i * inc] = v[(N - 1) * inc];
  v[(N - 1) * inc] = pivot;

  i = -1;
  j = N - 1;
  switch(type){
    case vec_order_INCREASING:
      while (i < j) {
        do {
          i++;
        } while (v[i * inc] < pivot);
        do {
          j--;
        } while (v[j * inc] > pivot && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING:
      while (i < j) {
        do {
          i++;
        } while (v[i * inc] > pivot);
        do {
          j--;
        } while (v[j * inc] < pivot && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_INCREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabs(v[i * inc]) < fabs(pivot));
        do {
          j--;
        } while (fabs(v[j * inc]) > fabs(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabs(v[i * inc]) > fabs(pivot));
        do {
          j--;
        } while (fabs(v[j * inc]) < fabs(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
  }
  svec_sort(i, v, inc, type);
  v[(N - 1) * inc] = v[i * inc];
  v[i * inc] = pivot;
  svec_sort(N - i - 1, v + (i + 1) * inc, inc, type);
}

void svec_reverse(int N, float* v, int inc) {
  int i;
  float t;
  for (i = 0 ; i < N / 2; i++)
  {
    SWAP_(v[i*inc], v[(N -i-1)*inc], t);
  }
}

void svec_shuffle(int N, float* v, int inc) {
  int i;
  float t;
  int gap;
  for (i = 0; i < N - 2; i++)
  {
    gap = rand() % (N - i);
    SWAP_(v[i*inc], v[(i + gap)*inc], t);
  }
}

void dvec_sort(int N, double* v, int inc, int type) {
  double pivot;
  double t;
  int i, j;
  if (N <= 1) {
    return;
  }
  i = (rand() % N);

  pivot = v[i * inc];
  v[i * inc] = v[(N - 1) * inc];
  v[(N - 1) * inc] = pivot;

  i = -1;
  j = N - 1;
  switch(type){
    case vec_order_INCREASING:
      while (i < j) {
        do {
          i++;
        } while (v[i * inc] < pivot);
        do {
          j--;
        } while (v[j * inc] > pivot && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING:
      while (i < j) {
        do {
          i++;
        } while (v[i * inc] > pivot);
        do {
          j--;
        } while (v[j * inc] < pivot && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_INCREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabs(v[i * inc]) < fabs(pivot));
        do {
          j--;
        } while (fabs(v[j * inc]) > fabs(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabs(v[i * inc]) > fabs(pivot));
        do {
          j--;
        } while (fabs(v[j * inc]) < fabs(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
  }
  dvec_sort(i, v, inc, type);
  v[(N - 1) * inc] = v[i * inc];
  v[i * inc] = pivot;
  dvec_sort(N - i - 1, v + (i + 1) * inc, inc, type);
}

void dvec_reverse(int N, double* v, int inc) {
  int i;
  double t;
  for (i = 0 ; i < N / 2; i++)
  {
    SWAP_(v[i*inc], v[(N -i-1)*inc], t);
  }
}

void dvec_shuffle(int N, double* v, int inc) {
  int i;
  double t;
  int gap;
  for (i = 0 ; i < N - 2; i++)
  {
    gap = rand() % (N - i);
    SWAP_(v[i*inc], v[(i + gap)*inc], t);
  }
}

void cvec_sort(int N, float complex* v, int inc, int type) {
  float complex pivot;
  float complex t;
  int i, j;
  if (N <= 1) {
    return;
  }
  i = (rand() % N);

  pivot = v[i * inc];
  v[i * inc] = v[(N - 1) * inc];
  v[(N - 1) * inc] = pivot;

  i = -1;
  j = N - 1;
  switch(type){
    case vec_order_INCREASING:
      while (i < j) {
        do {
          i++;
        } while (CREAL_(v[i * inc]) < CREAL_(pivot));
        do {
          j--;
        } while (CREAL_(v[j * inc]) > CREAL_(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING:
      while (i < j) {
        do {
          i++;
        } while (CREAL_(v[i * inc]) > CREAL_(pivot));
        do {
          j--;
        } while (CREAL_(v[j * inc]) < CREAL_(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_INCREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabsc(v[i * inc]) < fabsc(pivot));
        do {
          j--;
        } while (fabsc(v[j * inc]) > fabsc(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabsc(v[i * inc]) > fabsc(pivot));
        do {
          j--;
        } while (fabsc(v[j * inc]) < fabsc(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
  }
  cvec_sort(i, v, inc, type);
  v[(N - 1) * inc] = v[i * inc];
  v[i * inc] = pivot;
  cvec_sort(N - i - 1, v + (i + 1) * inc, inc, type);
}

void cvec_reverse(int N, float complex* v, int inc) {
  int i;
  float complex t;
  for (i = 0 ; i < N / 2; i++)
  {
    SWAP_(v[i*inc], v[(N -i-1) * inc], t);
  }
}

void cvec_shuffle(int N, float complex* v, int inc) {
  int i;
  float complex t;
  int gap;
  for (i = 0 ; i < N - 2; i++)
  {
    gap = rand() % (N - i);
    SWAP_(v[i * inc], v[(i + gap) * inc], t);
  }
}

void zvec_sort(int N, double complex* v, int inc, int type) {
  double complex pivot;
  double complex t;
  int i, j;
  if (N <= 1) {
    return;
  }
  i = (rand() % N);

  pivot = v[i * inc];
  v[i * inc] = v[(N - 1) * inc];
  v[(N - 1) * inc] = pivot;

  i = -1;
  j = N - 1;
  switch(type){
    case vec_order_INCREASING:
      while (i < j) {
        do {
          i++;
        } while (ZREAL_(v[i * inc]) < ZREAL_(pivot));
        do {
          j--;
        } while (ZREAL_(v[j * inc]) > ZREAL_(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING:
      while (i < j) {
        do {
          i++;
        } while (ZREAL_(v[i * inc]) > ZREAL_(pivot));
        do {
          j--;
        } while (ZREAL_(v[j * inc]) < ZREAL_(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_INCREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabsz(v[i * inc]) < fabsz(pivot));
        do {
          j--;
        } while (fabsz(v[j * inc]) > fabsz(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
    case vec_order_DECREASING_MAGNITUDE:
      while (i < j) {
        do {
          i++;
        } while (fabsz(v[i * inc]) > fabsz(pivot));
        do {
          j--;
        } while (fabsz(v[j * inc]) < fabsz(pivot) && i < j);

        if (i < j) {
          SWAP_(v[i * inc], v[j * inc], t);
        }
      }
      break;
  }
  zvec_sort(i, v, inc, type);
  v[(N - 1) * inc] = v[i * inc];
  v[i * inc] = pivot;
  zvec_sort(N - i - 1, v + (i + 1) * inc, inc, type);
}

void zvec_reverse(int N, double complex* v, int inc) {
  int i;
  double complex t;
  for (i = 0 ; i < N / 2; i++)
  {
    SWAP_(v[i*inc], v[(N -i-1) * inc], t);
  }
}

void zvec_shuffle(int N, double complex* v, int inc) {
  int i;
  double complex t;
  int gap;
  for (i = 0 ; i < N - 2; i++)
  {
    gap = rand() % (N - i);
    SWAP_(v[i * inc], v[(i + gap) * inc], t);
  }
}

float* svec_alloc(int N, int inc) {
  return (float*)calloc(N * inc, sizeof(float));
}

double* dvec_alloc(int N, int inc) {
  return (double*)calloc(N * inc, sizeof(double));
}

float complex* cvec_alloc(int N, int inc) {
  return (float complex*)calloc(N * inc, sizeof(float complex));
}

double complex* zvec_alloc(int N, int inc) {
  return (double complex*)calloc(N * inc, sizeof(double complex));
}

void vec_random_seed(void) {
  struct timeval st;
  gettimeofday( &st, NULL );
  srand48((long)(st.tv_usec + 1e6*st.tv_sec));
}

void svec_fill(int N, float* v, int inc, int type, float a, float b) {
  int i;
	float small = 1.0 / 1024.0;			// 2^-10
	float big   = 1024.0 * 32;		// 2^15
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
      for (i = 0; i < N; i++) {
        v[i*inc] = 1.0;
      }
      break;
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
      for (i = 0; i < N; i++) {
        v[i*inc] = (float)drand48() * (1+1e-4);
      }
      break;
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        v[i*inc] = (2 * (float)drand48() * (1+1e-4) - 1);
      }
      break;
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        v[i*inc] = (float)drand48() * (1+1e-4) + ((float)drand48() * (1+1e-4) - 1);
      }
      break;
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
      for (i = 0; i < N; i++) {
        double t1 = drand48();
        double t2 = drand48();
        v[i * inc] = (float)(sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2));
      }
      break;
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      for (i = 0; i < N; i++) {
        v[i*inc] = (float)sin(2.0 * M_PI * ((float)i / (float)N));
      }
      break;
    case vec_fill_RAND_COND:
      {
        int quart = (N / 8) & ~1;
        int mid = N - quart * 2;
        double c1, c2, c, f;
        int i;

        svec_fill(quart / 2, v, inc, vec_fill_RAND, 1.0e-10, 1.0);
        svec_fill(quart / 2, v + (quart / 2) * inc, inc, vec_fill_RAND, 1.0, 1.0);
        svec_fill(quart / 8, v, inc, vec_fill_RAND, 1e-20, 1.0);
        c1 = 0.0;
        for (i = 0; i < quart; i++) {
          c1 += v[i*inc];
          v[i + quart] = -v[i*inc];
        }
        svec_fill(quart, v + 2 * quart * inc, inc, vec_fill_RAND, 1.0, 1.0);
        svec_fill(mid - quart, v + 3 * quart * inc, inc, vec_fill_RAND, 1e-8, 1.0);
        c2 = 0.0; c = 0.0;
        for (i = 2 * quart; i < N; i++) {
          c2 += fabs(v[i*inc]);
          c  += v[i*inc];
        }

        f = 2 * c1 / (b * fabs(c) - c2);
        if (f < 1e-16) f = 1e-16;
        if (f > 1e16) f = 1e16;
        for (i = 2 * quart; i < N; i++) {
          v[i*inc] *= f;
        }

        svec_shuffle(N, v, inc);
      }
      break;
    case vec_fill_SMALL_PLUS_INCREASING_BIG:
      for (i = 0; i < N; i++) {
        v[i*inc] = small + (big - small) * i / N;
      }
      break;
    case vec_fill_SMALL_PLUS_RAND_BIG:
      for (i = 0; i < N; i++) {
        v[i*inc] = small + (big - small) * (float)drand48() * (1+1e-4);
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
          v[i*inc] *= a;
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
        v[i*inc] *= 1e-12;
      }
      break;
  }
}

void dvec_fill(int N, double* v, int inc, int type, double a, double b) {
  int i;
	double small = 1.0 / (1024.0 * 1024.0);			// 2^-20
	double big   = 1024.0 * 1024.0 * 1024.0 * 32;	// 2^35
  switch(type){
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
      for (i = 0; i < N; i++) {
        v[i*inc] = 1.0;
      }
      break;
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
      for (i = 0; i < N; i++) {
        v[i*inc] = drand48() * (1+1e-9);
      }
      break;
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        v[i*inc] = (2 * drand48() * (1+1e-9) - 1);
      }
      break;
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      for (i = 0; i < N; i++) {
        v[i*inc] = drand48() * (1+1e-9) + (drand48() * (1+1e-9) - 1);
      }
      break;
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
      for (i = 0; i < N; i++) {
        double t1 = drand48();
        double t2 = drand48();
        v[i * inc] = sqrt(-2.0 * log(t1)) * cos(2.0 * M_PI * t2);
      }
      break;
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
      for (i = 0; i < N; i++) {
        v[i*inc] = sin(2.0 * M_PI * ((double)i / (double)N));
      }
      break;
    case vec_fill_RAND_COND:
      {
        int quart = (N / 8) & ~1;
        int mid = N - quart * 2;
        double c1, c2, c, f;
        int i;

        dvec_fill(quart / 2, v, inc, vec_fill_RAND, 1.0e-10, 1.0);
        dvec_fill(quart / 2, v + inc*(quart / 2), inc, vec_fill_RAND, 1, 1.0);
        dvec_fill(quart / 8, v, inc, vec_fill_RAND, 1e-20, 1.0);
        c1 = 0.0;
        for (i = 0; i < quart; i++) {
          c1 += v[i*inc];
          v[(i + quart)*inc] = -v[i*inc];
        }
        dvec_fill(quart, v + 2 * quart*inc, inc, vec_fill_RAND, 1.0, 1.0);
        dvec_fill(mid - quart, v + 3 * quart *inc, inc, vec_fill_RAND, 1e-8, 1.0);
        c2 = 0.0; c = 0.0;
        for (i = 2 * quart; i < N; i++) {
          c2 += fabs(v[i*inc]);
          c  += v[i*inc];
        }

        f = 2 * c1 / (b * fabs(c) - c2);
        if (f < 1e-16) f = 1e-16;
        if (f > 1e16) f = 1e16;
        for (i = 2 * quart; i < N; i++) {
          v[i*inc] *= f;
        }

        dvec_shuffle(N, v, inc);
      }
      break;
    case vec_fill_SMALL_PLUS_INCREASING_BIG:
      for (i = 0; i < N; i++) {
        v[i*inc] = small + (big - small) * i / N;
      }
      break;
    case vec_fill_SMALL_PLUS_RAND_BIG:
      for (i = 0; i < N; i++) {
        v[i*inc] = small + (big - small) * drand48() * (1+1e-9);
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
          v[i*inc] *= a;
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
        v[i*inc] *= 1e-12;
      }
      break;
  }
}

void cvec_fill(int N, float complex* v, int inc, int type, float complex a, float complex b) {
  int i;
  svec_fill(N, (float*)v, inc * 2, type, 1.0, CREAL_(b));
  switch(type){
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      svec_fill(N, (float*)v + 1, inc * 2, type, 1.0, 1.0);
      break;
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
    case vec_fill_RAND_COND:
      svec_fill(N, (float*)v + 1, inc * 2, vec_fill_CONSTANT, 0.0, 1.0);
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
          v[i*inc] *= a;
        }
      }
      break;
  }
}

void zvec_fill(int N, double complex* v, int inc, int type, double complex a, double complex b) {
  int i;
  dvec_fill(N, (double*)v, inc * 2, type, 1.0, CREAL_(b));
  switch(type){
    case vec_fill_RAND_DROP:
    case vec_fill_RAND:
    case vec_fill_2_TIMES_RAND_MINUS_1_DROP:
    case vec_fill_2_TIMES_RAND_MINUS_1:
    case vec_fill_RAND_PLUS_RAND_MINUS_1_DROP:
    case vec_fill_RAND_PLUS_RAND_MINUS_1:
      dvec_fill(N, (double*)v + 1, inc * 2, type, 1.0, 1.0);
      break;
    case vec_fill_CONSTANT_DROP:
    case vec_fill_CONSTANT:
    case vec_fill_NORMAL_DROP:
    case vec_fill_NORMAL:
    case vec_fill_SINE_DROP:
    case vec_fill_SINE:
    case vec_fill_RAND_COND:
      dvec_fill(N, (double*)v + 1, inc * 2, vec_fill_CONSTANT, 0.0, 1.0);
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
          v[i*inc] *= a;
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
