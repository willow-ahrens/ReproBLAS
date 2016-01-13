#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <idxd.h>
#include <idxdBLAS.h>
#include <reproBLAS.h>

static struct timeval start;
static struct timeval end;

void tic() {
  gettimeofday( &start, NULL );
}

double toc( )
{
  gettimeofday( &end, NULL );

  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void doubledouble_plus_double(double* a, double b) {
  double bv;
  double s1, s2, t1, t2;

  // Add two hi words
  s1 = a[0] + b;
  bv = s1 - a[0];
  s2 = ((b - bv) + (a[0] - (s1 - bv)));

  t1 = a[1] + s2;
  bv = t1 - a[1];
  t2 = ((s2 - bv) + (a[1] - (t1 - bv)));

  s2 = t1;

  // Renormalize (s1, s2)  to  (t1, s2)
  t1 = s1 + s2;
  t2 += s2 - (t1 - s1);

  // Renormalize (t1, t2)
  a[0] = t1 + t2;
  a[1] = t2 - (a[0] - t1);
}

int main (int argc, char** args) {
  int n = 100000;
  double *x = malloc(n * sizeof(double));
  double *x_shuffled = malloc(n * sizeof(double));
  double sum;
  double sum_shuffled;
  double elapsed_time;

  printf("Sum of sin(2* M_PI * (i / (double)n - 0.5)).  n = %d.\n\n", n);

  // Set x to be a sine wave
  for(int i = 0; i < n; i++){
    x[i] = sin(2 * M_PI * (i / (double)n - 0.5));
  }

  // Shuffle x into x_shuffled
  for(int i = 0; i < n; i++){
    x_shuffled[i] = x[i];
  }
  double t;
  int r;
  for(int i = 0; i < n; i++){
    r = rand();
    t = x_shuffled[i];
    x_shuffled[i] = x_shuffled[i + (r % (n - i))];
    x_shuffled[i + (r % (n - i))] = t;
  }

  // Make a header
  printf("%10s : Time (s)\t: |Sum - Sum of Shuffled| = ?\n", "Sum Method");

  // First, we sum x using double precision
  tic();
  sum = 0;
  for(int i = 0; i < n; i++){
    sum += x[i];
  }
  elapsed_time = toc();

  // Next, we sum the shuffled x
  sum_shuffled = 0;
  for(int i = 0; i < n; i++){
    sum_shuffled += x_shuffled[i];
  }

  printf("%10s : %.3g\t: |%.17e - %.17e| = %g\n", "double", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  free(x);
  free(x_shuffled);
}
