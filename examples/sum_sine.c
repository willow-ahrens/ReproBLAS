#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <idxd.h>
#include <idxdBLAS.h>
#include <reproBLAS.h>

static struct timeval start;
static struct timeval end;

void tic(void) {
  gettimeofday( &start, NULL );
}

double toc(void) {
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

int main (int argc, char** argv) {
  int n = 1000000;
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
  printf("%15s : Time (s) : |Sum - Sum of Shuffled| = ?\n", "Sum Method");

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

  printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "double", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  // What if we sum x using double double precision?
  tic();
  double ddsum[2] = {0.0, 0.0};
  for(int i = 0; i < n; i++){
    doubledouble_plus_double(ddsum, x[i]);
  }
  sum = ddsum[0] + ddsum[1];
  elapsed_time = toc();

  // Next, we sum the shuffled x
  ddsum[0] = 0.0;
  ddsum[1] = 0.0;
  for(int i = 0; i < n; i++){
    doubledouble_plus_double(ddsum, x_shuffled[i]);
  }
  sum_shuffled = ddsum[0] + ddsum[1];

  printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "doubledouble", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  // Pretty soon we're gonna need some indexed types.
  double_indexed *isum = idxd_dialloc(3);
  double_indexed *itmp = idxd_dialloc(3);

  // Here, we sum x by converting it to an array of indexed types,
  // then we sum the indexed types. This is the most inefficient way to
  // apply the indexed summation algorithm, but it demonstrates the method
  // of adding indexed types.
  tic();
  idxd_disetzero(3, isum);
  for(int i = 0; i < n; i++){
    idxd_didconv(3, x[i], itmp);
    idxd_didiadd(3, itmp, isum);
  }
  sum = idxd_ddiconv(3, isum);
  elapsed_time = toc();

  // Next, we sum the shuffled x
  tic();
  idxd_disetzero(3, isum);
  for(int i = 0; i < n; i++){
    idxd_didconv(3, x_shuffled[i], itmp);
    idxd_didiadd(3, itmp, isum);
  }
  sum_shuffled = idxd_ddiconv(3, isum);
  elapsed_time = toc();

  printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "idxd_didiadd", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  // Here, we sum x using idxd primitives. This is less efficient than the
  // optimized reproBLAS_sum method, but might be useful if the data isn't
  // arranged in a vector.
  tic();
  idxd_disetzero(3, isum);
  for(int i = 0; i < n; i++){
    idxd_didadd(3, x[i], isum);
  }
  sum = idxd_ddiconv(3, isum);
  elapsed_time = toc();

  // Next, we sum the shuffled x
  tic();
  idxd_disetzero(3, isum);
  for(int i = 0; i < n; i++){
    idxd_didadd(3, x_shuffled[i], isum);
  }
  sum_shuffled = idxd_ddiconv(3, isum);
  elapsed_time = toc();

  printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "idxd_didadd", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  // Here, we sum x using idxdBLAS. This shows off the internal methods of
  // reproBLAS.
  tic();
  idxd_disetzero(3, isum);
  idxdBLAS_didsum(3, n, x, 1, isum);
  sum = idxd_ddiconv(3, isum);
  elapsed_time = toc();

  // Next, we sum the shuffled x
  tic();
  idxd_disetzero(3, isum);
  idxdBLAS_didsum(3, n, x_shuffled, 1, isum);
  sum_shuffled = idxd_ddiconv(3, isum);
  elapsed_time = toc();

  printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "idxdBLAS_didsum", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  // Here, we sum x using reproBLAS. This is the fastest and easiest method
  tic();
  sum = reproBLAS_dsum(n, x, 1);
  elapsed_time = toc();

  // Next, we sum the shuffled x
  tic();
  sum_shuffled = reproBLAS_dsum(n, x_shuffled, 1);
  elapsed_time = toc();

  printf("%15s : %-8g : |%.17e - %.17e| = %g\n", "reproBLAS_dsum", elapsed_time, sum, sum_shuffled, fabs(sum - sum_shuffled));

  free(x);
  free(x_shuffled);
}
