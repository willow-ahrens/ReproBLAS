#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

#include <IndexedFP.h>

#define TIMING

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

double sum(int n, double* v)
{
  int i = 0;
  Idouble s;
  int count = 0;

  dISetZero(s);
  for(; i < n; i++ ){
    dIUpdate_(s, fabs(v[i]));
    dIAddd_(s, v[i]);
    if (++count > 1024) {
      dIRenorm_(s);
      count = 0;
    }
  }

  return Iconv2d(s);

}

int main( int argc, char **argv)
{
  int n = 1000000;
  double elapsed = 0.0;
  int iters = 0;
  double s;
  int i;

  double* v = (double*) malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    v[i] = sin(M_PI - 2 * M_PI * i / (double)(n-1));
  }

  while (elapsed < 0.5 && iters < 1000) {
    tic();
    s = sum(n, v);
    elapsed += toc();
    iters++;
  }
  printf("Result :   %22.17g", s);
  printf(".   Running time: %e", elapsed / iters);
  printf("\n");
  

  return 0;

}


