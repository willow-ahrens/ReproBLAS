#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

#include <IndexedFP.h>

#define NB 1024

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
  int i = 0, j;
  Idouble s;
  double s1;
  int count = 0;
  int lN;

  dISetZero(s);
  for(; i < n; i += NB, v += NB ){
    lN = NB < (n-i) ? NB : (n-i);

    s1 = 0;
    for (j = 0; j < lN; j++)
      s1 += v[j];

#if 1
    dIAddd(&s, s1);
#else
    dIUpdate_(s, fabs(s1));
    dIAddd_(s, s1);
    if (++count > 1024) {
      dIRenorm_(s);
      count = 0;
    }
#endif
  }
  dIRenorm_(s);

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


