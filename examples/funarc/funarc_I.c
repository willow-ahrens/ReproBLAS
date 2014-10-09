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

double fun( double x )
{
  int k, n = 5;
  double t1, d1 = 1.0;

  t1 = x;

  for( k = 1; k <= n; k++ )
  {
    d1 = 2.0 * d1;
    t1 = t1 + sin (d1 * x) / d1;
  }

  return t1;
}

double arclength(int n)
{
  int i;
  double h = M_PI / n, t1, t2, arc;
  I_double s1; 

  t1 = 0.0;
  dISetZero(s1);

  for( i = 1; i <= n; i++ )
  {
    t2  = fun (i * h);
    arc = sqrt (h*h + (t2 - t1)*(t2 - t1));
    t1  = t2;
    dIAddd(&s1, arc);
  }

  return Iconv2d(s1);

}

int main( int argc, char **argv)
{
  int n = 1000000;
  double ans = 5.795776322412856L;
  double s; 

  double elapsed;
  int iters;

  elapsed = 0.0;
  iters = 0;
  while (elapsed < 0.5 && iters < 1000) {
    tic();
    s = arclength(n);
    elapsed += toc();
    iters++;
  }
  printf("Result :%22.17g.   Error: %5.1e", s, fabs(s - ans) / ans);
  printf(".   Running time: %g", elapsed / iters);
  printf("\n");
  

  return 0;

}


