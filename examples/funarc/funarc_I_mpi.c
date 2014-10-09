#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <mpi.h>

#include <MPIndexedFP.h>

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
  int i, j, k;
  double h, t1, t2, dppi;
  double tmp;
  Idouble I_s1, I_s; 

  int start, end;
  int rank, nprocs;

  // COMMUNICATOR INFORMATION
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // PARTITIONING
  int q, r;
  q = n / nprocs;
  r = n % nprocs;

  if (rank < r) {
    start = (q + 1) * rank;
    end   = start + q + 1;
  }
  else {
    start = (q + 1) * r + q * (rank - r);
    end   = start + q;
  }

  // LOCAL COMPUTATION
  t1 = -1.0;
  dppi = acos(t1);
  h = dppi / n;
  t1 = fun(start * h);

  dISetZero(I_s1);
  for( i = start + 1; i <= end; i++ )
  {
    t2 = fun (i * h);
    tmp = sqrt (h*h + (t2 - t1)*(t2 - t1));
    t1 = t2;
    dIAddd(&I_s1, tmp);
  }

  // REDUCE FINAL RESULT
  MPI_Reduce(&I_s1, &I_s, 1, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);

  return Iconv2d(I_s);

}

int main( int argc, char **argv)
{
  int n = 1000000;
  double ans = 5.795776322412856L;
  double s; 

  double elapsed;
  int iters;
  int rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  RMPI_Init();

  elapsed = 0.0;
  iters = 0;
  while (elapsed < 0.5 && iters < 1000) {
    tic();
    s = arclength(n);
    elapsed += toc();
    iters++;
  }

  if (rank == 0) {
    printf("Result :%22.17g.   Error: %5.1e", s, fabs(s - ans) / ans);
    printf(".   Running time: %g", elapsed / iters);
    printf("\n");
  }

  MPI_Finalize();
  return 0;

}


