#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <mpi.h>
#include <rblas_mpi.h>

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
  return prdsum(MPI_COMM_WORLD, 0, n, v, 1);

}

int main( int argc, char **argv)
{
  int N = 1000000;
  int n;
  double elapsed = 0.0;
  int iters = 0;
  double s;
  int i;
  double* V;
  double* v;

  int start, end;
  int rank, nprocs;

  MPI_Init(&argc, &argv);
  RMPI_Init();

  // COMMUNICATOR INFORMATION
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (rank == 0) {
    V =  (double*) malloc(N * sizeof(double));
    for (i = 0; i < N; i++) {
      V[i] = sin(M_PI - 2 * M_PI * i / (double)(N-1));
    }
  }
  // PARTITIONING
  int q, r;
  q = N / nprocs;
  r = N % nprocs;


  int* sendcnts = (int*) malloc(nprocs * sizeof(int));
  int* displs   = (int*) malloc(nprocs * sizeof(int));
  for (i = 0; i < r; i++) {
    sendcnts[i] = q + 1;
    displs[i]   = (q + 1) * i;
  }
  for (; i < nprocs; i++) {
    sendcnts[i] = q;
    displs[i]   = (q + 1) * r + q * (i - r);
  }

  n = sendcnts[rank];

  // local buffer
  v = (double*) malloc(n * sizeof(double));
  // DISTRIBUTE
  MPI_Scatterv(V, sendcnts, displs, MPI_DOUBLE, v, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  while (iters < 100) {
    tic();
    s = sum(n, v);
    elapsed += toc();
    iters++;
  }

  if (rank == 0) {
    printf("Result :   %22.17g", s);
    printf(".   Running time: %e", elapsed / iters);
    printf("\n");
  }

  free(sendcnts);
  free(displs);
  free(v);
  if ( rank == 0 )
    free(V);

  MPI_Finalize();

  return 0;

}


