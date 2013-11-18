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

static struct timeval start;
static struct timeval end;

#define NB 1024

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
  int i = 0, j, lN, count = 0;
  Idouble s1;
  Idouble s;
  double s2;

  dISetZero(s1);
  for( ; i < n; i+= NB, v += NB ){
    lN = NB < (n-i) ? NB : (n-i);

    s2 = 0;
    for (j = 0; j < lN; j++)
      s2 += v[j];

//    dIAddd(&s1, s2);
    dIUpdate_(s1, fabs(s2));
    dIAddd_(s1, s2);
    if (++count > 1024) {
      dIRenorm_(s1);
      count = 0;
    }
  }

  MPI_Reduce(&s1, &s, 1, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);

  return Iconv2d(s);

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
  int nnb = (N + NB - 1) / NB;
  q = nnb / nprocs;
  r = nnb % nprocs;


  int* sendcnts = (int*) malloc(nprocs * sizeof(int));
  int* displs   = (int*) malloc(nprocs * sizeof(int));
  for (i = 0; i < r; i++) {
    sendcnts[i] = (q + 1) * NB;
    displs[i]   = ((q + 1) * i) * NB;
    if (displs[i] > N)
      displs[i] = N;
    if (sendcnts[i] + displs[i] > N)
      sendcnts[i] = N - displs[i];
  }
  for (; i < nprocs; i++) {
    sendcnts[i] = q * NB;
    displs[i]   = ((q + 1) * r + q * (i - r)) * NB;
    if (displs[i] > N)
      displs[i] = N;
    if (sendcnts[i] + displs[i] > N)
      sendcnts[i] = N - displs[i];
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

  MPI_Barrier(MPI_COMM_WORLD);
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


