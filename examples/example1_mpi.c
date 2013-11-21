#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include <IndexedFP.h>
#include <MPIndexedFP.h>

double v(int i, int n) {
	return  sin(2*M_PI * (i / (double)n - 0.5));
}

int main (int argc, char** args) {
	int n = 1000000;
	int i;
	double s1, s;
	I_double I_s1, I_s;

	int rank, nprocs;

	MPI_Init(&argc, &args);
	RMPI_Init();

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int start, end;
	int q, r;
	q = n / nprocs;
	r = n % nprocs;


	if (rank == 0) {
		printf("Number of processors: %d\n", nprocs);
		printf("Input size          : %d\n", n);
	}

	if (rank < r) {
		start = rank * (q + 1);
		end   = start + q + 1;
	}
	else {
		start = r * (q + 1) + (rank - r) * q;
		end   = start + q;
	}

	// SUMMATION USING PRIMITIVE TYPE
	s1 = 0;
	// local computation
	for (i = start; i < end; i++) {
		s1 += v(i,n);
	}
	// REDUCTION
	MPI_Reduce(&s1, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// COMPARE WITH 1 PROC
	if (rank == 0) {
		s1 = 0.0;
		for (i = 0; i < n; i++) {
			s1 += v(i,n);
		}
		printf("\nNormal summation  \n");
		printf("  >> Sequential  result: %24.17g \n", s1);
		printf("  >>    Parallel result: %24.17g \n", s); 
		printf("  >>               Diff: %24g\n", s1 - s);
	}


	// SUMMATION USING REPRODUCIBLE TYPE
	dISetZero(I_s1);
	for (i = start; i < end; i++) {
		dIAddd(&I_s1, v(i,n));
	}
	// REDUCTION
	MPI_Reduce(&I_s1, &I_s, 1, MPI_IDOUBLE, MPI_RSUM, 0, MPI_COMM_WORLD);
	s = Iconv2d(I_s);

	// COMPARE WITH 1 PROC
	if (rank == 0) {
		dISetZero(I_s1);
		for (i = 0; i < n; i++) {
			dIAddd(&I_s1, v(i,n));
		}
		s1 = Iconv2d(I_s1);
		printf("\nReproducible summation \n");
		printf("  >> Sequential  result: %24.17g \n", s1);
		printf("  >>    Parallel result: %24.17g \n", s); 
		printf("  >>               Diff: %24g\n", s1 - s);
	}

	MPI_Finalize();
}
