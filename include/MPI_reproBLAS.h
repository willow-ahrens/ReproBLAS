#ifndef _REPRODUCIBLE_PAR_BLAS__H_
#define _REPRODUCIBLE_PAR_BLAS__H_
#include "binned.h"
#include "binnedBLAS.h"
#include <mpi.h>

//---- double precision
extern void   prdsumI (MPI_Comm comm, int root, int N, double* x, int incx, I_double* s);
extern void   prdasumI(MPI_Comm comm, int root, int N, double* x, int incx, I_double* s);
extern void   prdnrm2I(MPI_Comm comm, int root, int N, double* x, int incx, double* s);
extern void   prddotI (MPI_Comm comm, int root, int N, double* x, int incx,
                double* y, int incy, I_double* s);

extern double prdsum (MPI_Comm comm, int root, int N, double* x, int incx);
extern double prdasum(MPI_Comm comm, int root, int N, double* x, int incx);
extern double prdnrm2(MPI_Comm comm, int root, int N, double* x, int incx);
extern double prddot (MPI_Comm comm, int root, int N, double* x, int incx,
               double* y, int incy);
//---- single precision
extern void prsasumI(MPI_Comm comm, int root, int N, float* x, int incx, I_float* s);
extern void prssumI (MPI_Comm comm, int root, int N, float* x, int incx, I_float* s);
extern void prsnrm2I(MPI_Comm comm, int root, int N, float* x, int incx, void* s);
extern void prsdotI (MPI_Comm comm, int root, int N, float* x, int incx,
              float* y, int incy, I_float* s);

extern float prsasum(MPI_Comm comm, int root, int N, float* x, int incx);
extern float prssum (MPI_Comm comm, int root, int N, float* x, int incx);
extern float prsnrm2(MPI_Comm comm, int root, int N, float* x, int incx);
extern float prsdot (MPI_Comm comm, int root, int N, float* x, int incx, float* y, int incy);
//====
// Complex
extern void prcsumI  (MPI_Comm comm, int root, int N, float complex* x, int incx, I_float_Complex* sum);
extern void prscasumI(MPI_Comm comm, int root, int N, float complex* x, int incx, I_float* sum);
extern void prscnrm2I(MPI_Comm comm, int root, int N, float complex* x, int incx, float* sum);
extern void prcdotuI (MPI_Comm comm, int root, int N, float complex* x, int incx,
                      float complex* y, int incy, I_float_Complex* sum);
extern void prcdotcI (MPI_Comm comm, int root, int N, float complex* x, int incx,
                      float complex* y, int incy, I_float_Complex* sum);

extern float complex prcsum(MPI_Comm comm, int root, int N, float complex* x, int incx);
extern float prscasum(MPI_Comm comm, int root, int N, float complex* x, int incx);
extern float prscnrm2(MPI_Comm comm, int root, int N, float complex* x, int incx);
extern float complex prcdotu(MPI_Comm comm, int root, int N, float complex* x, int incx,
                     float complex* y, int incy);
extern float complex prcdotc(MPI_Comm comm, int root, int N, float complex* x, int incx,
                     float complex* y, int incy);

//---- double complex
extern void przsumI  (MPI_Comm comm, int root, int N, double complex* x, int incx, I_double_Complex* sum);
extern void prdzasumI(MPI_Comm comm, int root, int N, double complex* x, int incx, I_double* sum);
extern void prdznrm2I(MPI_Comm comm, int root, int N, double complex* x, int incx, double* sum);
extern void przdotuI (MPI_Comm comm, int root, int N, double complex* x, int incx,
                      double complex* y, int incy, I_double_Complex* sum);
extern void przdotcI (MPI_Comm comm, int root, int N, double complex* x, int incx,
                      double complex* y, int incy, I_double_Complex* sum);

extern double complex przsum(MPI_Comm comm, int root, int N, double complex* x, int incx);
extern double prdzasum(MPI_Comm comm, int root, int N, double complex* x, int incx);
extern double prdznrm2(MPI_Comm comm, int root, int N, double complex* x, int incx);
extern double complex przdotu(MPI_Comm comm, int root, int N, double complex* x, int incx,
                     double complex* y, int incy);
extern double complex przdotc(MPI_Comm comm, int root, int N, double complex* x, int incx,
                     double complex* y, int incy);

extern void prddotI2(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	double* y, int incy,
	int fold, double* sum,
	double* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
extern void prdsumI2(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	int fold, double* sum,
	double* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
extern void prdasumI2(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	int fold, double* sum,
	double* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
extern void prdnrm2I2(
	MPI_Comm comm, int root,
	int N,
	double* x, int incx,
	int fold, double* sum,
	double* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
//--- single precision ---
extern void prsdotI2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	float* y, int incy,
	int fold, 
	void* sum,
	void* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
extern void prsasumI2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	int fold, 
	void* sum,
	void* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
extern void prssumI2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	int fold, 
	void* sum,
	void* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
extern void prsnrm2I2(
	MPI_Comm comm, int root,
	int N,
	float* x, int incx,
	int fold, 
	void* sum,
	void* local_sum	// WORKING BUFFER TO STORE LOCAL SUM
);
//====
extern double prdasum2(MPI_Comm comm, int root, int N, double* x, int incx,
                int fold);
extern double prdsum2 (MPI_Comm comm, int root, int N, double* x, int incx,
                int fold);
extern double prdnrm22(MPI_Comm comm, int root, int N, double* x, int incx,
                int fold);
extern double prddot2 (MPI_Comm comm, int root, int N,
                double* x, int incx, double* y, int incy,
                int fold);
//--- single precision ---
extern float prsasum2(MPI_Comm comm, int root,
                int N, float* x, int incx, int fold);
extern float prssum2 (MPI_Comm comm, int root,
                int N, float* x, int incx, int fold);
extern float prsnrm22(MPI_Comm comm, int root, int N, float* x, int incx,
                int fold);
extern float prsdot2 (MPI_Comm comm, int root, int N,
                float* x, int incx, float* y, int incy,
                int fold);
//====

void prdgemv(int rank, int nprocs, rblas_order_t Order, rblas_transpose_t TransA, int M, int N, double *myA, int lda, double *myX, int incX, double *myY, int incY);
void prbdgemv(int rank, int nprocs, rblas_order_t Order, rblas_transpose_t TransA, int M, int N, double *myA, int lda, double *myX, int incX, double *myY, int incY);
#endif
