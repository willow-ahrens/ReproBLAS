/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "Indexed.h"
#include "dIndexedMPI.h"

// DEFAULT FOLD (3)
void dIAdd_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
) {
	int i;
	double* pin;
	double* pout;

	int size;
	int fold;

	// SIZE OF TYPE: 2 * FOLD * SIZEOF(DOUBLE)
	MPI_Type_size(*datatype, &size);

	fold = size / (2 * sizeof(double));

	for (i = 0; i < *len; i++) {
		pin  = (double*) invec;
		pout = (double*) inoutvec;

		dIAdd1(fold, pout, pout + fold, 1, pin, pin + fold, 1);
		dIRenorm1(fold, pout, pout + fold, 1);

		invec    += size;
		inoutvec += size;
	}
}

void zIAdd_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
) {
	int i;
	double complex* pin;
	double complex* pout;

	int size;
	int fold;

	// SIZE OF TYPE: 2 * FOLD * SIZEOF(DOUBLE)
	MPI_Type_size(*datatype, &size);

	fold = size / (2 * sizeof(double complex));

	for (i = 0; i < *len; i++) {
		pin  = (double complex*) invec;
		pout = (double complex*) inoutvec;

		zIAdd1(fold, pout, pout + fold, 1, pin, pin + fold, 1);
		zIRenorm1(fold, pout, pout + fold, 1);

		invec    += size;
		inoutvec += size;
	}
}

void dINrm2_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
) {
	int i, j;
	double* pin;
	double* pout;
	int fold;
	double scale1, scale2;

	int size;

	// SIZE OF TYPE: 2 * FOLD * SIZEOF(DOUBLE)
	MPI_Type_size(*datatype, &size);

	fold = (size - sizeof(double)) / (2 * sizeof(double));

	for (i = 0; i < *len; i++) {
		pin  = (double*) invec;
		pout = (double*) inoutvec;

		scale1 = pout[0];
		scale2 = pin[0];

		if (scale1 < scale2) {
			scale1 = scale1 / scale2;
			scale1 = scale1 * scale1;
			for (j = 1; j < 1 + fold; j++)
				pout[j] *= scale1;
			pout[0] = scale2;
		}
		if (scale1 > scale2) {
			scale2 = scale2 / scale1;
			scale2 = scale2 * scale2;
			for (j = 1; j < 1 + fold; j++)
				pin[j] *= scale2;
			pin[0] = scale1;
		}

		dIAdd1(fold, pout + 1, pout + 1 + fold, 1, pin + 1, pin + 1 + fold, 1);
		dIRenorm1(fold, pout + 1, pout + fold + 1, 1);

		invec    += size;
		inoutvec += size;
	}
}

// CREATE TYPE
int dIMPICreate_(MPI_Datatype* newtype, int size) {
   	int ret = MPI_Type_contiguous(size, MPI_DOUBLE, newtype); 
   	MPI_Type_commit( newtype);

	if (ret == MPI_SUCCESS)
		return MPI_SUCCESS;

	MPI_Abort(MPI_COMM_WORLD, ret);

	return ret;
}

int dIMPICreate(int fold, MPI_Datatype* newtype) {
	return dIMPICreate_(newtype, 2 * fold);
}

int dIMPIScaleCreate(int fold, MPI_Datatype* newtype) {
	return dIMPICreate_(newtype, 1 + 2 * fold);
}

int zIMPICreate(int fold, MPI_Datatype* newtype) {
   	int ret = MPI_Type_contiguous(2*fold, MPI_DOUBLE_COMPLEX, newtype); 
   	MPI_Type_commit( newtype);

	if (ret == MPI_SUCCESS)
		return MPI_SUCCESS;

	MPI_Abort(MPI_COMM_WORLD, ret);

	return ret;
}

