/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "Indexed.h"
#include "sIndexedMPI.h"
#include <stdio.h>

void sIAdd_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
) {
	int i;
	int* pfold;
	int fold;
	int size;

	int NI, NA, ND, combiner;
	int array_of_integers[2];
	MPI_Aint array_of_addresses[2];
	MPI_Datatype array_of_datatypes[2];

	MPI_Type_get_envelope(*datatype, &NI, &NA, &ND, &combiner);
	MPI_Type_get_contents(*datatype, NI, NA, ND,
		array_of_integers, array_of_addresses, array_of_datatypes);

	MPI_Type_size(*datatype, &size);

	fold = size / (sizeof(float) + sizeof(F_CARRY_T));

	for (i = 0; i < *len; i++) {

		sIAdd1(fold,
			(float*)inoutvec, (F_CARRY_T*)(inoutvec + array_of_addresses[1]), 1,
			(float*)invec, (F_CARRY_T*)(invec + array_of_addresses[1]), 1
		);
		sIRenorm1(fold,
			(float*)inoutvec, (F_CARRY_T*)(inoutvec + array_of_addresses[1]), 1
		);

		invec    = (char*)invec + size;
		inoutvec = (char*)inoutvec + size;
	}
}

void cIAdd_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
) {
	int i;
	int* pfold;
	int fold;
	int size;

	int NI, NA, ND, combiner;
	int array_of_integers[2];
	MPI_Aint array_of_addresses[2];
	MPI_Datatype array_of_datatypes[2];

	MPI_Type_get_envelope(*datatype, &NI, &NA, &ND, &combiner);
	MPI_Type_get_contents(*datatype, NI, NA, ND,
		array_of_integers, array_of_addresses, array_of_datatypes);

	MPI_Type_size(*datatype, &size);

	fold = size / (sizeof(float complex) + 2 * sizeof(F_CARRY_T));

	for (i = 0; i < *len; i++) {

		cIAdd1(fold,
			(float complex*)inoutvec,
			(F_CARRY_T*)(inoutvec + array_of_addresses[1]), 1,
			(float complex*)invec,
			(F_CARRY_T*)(invec + array_of_addresses[1]), 1
		);
		cIRenorm1(fold,
			(float complex*)inoutvec,
			(F_CARRY_T*)(inoutvec + array_of_addresses[1]), 1
		);

		invec    = (char*)invec + size;
		inoutvec = (char*)inoutvec + size;
	}
}

void sINrm2_MPI(
	void *invec, void* inoutvec,
	int* len, MPI_Datatype *datatype
) {
	int i, j;
	float* pin;
	float* pout;
	int* pfold;
	int fold;
	float scale1, scale2;

	int size;

	int NI, NA, ND, combiner;
	int array_of_integers[2];
	MPI_Aint array_of_addresses[2];
	MPI_Datatype array_of_datatypes[2];

	MPI_Type_get_envelope(*datatype, &NI, &NA, &ND, &combiner);
	MPI_Type_get_contents(*datatype, NI, NA, ND,
		array_of_integers, array_of_addresses, array_of_datatypes);

	MPI_Type_size(*datatype, &size);

	fold = (size - sizeof(float)) / (sizeof(float) + sizeof(F_CARRY_T));

	for (i = 0; i < *len; i++) {
		pin  = (float*) invec;
		pout = (float*) inoutvec;

		scale1 = pout[0];
		scale2 = pin[0];

		if (scale1 < scale2) {
			scale1 = scale1 / scale2;
			scale1 = scale1 * scale1;
			for (j = 1; j < 1+fold; j++) {
				pout[j] *= scale1;
			}
			pout[0] = scale2;
		}
		if (scale1 > scale2) {
			scale2 = scale2 / scale1;
			scale2 = scale2 * scale2;
			for (j = 1; j < 1+fold; j++) {
				pin[j] *= scale2;
			}
			pin[0] = scale1;
		}

		sIAdd1(fold,
			pout + 1, (F_CARRY_T*)(inoutvec + array_of_addresses[1]), 1,
			pin  + 1, (F_CARRY_T*)(invec + array_of_addresses[1]), 1
		);

		invec    = (char*)invec + size;
		inoutvec = (char*)inoutvec + size;

	}
}

int sIMPICreate_(MPI_Datatype* newtype, int size1, int size2) {
	int array_of_blocklengths[2];
	MPI_Aint array_of_displacements[2];
	MPI_Datatype array_of_types[2];

	array_of_blocklengths[0] = size1;
	array_of_blocklengths[1] = size2;

	array_of_displacements[0] = 0;
	array_of_displacements[1] = size1 * sizeof(float);

	array_of_types[0] = MPI_FLOAT;
#if F_CARRY_OP == 0
	array_of_types[1] = MPI_FLOAT;
#elif F_CARRY_OP == 1
	array_of_types[1] = MPI_DOUBLE;
#elif F_CARRY_OP == 2
	array_of_types[1] = MPI_INT;
#elif F_CARRY_OP == 3
	array_of_types[1] = MPI_LONG;
#endif

	int ret = MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,array_of_types, newtype);
	if (ret == MPI_SUCCESS) {
   		ret = MPI_Type_commit( newtype );
	}

	if (ret == MPI_SUCCESS)
		return MPI_SUCCESS;

	if (ret == MPI_ERR_ARG) {
		printf("[%s.%d] Invalid argument.", __FILE__, __LINE__);
	}
	if (ret == MPI_ERR_TYPE) {
		printf("[%s.%d] Invalid datatype argument.", __FILE__, __LINE__);
	}

	MPI_Abort(MPI_COMM_WORLD, ret);

	return ret;
}

int sIMPICreate(int fold, MPI_Datatype* newtype) {
	return sIMPICreate_(newtype, fold, fold);
}

int sIMPIScaleCreate(int fold, MPI_Datatype* newtype) {
	return sIMPICreate_(newtype, 1+fold, fold);
}

int cIMPICreate(int fold, MPI_Datatype* newtype) {
	int array_of_blocklengths[2];
	MPI_Aint array_of_displacements[2];
	MPI_Datatype array_of_types[2];

	array_of_blocklengths[0] = fold;
	array_of_blocklengths[1] = 2 * fold;

	array_of_displacements[0] = 0;
	array_of_displacements[1] = fold * sizeof(float complex);

	array_of_types[0] = MPI_COMPLEX;
#if F_CARRY_OP == 0
	array_of_types[1] = MPI_FLOAT;
#elif F_CARRY_OP == 1
	array_of_types[1] = MPI_DOUBLE;
#elif F_CARRY_OP == 2
	array_of_types[1] = MPI_INT;
#elif F_CARRY_OP == 3
	array_of_types[1] = MPI_LONG;
#endif

	int ret = MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,array_of_types, newtype);
	if (ret == MPI_SUCCESS) {
   		ret = MPI_Type_commit( newtype );
	}

	if (ret == MPI_SUCCESS)
		return MPI_SUCCESS;

	if (ret == MPI_ERR_ARG) {
		printf("[%s.%d] Invalid argument.", __FILE__, __LINE__);
	}
	if (ret == MPI_ERR_TYPE) {
		printf("[%s.%d] Invalid datatype argument.", __FILE__, __LINE__);
	}

	MPI_Abort(MPI_COMM_WORLD, ret);

	return ret;
}

