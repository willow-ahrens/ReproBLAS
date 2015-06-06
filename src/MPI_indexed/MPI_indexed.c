/*
 *  Created   13/10/25   H.D. Nguyen
 */

#include "MPI_indexed.h"
#include "MPI_dIndexed.h"
#include "MPI_sIndexed.h"
#include "indexed.h"
#include "../../config.h"
#include <stdio.h>
#include <complex.h>

MPI_Datatype MPI_IDOUBLE;
MPI_Datatype MPI_IFLOAT;
MPI_Datatype MPI_IDOUBLE_SCALE; // FOR SECOND NORM 
MPI_Datatype MPI_IFLOAT_SCALE;  // FOR SECOND NORM
MPI_Datatype MPI_ICOMPLEX;
MPI_Datatype MPI_IDOUBLE_COMPLEX;

MPI_Op MPI_RSUM  = MPI_OP_NULL;
MPI_Op MPI_RNRM2 = MPI_OP_NULL;

static int MPI_INDEXED_INIT = -1;

int MPI_IFLOAT_Create(MPI_Datatype* newtype) {
	int array_of_blocklengths[2];
	MPI_Aint array_of_displacements[2];
	MPI_Datatype array_of_types[2];
	I_float sample;
	MPI_Aint base_address;
	MPI_Aint mantissa_address;
	MPI_Aint carry_address;

	MPI_Get_address(&sample, &base_address);
	MPI_Get_address(&sample.m, &mantissa_address);
	MPI_Get_address(&sample.c, &carry_address);

	array_of_blocklengths[0] = (carry_address - mantissa_address) / sizeof(float);
	array_of_blocklengths[1] = DEFAULT_FOLD;

	array_of_displacements[0] = mantissa_address - base_address;
	array_of_displacements[1] = carry_address - base_address;

	array_of_types[0] = MPI_FLOAT;
	array_of_types[1] = MPI_FLOAT;

	int ret = MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,array_of_types, newtype);
	if (ret == MPI_SUCCESS) {
   		ret = MPI_Type_commit( newtype );
	}

	if (ret == MPI_SUCCESS)
		return MPI_SUCCESS;

	if (ret == MPI_ERR_ARG) {
		fprintf(stderr, "[%s.%d] Invalid argument.", __FILE__, __LINE__);
	}
	if (ret == MPI_ERR_TYPE) {
		fprintf(stderr, "[%s.%d] Invalid datatype argument.", __FILE__, __LINE__);
	}

	MPI_Abort(MPI_COMM_WORLD, ret);

	return ret;
}

int RMPI_Init() {
	if (MPI_INDEXED_INIT == 0)
		return MPI_INDEXED_INIT;

	// CREATE DATA TYPE
	int status;
	int array_of_blocklengths[2];
	MPI_Aint array_of_displacements[2];
	MPI_Datatype array_of_types[2];

	I_float  sample;
	I_float_Complex  csample;
	MPI_Aint base_address;
	MPI_Aint mantissa_address;
	MPI_Aint carry_address;

	// I_DOUBLE
   	status = dIMPICreate(DEFAULT_FOLD, &MPI_IDOUBLE );
	if (status != MPI_SUCCESS) {
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}

	// I_DOUBLE_COMPLEX
   	status = zIMPICreate(DEFAULT_FOLD, &MPI_IDOUBLE_COMPLEX );
	if (status != MPI_SUCCESS) {
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}

	// Double scale
   	dIMPIScaleCreate(DEFAULT_FOLD, &MPI_IDOUBLE_SCALE );

	// I_FLOAT
	MPI_Get_address(&sample, &base_address);
	MPI_Get_address(&sample.m, &mantissa_address);
	MPI_Get_address(&sample.c, &carry_address);

	array_of_blocklengths[0] = (carry_address - mantissa_address) / sizeof(float);
	array_of_blocklengths[1] = DEFAULT_FOLD;

	array_of_displacements[0] = mantissa_address - base_address;
	array_of_displacements[1] = carry_address - base_address;

	array_of_types[0] = MPI_FLOAT;
	array_of_types[1] = MPI_FLOAT;

	status = MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,array_of_types, &MPI_IFLOAT);
	if (status == MPI_SUCCESS) {
   		status = MPI_Type_commit( &MPI_IFLOAT );
	}

	if (status == MPI_ERR_ARG) {
		fprintf(stderr, "[%s.%d] Invalid argument.", __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}
	if (status == MPI_ERR_TYPE) {
		fprintf(stderr, "[%s.%d] Invalid datatype argument.", __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}

	// I_FLOAT_SCALE
	array_of_blocklengths[0] += 1;
	array_of_displacements[1] += sizeof(float);

	status = MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,array_of_types, &MPI_IFLOAT_SCALE);
	if (status == MPI_SUCCESS) {
   		status = MPI_Type_commit( &MPI_IFLOAT_SCALE );
	}

	if (status == MPI_ERR_ARG) {
		fprintf(stderr, "[%s.%d] Invalid argument.", __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}
	if (status == MPI_ERR_TYPE) {
		fprintf(stderr, "[%s.%d] Invalid datatype argument.", __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}

	// I_COMPLEX
	MPI_Get_address(&csample, &base_address);
	MPI_Get_address(&csample.m, &mantissa_address);
	MPI_Get_address(&csample.c, &carry_address);

	array_of_blocklengths[0] = (carry_address - mantissa_address) / sizeof(float complex);
	array_of_blocklengths[1] = 2 * DEFAULT_FOLD;

	array_of_displacements[0] = mantissa_address - base_address;
	array_of_displacements[1] = carry_address - base_address;

	array_of_types[0] = MPI_COMPLEX;

	status = MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,array_of_types, &MPI_ICOMPLEX);
	if (status == MPI_SUCCESS) {
   		status = MPI_Type_commit( &MPI_ICOMPLEX);
	}

	if (status == MPI_ERR_ARG) {
		fprintf(stderr, "[%s.%d] Invalid argument.", __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}
	if (status == MPI_ERR_TYPE) {
		fprintf(stderr, "[%s.%d] Invalid datatype argument.", __FILE__, __LINE__);
		MPI_Abort(MPI_COMM_WORLD, status);
		return status;
	}

	// CREATE OPERATION
	MPI_Op_create(&IAdd_MPI, 1, &MPI_RSUM);
	MPI_Op_create(&INrm2_MPI, 1, &MPI_RNRM2);

	MPI_INDEXED_INIT = 0;
	return 0;
}

void IAdd_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype) {

	// PREDEFINED DATATYPES
	if (*datatype == MPI_IDOUBLE) {
		dIAdd_MPI(in, inout, len, datatype);
		return;
	}

	if (*datatype == MPI_IFLOAT) {
		sIAdd_MPI(in, inout, len, datatype);
		return;
	}

	if (*datatype == MPI_ICOMPLEX) {
		cIAdd_MPI(in, inout, len, datatype);
		return;
	}

	if (*datatype == MPI_IDOUBLE_COMPLEX) {
		zIAdd_MPI(in, inout, len, datatype);
		return;
	}

	// CUSTOM DATATYPES WITH DIFFERENT FOLD
	int NI, NA, ND, combiner;
	int array_of_integers[3];
	MPI_Aint array_of_addresses[3];
	MPI_Datatype array_of_datatypes[3];

	MPI_Type_get_envelope(*datatype, &NI, &NA, &ND, &combiner);
	MPI_Type_get_contents(*datatype, NI, NA, ND, array_of_integers, array_of_addresses, array_of_datatypes);
	
	if (array_of_datatypes[0] == MPI_DOUBLE) {
		dIAdd_MPI(in, inout, len, datatype);
		return;
	}

	if (array_of_datatypes[0] == MPI_FLOAT) {
		sIAdd_MPI(in, inout, len, datatype);
		return;
	}

	if (array_of_datatypes[0] == MPI_COMPLEX) {
		cIAdd_MPI(in, inout, len, datatype);
		return;
	}

	if (array_of_datatypes[0]== MPI_DOUBLE_COMPLEX) {
		zIAdd_MPI(in, inout, len, datatype);
		return;
	}

}

void INrm2_MPI(void *in, void* inout,
	int* len, MPI_Datatype *datatype) {

	// PREDEFINED DATATYPES
	if (*datatype == MPI_IDOUBLE_SCALE) {
		dINrm2_MPI(in, inout, len, datatype);
		return;
	}

	if (*datatype == MPI_IFLOAT_SCALE) {
		sINrm2_MPI(in, inout, len, datatype);
		return;
	}

	// CUSTOM DATATYPES WITH DIFFERENT FOLD
	int NI, NA, ND, combiner;
	int array_of_integers[3];
	MPI_Aint array_of_addresses[3];
	MPI_Datatype array_of_datatypes[3];

	MPI_Type_get_envelope(*datatype, &NI, &NA, &ND, &combiner);
	MPI_Type_get_contents(*datatype, NI, NA, ND, array_of_integers, array_of_addresses, array_of_datatypes);
	
	if (array_of_datatypes[0] == MPI_DOUBLE) {
		dINrm2_MPI(in, inout, len, datatype);
		return;
	}

	if (array_of_datatypes[0] == MPI_FLOAT) {
		sINrm2_MPI(in, inout, len, datatype);
		return;
	}
}
