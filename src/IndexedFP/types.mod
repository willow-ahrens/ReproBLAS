#include <complex.h>
#include <math.h>

//====================================//
//  DEFINE TYPES                      //
//====================================//

#if (F_CARRY_OP == 1)
#define F_CARRY_T  double
#elif (F_CARRY_OP == 0)
#define F_CARRY_T  float
#elif (F_CARRY_OP == 2)
#define F_CARRY_T  int
#elif (F_CARRY_OP == 3)
#define F_CARRY_T  int
#endif

typedef struct Idouble_ {
	double m[DEFAULT_FOLD];
	double c[DEFAULT_FOLD];
} I_double;

#define Idouble I_double
#define Ifloat  I_float
#define dIcomplex I_double_Complex
#define sIcomplex I_float_Complex

typedef struct dIComplex_ {
	double m[2*DEFAULT_FOLD];
	double c[2*DEFAULT_FOLD];
} I_double_Complex;

typedef struct Ifloat_ {
	float  m[DEFAULT_FOLD];
	F_CARRY_T c[DEFAULT_FOLD];
} I_float;

typedef struct sIComplex_{
	float m[2*DEFAULT_FOLD];
	F_CARRY_T c[2*DEFAULT_FOLD];
} I_float_Complex;

#define zISize(K) (2*K*sizeof(double complex))
#define dISize(K) (2*K*sizeof(double))
#define sISize(K) (K*(sizeof(float)+sizeof(F_CARRY_T)))
#define cISize(K) (K*sizeof(float complex)+2*K*sizeof(F_CARRY_T))

// SET ZERO
#define ISetZero_(K,M,C) {	\
	int i;	\
	for (i = 0; i < K; i++) { \
		M[i] = 0;	\
		C[i] = 0;	\
	}	\
}
#define dISetZero(I) ISetZero_(DEFAULT_FOLD, (I).m, (I).c)
#define sISetZero(I) ISetZero_(DEFAULT_FOLD, (I).m, (I).c)
#define zISetZero(I) ISetZero_(2*DEFAULT_FOLD, (I).m, (I).c)
#define cISetZero(I) ISetZero_(2*DEFAULT_FOLD, (I).m, (I).c)

// DI = SI
#define ISet_(K,DI,SI) {	\
	int i;	\
	for (i = 0; i <  K; i++) {	\
		(DI).m[i] = (SI).m[i]; \
		(DI).c[i] = (SI).c[i]; \
	}	\
}
#define dISet(DI,SI) ISet_(DEFAULT_FOLD, DI, SI)
#define sISet(DI,SI) ISet_(DEFAULT_FOLD, DI, SI)
#define zISet(DI,SI) ISet_(2*DEFAULT_FOLD, DI, SI)
#define cISet(DI,SI) ISet_(2*DEFAULT_FOLD, DI, SI)

// ZI = DI
#define zdISet_(K,ZI,DI) {	\
	ISetZero_(2 * K, (ZI).m, (ZI).c)\
	int i;	\
	for (i = 0; i <  K; i++) {	\
		(ZI).m[2 * i] = (DI).m[i]; \
		(ZI).c[2 * i] = (DI).c[i]; \
	}	\
}
#define zdISet(ZI, DI) zdISet_(DEFAULT_FOLD, ZI, DI)
#define csISet(CI, SI) zdISet_(DEFAULT_FOLD, CI, SI)
