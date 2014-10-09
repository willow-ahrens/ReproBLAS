/*
 *  Created   13/10/25   H.D. Nguyen
 */

#ifndef _REPRODUCIBLE_COMMON_UTIL__H_
#define _REPRODUCIBLE_COMMON_UTIL__H_

#include <stdint.h>

#ifndef MAX
#define MAX(A,B) (A>B?A:B)
#endif

#ifndef MIN
#define MIN(A,B) (A<B?A:B)
#endif

#define RELERR(A,B) (A==B)?0:fabs(A - B)/MAX(fabs(A),fabs(B))

#define DIMM1(TRANS,M,N) ((TRANS=='T'||TRANS=='t')?N:M)
#define DIMM2(TRANS,M,N) ((TRANS=='T'||TRANS=='t')?M:N)

typedef union l_double_ {
	double d;
	long   l;
} l_double;

typedef union i_float_ {
	float f;
	int   i;
} i_float;

void dmovv(int N, double* dst, int incD, double* src, int incS);
void dmovm(int M, int N, double* dst, int ldD, double* src, int ldS);
void dsetv(int N, double* dst, int incD, double x);
void dsetm(int M, int N, double* dst, int ldD, double x);
void dshowm(double* M, int ldM, int frM, int toM, int frN, int toN);

// CHECK IF AN ADDRESS IS ALIGNED TO 16-BYTE BOUNDARE OR NOT
#define IS_ALIGNED(ADDR, ALIGN) ((((uintptr_t)ADDR) & (ALIGN - 1)) == 0)
#define IS_UNALIGNED(ADDR, ALIGN) ((((uintptr_t)ADDR) & (ALIGN - 1)) > 0)

#define SCALAR_SPLIT(M,X,Q)	\
	Q = (M + X);		\
	Q = (Q - M);		\
	X = (X - Q);

#define SCALAR_EXTRACT(M,X)	\
	X = (M + X);		\
	X = (X - M);

#define SCALAR_SPLIT_SIMD(M,X,Q)	\
	Q = _mm_add_pd(M, X);		\
	Q = _mm_sub_pd(Q, M);		\
	X = _mm_sub_pd(X, Q);

#define SCALAR_EXTRACT_SIMD(M,X)	\
	X = _mm_add_pd(M, X);		\
	X = _mm_sub_pd(X, M);

#define SCALAR_ACC(M,X,Q)	\
	Q = M;				\
	M = (M + X);		\
	Q = (Q - M);		\
	X = X + Q;

#define SCALAR_SPLIT_FNR(M,X,Q,MASK)	\
	X.l |= MASK;		\
	Q = (M + X.d);		\
	Q = (Q - M);		\
	X.d = (X.d - Q);

#define SCALAR_EXTRACT_FNR(M,X,MASK)	\
	X.l |= MASK;		\
	X.d = (M + X);		\
	X.d = (X - M);

#define SCALAR_SPLIT_SIMD_FNR(M,X,Q,MASK)	\
	X = _mm_and_pd(X, MASK);	\
	Q = _mm_add_pd(M, X);		\
	Q = _mm_sub_pd(Q, M);		\
	X = _mm_sub_pd(X, Q);

#define SCALAR_EXTRACT_SIMD_FNR(M,X)	\
	X = _mm_and_pd(X, MASK);	\
	X = _mm_add_pd(M, X);		\
	X = _mm_sub_pd(X, M);

#define SSE_ABS_MASKS(MASK)		\
{	\
	__m128 r1__;					\
	r1__ = _mm_set1_ps(1);			\
	MASK = _mm_set1_ps(-1);			\
	MASK = _mm_xor_ps(MASK, r1__);	\
	r1__ = _mm_cmpeq_ps(r1__, r1__);\
	MASK = _mm_xor_ps(MASK, r1__);	\
}

#define SSE_CONJ_MASKS(MASK)		\
{	\
	__m128 r1__;					\
	r1__ = _mm_set_ps(1, 0, 1, 0);			\
	MASK = _mm_set_ps(-1, 0, -1, 0);			\
	MASK = _mm_xor_ps(MASK, r1__);	\
}

#define SSE_CONJ_MASKD(MASK)		\
{	\
	__m128d r1__;					\
	r1__ = _mm_set_pd(1, 0);			\
	MASK = _mm_set_pd(-1, 0);			\
	MASK = _mm_xor_pd(MASK, r1__);	\
}

#define AVX_CONJ_MASKS(MASK)		\
{	\
	__m256 r1__;					\
	r1__ = _mm256_set_ps(1, 0, 1, 0, 1, 0, 1, 0);			\
	MASK = _mm256_set_ps(-1, 0, -1, 0, -1, 0, -1, 0);			\
	MASK = _mm256_xor_ps(MASK, r1__);	\
}

#define AVX_CONJ_MASKD(MASK)		\
{	\
	__m256d r1__;					\
	r1__ = _mm256_set_pd(1, 0, 1, 0);			\
	MASK = _mm256_set_pd(-1, 0, -1, 0);			\
	MASK = _mm256_xor_pd(MASK, r1__);	\
}

#define SSE_NCONJ_MASKS(MASK)		\
{	\
	__m128 r1__;					\
	r1__ = _mm_set_ps(0, 1, 0, 1);			\
	MASK = _mm_set_ps(0, -1, 0, -1);			\
	MASK = _mm_xor_ps(MASK, r1__);	\
}

#define SSE_NCONJ_MASKD(MASK)		\
{	\
	__m128d r1__;					\
	r1__ = _mm_set_pd(0, 1);			\
	MASK = _mm_set_pd(0, -1);			\
	MASK = _mm_xor_pd(MASK, r1__);	\
}

#define AVX_NCONJ_MASKS(MASK)		\
{	\
	__m256 r1__;					\
	r1__ = _mm256_set_ps(0, 1, 0, 1, 0, 1, 0, 1);			\
	MASK = _mm256_set_ps(0, -1, 0, -1, 0, -1, 0, -1);			\
	MASK = _mm256_xor_ps(MASK, r1__);	\
}

#define AVX_NCONJ_MASKD(MASK)		\
{	\
	__m256d r1__;					\
	r1__ = _mm256_set_pd(0, 1, 0, 1);			\
	MASK = _mm256_set_pd(0, -1, 0, -1);			\
	MASK = _mm256_xor_pd(MASK, r1__);	\
}

#define SSE_ABS_MASKD(MASK)		\
{	\
	__m128d r1__;					\
	r1__ = _mm_set1_pd(1);			\
	MASK = _mm_set1_pd(-1);			\
	MASK = _mm_xor_pd(MASK, r1__);	\
	r1__ = _mm_cmpeq_pd(r1__, r1__);\
	MASK = _mm_xor_pd(MASK, r1__);	\
}

#define AVX_ABS_MASKS(MASK)		\
{	\
	__m256 r1__;					\
	r1__ = _mm256_set1_ps(1);			\
	MASK = _mm256_set1_ps(-1);			\
	MASK = _mm256_xor_ps(MASK, r1__);	\
	r1__ = _mm256_cmp_ps(r1__, r1__, 0);\
	MASK = _mm256_xor_ps(MASK, r1__);	\
}

#define AVX_ABS_MASKD(MASK)		\
{	\
	__m256d r1__;					\
	r1__ = _mm256_set1_pd(1);			\
	MASK = _mm256_set1_pd(-1);			\
	MASK = _mm256_xor_pd(MASK, r1__);	\
	r1__ = _mm256_cmp_pd(r1__, r1__, 0);\
	MASK = _mm256_xor_pd(MASK, r1__);	\
}

#define SSE_BLP_MASKS(MASK)		\
{	\
	__m128 r1__;					\
	MASK = _mm_set1_ps(1.0);		\
	r1__ = _mm_set1_ps(1.0 + (FLT_EPSILON * 1.0001));	\
	MASK = _mm_xor_ps(MASK, r1__);	\
}

#define SSE_BLP_MASKD(MASK)		\
{	\
	__m128d r1__;					\
	MASK = _mm_set1_pd(1.0);		\
	r1__ = _mm_set1_pd(1.0 + (DBL_EPSILON * 1.0001));	\
	MASK = _mm_xor_pd(MASK, r1__);	\
}

#define AVX_BLP_MASKS(MASK)		\
{	\
	__m256 r1__;					\
	MASK = _mm256_set1_ps(1.0);		\
	r1__ = _mm256_set1_ps(1.0 + (FLT_EPSILON * 1.0001));	\
	MASK = _mm256_xor_ps(MASK, r1__);	\
}

#define AVX_BLP_MASKD(MASK)		\
{	\
	__m256d r1__;					\
	MASK = _mm256_set1_pd(1.0);		\
	r1__ = _mm256_set1_pd(1.0 + (DBL_EPSILON * 1.0001));	\
	MASK = _mm256_xor_pd(MASK, r1__);	\
}

#define AVX_0001_MASK(MASK)		\
{ \
  MASK = _mm256_set_epi64x(0, 0, 0, 0xFFFFFFFFFFFFFFFFUL);\
}

#define SET_DAZ_FLAG \
	unsigned int _old_csr_ = _mm_getcsr();	\
	unsigned int _new_csr_ = _old_csr_ | 0x8040; \
	if (_new_csr_ != _old_csr_) _mm_setcsr(_new_csr_);

#define RESET_DAZ_FLAG \
	if (_new_csr_ != _old_csr_) _mm_setcsr(_old_csr_);

#endif
