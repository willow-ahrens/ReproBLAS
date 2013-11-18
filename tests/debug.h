#ifndef UTIL_DEBUG__H
#define UTIL_DEBUG__H

#include <string.h>
#include "../src/Common/Common.h"
#include "../src/types.h"
#include <complex.h>

#define CALL_RSDOT
#define CALL_RSASUM
#define CALL_RSSUM
#define CALL_RSNRM2

#define CALL_PRSASUM
#define CALL_PRSSUM
#define CALL_PRSDOT
#define CALL_PRSNRM2

#define ASUM_BIT 0
#define RASUM_BIT 1
#define RSUM_BIT 2
#define NRM2_BIT 3
#define RNRM2_BIT 4
#define DOT_BIT 5
#define RDOT_BIT 6

#include "blas_intf.h"

extern double TwoProd(double, double, double*);
extern double TwoSum(double, double, double*);

double read_timer();
int    find_option( int argc, char **argv, const char *option );
int    read_int( int argc, char **argv, const char *option, int default_value );
void show_opt(char* name, char* desc);
double read_double( int argc, char **argv, const char *option, double default_value );
float read_float( int argc, char **argv, const char *option, float default_value );
double ddiff(int N, double* v1, double* v2, int* ind);
double derrs(int N, double* v1, double* v2, int* ind, double* abs);
void drandm(double* A, int M, int N, int lda, double from, double to);
void drandomseed();
void dgenvec_(int N, double* v, int inc, int type, double factor);
void dgenvec (int N, double* v, int type, double factor);
void sgenvec_(int N, float*  v, int inc,  int type, float  factor);
void sgenvec (int N, float*  v, int type, float  factor);
void dsort_bubble(int N, double* v, int ord);
void dsort_merge (int N, double* v, int ord);
void ssort_bubble(int N, float* v, int ord);
void ssort_merge (int N, float* v, int ord);
double dsum(int N, double* v);
float  ssum(int N, float*  v);
double pdsum(int N, double* v, int procs);
int read_irange(int argc, char **argv, const char *option,
	int* pstart, int* pstop, int* pstep);
int rsign();
void dcheck_reproducibility(int n, double* x, double* y, int* status, double* result);
void scheck_reproducibility(int n, float* x, float* y, int* status, float* ref);
void zcheck_reproducibility(int n, dcomplex* x, int incx, dcomplex* y, int incy,
	int* status, dcomplex* ref);
void ccheck_reproducibility(int n, float complex* x, int incx, float complex* y, int incy,
	int* status, float complex* ref);

int fileExists(const char *fname);
void readFile(const char* fname, int* len, void** data, int eleSize);
void writeFile(const char* fname, int len, void* data, int eleSize);


void dmerge_(
	int N,
	int N1, double* v1, int* ptr_id1,
	int N2, double* v2, int* ptr_id2,
	double* v,
	int ord
);
void dmerge_sort_(
	int N1, double* v1,
	int N2, double* v2,
	int ord,
	double* buffer
);
void dsort_merge_(
	int N,
	double* v,
	int ord,
	double* buffer
);

#define WHERESTR  "[F %s, L%3d]:"
#define WHEREARG  __FILE__, __LINE__
#define DEBUGPRINT2(...)  printf(__VA_ARGS__)
#define debug(_fmt, ...)                                \
 DEBUGPRINT2(WHERESTR _fmt, WHEREARG, __VA_ARGS__)
#define Debug(str)  DEBUGPRINT2(WHERESTR "%s", WHEREARG, str)
#define MARK DEBUGPRINT2(WHERESTR, WHEREARG)

#define BREAKDOWN_N 4

void timing_breakdown(double* timing); 
void clock_start();
void tictoc(int);
void dreverse_(int N, double* x, int inc);
void dreverse (int N, double* x);
void dshuffle_(int N, double* x, int inc);
void dshuffle (int N, double* x);
void sreverse_(int N, float* x, int inc);
void sreverse (int N, float* x);
void sshuffle_(int N, float* x, int inc); 
void sshuffle (int N, float* x); 

void creverse (int N, float complex* x, int inc);
void zreverse (int N, dcomplex* x, int inc);
////////////////////////////////
// SORTING

void ccheck_reproducibility(int n, float complex* x, int inc, float complex* y, int incy,
	int* status, float complex* ref);
void zcheck_reproducibility(int n, dcomplex* x, int inc, dcomplex* y, int incy,
	int* status, dcomplex* ref);

void csort_bubble(int N, float complex* v, int inc, int ord);
void csort_merge(
	int N,
	float complex* v, int inc,
	int ord
);


void zsort_bubble(int N, dcomplex* v, int inc, int ord);
void zsort_merge(
	int N,
	dcomplex* v, int inc,
	int ord
);

#endif
