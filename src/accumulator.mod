//-==============================-//
//  Accumulator: to improve performance
//  of single aggregation by using
//  Iblas functions
//-==============================-//


typedef struct dIAccumulator {
	Idouble v;
	int counter;
	int size;
	double* BUFFER;
} dIAccum;

typedef struct sIAccumulator {
	Ifloat v;
	int counter;
	int size;
	float* BUFFER;
} sIAccum;

typedef struct zIAccumulator {
	dIcomplex v;
	int counter;
	int size;
	double complex* BUFFER;
} zIAccum;

typedef struct cIAccumulator {
	sIcomplex v;
	int counter;
	int size;
	float complex* BUFFER;
} cIAccum;

extern void dIAccReset(dIAccum *acc);
extern void dIAccInit(dIAccum *acc, int BUFFER_SIZE);
extern void dIAccDestroy(dIAccum *acc);
extern void dIAccumulate(dIAccum *acc, double x);
extern void dIAccumulates(dIAccum *acc, int n, double* x, int inc);
extern double dIAccExtract(dIAccum *acc);

extern void sIAccReset(sIAccum *acc);
extern void sIAccInit(sIAccum *acc, int BUFFER_SIZE);
extern void sIAccDestroy(sIAccum *acc);
extern void sIAccumulate(sIAccum *acc, float x);
extern void sIAccumulates(sIAccum *acc, int n, float* x, int inc);
extern float sIAccExtract(sIAccum *acc);

extern void zIAccReset(zIAccum *acc);
extern void zIAccInit(zIAccum *acc, int BUFFER_SIZE);
extern void zIAccDestroy(zIAccum *acc);
extern void zIAccumulate(zIAccum *acc, double complex x);
extern void zIAccumulates(zIAccum *acc, int n, double complex* x, int inc);
extern double complex zIAccExtract(zIAccum *acc);

extern void cIAccReset(cIAccum *acc);
extern void cIAccInit(cIAccum *acc, int BUFFER_SIZE);
extern void cIAccDestroy(cIAccum *acc);
extern void cIAccumulate(cIAccum *acc, float complex x);
extern void cIAccumulates(cIAccum *acc, int n, float complex* x, int inc);
extern float complex cIAccExtract(cIAccum *acc);

