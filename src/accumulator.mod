typedef struct dIAccumulator {
	Idouble v;
	int counter;
	int size;
	double* BUFFER;
} dIAccum;

extern void dResetAcc(dIAccum *acc);
extern void dInitAccum(dIAccum *acc, int BUFFER_SIZE);
extern void dDestroyAcc(dIAccum *acc);
extern void dAccumulate(dIAccum *acc, double x);
extern void dAccumulates(dIAccum *acc, int n, double* x, int inc);
extern double dExtractAcc(dIAccum *acc);

