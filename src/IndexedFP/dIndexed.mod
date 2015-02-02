
extern int dIWidth();
extern int dICapacity();
extern void dIprint1(int n, double *x, double* c, int inc);
extern void zIprint1(int n, double complex* x, double complex* carry, int inc);

extern void dIAdd1(int K,
		double* x, double* xc, int incx,
		double* y, double* yc, int incy);
extern void dIAdd(I_double* X, I_double Y);
extern void zIAdd1(int K,
		double complex* x, double complex* xc, int incx,
		double complex* y, double complex* yc, int incy);
extern void zIAdd(I_double_Complex* X, I_double_Complex Y);

// NEGATION
extern void dINeg1(int fold, double* x, double* c, int inc);
extern void zINeg1(int fold, double complex* x, double complex* c, int inc);

// ADD A DOUBLE TO AN INDEXED FP
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    Y     : FP TO BE ADDED TO X
// [IN/OUTPUT]
//    X     : INDEXED FP, AT RETURN X = X + Y
extern void dIAddd1(int fold, double* x, int inc, double y);
extern void dIAddd(I_double* X, double Y);

extern void zIAddz1(int fold, double complex* x, int inc, double complex y);
extern void zIAddz(I_double_Complex* X, double complex Y);

// UPDADATE INDEX FP BY NEW MAXIMUM ABSOLUTE VALUE
extern void dIUpdate1 (int fold, int W, double* x, double* c, int ldx, double y);
extern void zIUpdates1(int fold, int W,
		double complex* x, double complex* c, int ldx, double y);

extern void zIUpdate1 (int fold, int W,
		double complex* x, double complex* c, int ldx,double complex y);

// COMPUTE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    W     : WIDTH OF EACH BINS
//    MAX   : MAXIMUM ABSOLUTE VALUE OF INPUT VALUES TO BE SUMMED
// [OUTPUT]
//    M     : PRECOMPUTED BOUNDARIES
extern int  dIBoundary(int fold, int W, double max, double* M, int inc);

// COMPUTE THE UNIT IN TH FIRST PLACE
extern double ufp(double x);

// RENORMALIZATION TO AVOID OVERFLOW
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
// [IN/OUTPUT]
//    X     : INDEXED FP
extern void dIRenorm1(int fold, double* X, double* C, int inc);
extern void zIRenorm1(int fold, double complex* rep, double complex* c, int inc);

// CONVERT A DOUBLE TO INDEXED FORMAT
extern void dconv2I1(int fold, int W, double x, double* rep, double* C, int inc);
extern I_double dconv2I(double x);

// CONVERT AN INDEXED FP BACK TO DOUBLE
extern double Iconv2d1(int fold, double* m, double* c, int inc);

// CONVERT AN INDEXED COMPLEX TO COMPLEX
extern double complex Iconv2z1(int fold,
		double complex* rep, double complex* carry, int inc);

// CONVERT A DOUBLE COMPLEX TO INDEXED FORMAT
extern void zconv2I1(int fold, int W, double complex X,
		double complex* rep, double complex* C, int inc);
extern I_double_Complex zconv2I(double complex x);

