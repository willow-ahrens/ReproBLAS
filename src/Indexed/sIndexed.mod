
extern int sICapacity();
extern int sIWidth();

extern void sIprint1(int n, float* x, F_CARRY_T* carry, int inc);
extern void cIprint1(int n, float complex* x, F_CARRY_T* carry, int inc);

// ADD A FLOAT TO AN INDEXED FP
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    Y     : FP TO BE ADDED TO X
// [IN/OUTPUT]
//    X     : INDEXED FP, AT RETURN X = X + Y
extern void sIAddf1(int fold, float* x, int inc, float y);
extern void sIAddf(I_float* X, float Y);
extern void cIAddc1(int fold, float complex* x, int inc, float complex Y);
extern void cIAddc(I_float_Complex* X, float complex Y);

// NEGATION
extern void sINeg1(int fold, float* x, F_CARRY_T* c, int inc);
extern void cINeg1(int fold, float complex* x, F_CARRY_T* c, int inc);

// COMPUTE BOUNDARIES BASED ON MAXIMUM ABSOLUTE VALUE
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
//    W     : WIDTH OF EACH BINS
//    MAX   : MAXIMUM ABSOLUTE VALUE OF INPUT VALUES TO BE SUMMED
// [OUTPUT]
//    M     : PRECOMPUTED BOUNDARIES
extern int  sIBoundary_(int fold, int W, float max, float* M, int inc);

// UNIT IN THE FIRST PLACE
extern float ufpf(float x);

// X += Y
extern void sIAdd1(int fold, float* x,
	F_CARRY_T* xc, int incx, float* y, F_CARRY_T* yc, int incy);
extern void sIAdd(I_float* X, I_float Y);

extern void cIAdd1(int fold, float complex* x, F_CARRY_T* xc, int incx,
	float complex* y, F_CARRY_T* yc, int incy);
extern void cIAdd(I_float_Complex* X, I_float_Complex Y);

//====
extern void sIUpdate1(int fold, int W, float y, float* x, F_CARRY_T* C, int inc);

// UPDATE A COMPLEX USING A FLOAT
extern void cIUpdates1(int K, int W,
	float complex* X, F_CARRY_T* C, int INC, float Y);

// UPDATE A COMPLEX USING A COMPLEX
extern void cIUpdate1(int K, int W,
	float complex* X, F_CARRY_T* C,int INC, float complex Y);

// RENORMALIZATION TO AVOID OVERFLOW
// [INPUT]
//    FOLD  : NB OF BINS OF LEADING BITS (NB OF M TO BE COMPUTED)
// [IN/OUTPUT]
//    X     : INDEXED FP
extern void sIRenorm1(int fold, float* X, F_CARRY_T* C, int inc);

extern void cIRenorm1(int fold, float complex* sum, F_CARRY_T* C, int inc);

// CONVERSION FROM INDEXED FORMAT TO FLOAT
extern float   Iconv2f1(int fold, float* m, F_CARRY_T* c, int inc);

extern float complex Iconv2c1(int fold, float complex* m, F_CARRY_T* c, int inc);
extern float complex Iconv2c_(int fold, float complex* rep);

// CONVERT FROM FLOAT TO INDEXED FORMAT
extern void fconv2I1(int fold, int W, float x,
              float* rep, F_CARRY_T* carry, int inc);
extern void cconv2I1(int fold, int W,
              float complex x, float complex* m, F_CARRY_T* c, int inc);

extern I_float         fconv2I(float x);
extern I_float_Complex cconv2I(float complex x);

