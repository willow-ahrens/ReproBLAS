
/*******************************************/
/* WRAPPER FOR DOUBLE PRECISION            */
/*******************************************/

#define dIprint(X) dIprint1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define zIprint(X) zIprint1(DEFAULT_FOLD, \
	(double complex*)((X).m), (double complex*)((X).c), 1)

//====================================//
// ADDITION
//====================================//


// ADDING A NATIVE FP TO AN INDEXED FP
#define dIUpdate_(X,Y) dIUpdate1(DEFAULT_FOLD, 0, (X).m, (X).c, 1, fabs(Y))
#define dIAddd_(X,Y)   dIAddd1(DEFAULT_FOLD, (X).m, 1, Y)
#define dIRenorm_(X)   dIRenorm1(DEFAULT_FOLD, (X).m, (X).c, 1)

#define zIUpdate_(X,Y) zIUpdate1(DEFAULT_FOLD, 0,	\
		(double complex*)((X).m),(double complex*)((X).c),1,Y)
#define zIAddz_(X,Y) zIAddz1(DEFAULT_FOLD, (double complex*)((X).m), 1, Y)
#define zIRenorm_(X) zIRenorm1(DEFAULT_FOLD,	\
		(double complex*)((X).m), (double complex*)((X).c), 1)

//====================================//
// NEGATION
//====================================//
#define dINeg(X) dINeg1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define zINeg(X) zINeg1(DEFAULT_FOLD, (double complex*)(X).m, (double complex*)(X).c, 1)

//====================================//
// CONVERSION
//====================================//
#define Iconv2d(X)  Iconv2d1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define Iconv2z(X) Iconv2z1(DEFAULT_FOLD,	\
	(double complex*)(X).m, (double complex*)(X).c, 1)

