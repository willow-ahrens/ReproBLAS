
/*******************************************/
/* WRAPPER FOR SINGLE PRECISION            */
/*******************************************/

#define sIprint(X) sIprint1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define cIprint(X) cIprint1(DEFAULT_FOLD, (float complex*)(X).m, (X).c, 1)

//====================================//
// ADDITION
//====================================//


// ADDING A NATIVE FP TO AN INDEXED FP
#define sIUpdate_(X,Y) sIUpdate1(DEFAULT_FOLD, 0, fabs(Y), (X).m, (X).c, 1)
#define sIAddf_(X,Y)   sIAddf1(DEFAULT_FOLD, (X).m, 1, Y)
#define sIRenorm_(X)   sIRenorm1(DEFAULT_FOLD, (X).m, (X).c, 1)

#define cIUpdate_(X,Y) cIUpdate1(DEFAULT_FOLD,0,(float complex*)((X).m),(X).c,1,Y)
#define cIAddc_(X,Y) cIAddc1(DEFAULT_FOLD, (float complex*)((X).m), 1, Y)
#define cIRenorm_(X) cIRenorm1(DEFAULT_FOLD,(float complex*)((X).m), (X).c, 1)

//====================================//
// NEGATION
//====================================//
#define sINeg(X) dINeg1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define cINeg(X) zINeg1(DEFAULT_FOLD, (float complex*)(X).m, (X).c, 1)

//====================================//
// CONVERSION
//====================================//
#define Iconv2f(X) Iconv2f1(DEFAULT_FOLD, (X).m, (X).c, 1)
#define Iconv2c(X) Iconv2c1(DEFAULT_FOLD, (float complex*)(X).m, (X).c, 1)

