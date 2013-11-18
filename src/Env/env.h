#ifndef _RBLAS_ENV__H
#define _RBLAS_ENV__H

#define DGEMM_NB_M 32
#define DGEMM_NB_N 32
#define DGEMM_NB_K 32

#define DAXPY_NB   128

#define Env_Get(FCT,PARAM) FCT##_##PARAM
#define Env_Set(FCT,PARAM) {}

#endif
