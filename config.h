// DEFAULT NUMBER OF K
#define DEFAULT_FOLD 3
// MAXIMUM NUMBER OF K
#define MAX_FOLD 4

#ifdef __SSE__
  #include <xmmintrin.h>
#endif

#ifdef __STDC_IEC_559__
  #include <fenv.h>
#endif

struct env{
  #ifdef __SSE__
    unsigned int csr;
  #endif
  #ifdef __STDC_IEC_559__
    fenv_t fenv;
  #elif defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
    unsigned int controlfp;
  #endif
};

static inline struct env set_env(void){
  (void)set_env;
  struct env env;
  #ifdef __SSE__
    #ifndef _MM_DENORMALS_ZERO_MASK
      #define _MM_DENORMALS_ZERO_MASK 0x0040
    #endif
    #ifndef _MM_DENORMALS_ZERO_ON
      #define _MM_DENORMALS_ZERO_ON   0x0040
    #endif
    #ifndef _MM_DENORMALS_ZERO_OFF
      #define _MM_DENORMALS_ZERO_OFF  0x0040
    #endif
    unsigned int newcsr;
    env.csr = _mm_getcsr();
    newcsr = env.csr;
    newcsr &= ~_MM_DENORMALS_ZERO_MASK;
    newcsr |= _MM_DENORMALS_ZERO_ON;
    newcsr &= ~_MM_FLUSH_ZERO_MASK;
    newcsr |= _MM_FLUSH_ZERO_ON;
    newcsr &= ~_MM_EXCEPT_DENORM;
    newcsr &= ~_MM_EXCEPT_UNDERFLOW;
    newcsr &= ~_MM_ROUND_MASK;
    newcsr |= _MM_ROUND_NEAREST;
    _mm_setcsr(newcsr);
  #endif
  #ifdef __STDC_IEC_559__
    #pragma STDC ACCESS_FENV ON
    fegetenv(&env.env);
    fesetround(FE_TONEAREST);
    feclearexcept(FE_UNDERFLOW);
  #elif defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
    env.controlfp = _controlfp(0, 0xffffffff);
    unsigned int newcontrolfp = env.controlfp;
    newcontrolfp &= ~_MCW_DN;
    newcontrolfp |= _DN_FLUSH;
    newcontrolfp &= ~_EM_DENORMAL;
    newcontrolfp &= ~_EM_UNDERFLOW;
    newcontrolfp &= ~_MCW_RC;
    newcontrolfp |= _RC_NEAR;
    _controlfp(newcontrolfp, 0xffffffff);
  #endif
  return env;
}

static inline void reset_env(struct env env){
  (void)reset_env;
  #ifdef __STDC_IEC_559__
    fesetenv(&env.env);
  #elif defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
    _controlfp(env.controlfp, 0xffffffff);
  #endif
  #ifdef __SSE__
    _mm_setcsr(env.csr);
  #endif
}
