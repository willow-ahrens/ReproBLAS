/**
 * @file  config.h
 * @brief config.h is the configuration file for c code in the library.
 *
 * Feel free to modify anything here, but be sure to read the comments describing all of the fields
 */

/**
 * @brief The fold of indexed types used anywhere fold is not specified.
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define DEFAULT_FOLD 3

/**
 * @brief The maximum fold supported by the library.
 *
 * Increases to #MAX_FOLD will allow users to use larger folds (meaning greater accuracy), but will also increase suntime (by a small amount) and the time it takes to autotune the library (by a large amount).
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define MAX_FOLD 4

#ifdef __SSE__
  #include <xmmintrin.h>
#endif

#ifdef __STDC_IEC_559__
  #include <fenv.h>
#endif

/**
 * @brief Used to store the floating point environment before it is modified by ReproBLAS.
 *
 * If you are modifying #set_env or #get_env, feel free to add any appropriate fields
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
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

/**
 * @brief Set up the floating point environment for ReproBLAS
 *
 * #set_env needs to set both the SIMD (vector) and SISD (single) floating point unit rounding to "to nearest," and disable underflow exceptions for both units. Optionally, #set_env can set the "denormals are zero" and "flush to zero" modes, which will increase performance without affecting the result.
 *
 * @return the floating point environment before modification
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
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

/**
 * @brief Restore the floating point environment after ReproBLAS uses it.
 *
 * #reset_env has the easier (and actually optional) job of restoring the floating point environments after calls to #set_env.
 *
 * @return the floating point environment before modification
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
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
