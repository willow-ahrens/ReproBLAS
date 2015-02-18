#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>

/*[[[cog
import cog
import sys, os
from gen import generate
from gen import dataTypes
from gen import vectorizations
import amax
]]]*/
//[[[end]]]

#if defined( __AVX__ )
  float complex camax(int n, float complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.FloatComplex, vectorizations.AVX), args, params))
    ]]]*/
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  float complex camax(int n, float complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.FloatComplex, vectorizations.SSE), args, params))
    ]]]*/
    //[[[end]]]
  }
#else
  float complex camax(int n, float complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.FloatComplex, vectorizations.SIMD), args, params))
    ]]]*/
    //[[[end]]]
  }
#endif
