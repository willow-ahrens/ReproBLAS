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
import sumI2
]]]*/
//[[[end]]]

#if defined( __AVX__ )
  void zsumI2(int n, double complex* v, int incv, int fold, double complex* sum){
    /*[[[cog
    cog.out(generate.generate(sumI2.SumI2(dataTypes.FloatComplex, vectorizations.AVX), args, params))
    ]]]*/
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  void zsumI2(int n, double complex* v, int incv, int fold, double complex* sum){
    /*[[[cog
    cog.out(generate.generate(sumI2.SumI2(dataTypes.FloatComplex, vectorizations.SSE), args, params))
    ]]]*/
    //[[[end]]]
  }
#else
  void zsumI2(int n, double complex* v, int incv, int fold, double complex* sum){
    /*[[[cog
    cog.out(generate.generate(sumI2.SumI2(dataTypes.FloatComplex, vectorizations.SISD), args, params))
    ]]]*/
    //[[[end]]]
  }
#endif
