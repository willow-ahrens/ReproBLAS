#include <float.h>
#include <stdio.h>
#include <stdlib.h>
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
  void ssumI2(int n, float* v, int incv, int fold, float* sum){
    /*[[[cog
    cog.out(generate.generate(sumI2.SumI2(dataTypes.Float, vectorizations.AVX), args, params))
    ]]]*/
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  void ssumI2(int n, float* v, int incv, int fold, float* sum){
    /*[[[cog
    cog.out(generate.generate(sumI2.SumI2(dataTypes.Float, vectorizations.SSE), args, params))
    ]]]*/
    //[[[end]]]
  }
#else
  void ssumI2(int n, float* v, int incv, int fold, float* sum){
    /*[[[cog
    cog.out(generate.generate(sumI2.SumI2(dataTypes.Float, vectorizations.SISD), args, params))
    ]]]*/
    //[[[end]]]
  }
#endif
