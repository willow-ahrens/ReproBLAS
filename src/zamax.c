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
  double complex zamax(int n, double complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.DoubleComplex, vectorizations.AVX), args, params))
    ]]]*/
    __m256d mask_ABS; AVX_ABS_MASKD(mask_ABS);
    double tmp_max[4] __attribute__((aligned(32)));
    int i;
    double complex max;

    double* v_base = (double*) v;
    __m256d v_0, v_1, v_2, v_3;
    __m256d m_0;
    m_0 = _mm256_setzero_pd();

    if(incv == 1){

      for(i = 0; i + 8 <= n; i += 8, v_base += 16){
        v_0 = _mm256_and_pd(_mm256_loadu_pd(v_base), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_loadu_pd(v_base + 4), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_loadu_pd(v_base + 8), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_loadu_pd(v_base + 12), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
      }
      if(i + 4 <= n){
        v_0 = _mm256_and_pd(_mm256_loadu_pd(v_base), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_loadu_pd(v_base + 4), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 4, v_base += 8;
      }
      if(i + 2 <= n){
        v_0 = _mm256_and_pd(_mm256_loadu_pd(v_base), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        i += 2, v_base += 4;
      }
      if(i < n){
        v_0 = _mm256_and_pd(_mm256_set_pd(0, 0, v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
      }
    }else{

      for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
        v_0 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
      }
      if(i + 4 <= n){
        v_0 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 4, v_base += (incv * 8);
      }
      if(i + 2 <= n){
        v_0 = _mm256_and_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        i += 2, v_base += (incv * 4);
      }
      if(i < n){
        v_0 = _mm256_and_pd(_mm256_set_pd(0, 0, v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
      }
    }
    _mm256_store_pd(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[1] = (tmp_max[1] > tmp_max[3] ? tmp_max[1]: tmp_max[3]);
    (&max)[0] = ((double complex*)tmp_max)[0];
    return max;
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  double complex zamax(int n, double complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.DoubleComplex, vectorizations.SSE), args, params))
    ]]]*/
    __m128d mask_ABS; SSE_ABS_MASKD(mask_ABS);
    double tmp_max[2] __attribute__((aligned(16)));
    int i;
    double complex max;

    double* v_base = (double*) v;
    __m128d v_0, v_1, v_2, v_3;
    __m128d m_0;
    m_0 = _mm_setzero_pd();

    if(incv == 1){

      for(i = 0; i + 4 <= n; i += 4, v_base += 8){
        v_0 = _mm_and_pd(_mm_loadu_pd(v_base), mask_ABS);
        v_1 = _mm_and_pd(_mm_loadu_pd(v_base + 2), mask_ABS);
        v_2 = _mm_and_pd(_mm_loadu_pd(v_base + 4), mask_ABS);
        v_3 = _mm_and_pd(_mm_loadu_pd(v_base + 6), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
      }
      if(i + 2 <= n){
        v_0 = _mm_and_pd(_mm_loadu_pd(v_base), mask_ABS);
        v_1 = _mm_and_pd(_mm_loadu_pd(v_base + 2), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 2, v_base += 4;
      }
      if(i + 1 <= n){
        v_0 = _mm_and_pd(_mm_loadu_pd(v_base), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        i += 1, v_base += 2;
      }
    }else{

      for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
        v_0 = _mm_and_pd(_mm_loadu_pd(v_base), mask_ABS);
        v_1 = _mm_and_pd(_mm_loadu_pd(v_base + (incv * 2)), mask_ABS);
        v_2 = _mm_and_pd(_mm_loadu_pd(v_base + (incv * 4)), mask_ABS);
        v_3 = _mm_and_pd(_mm_loadu_pd(v_base + (incv * 6)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
      }
      if(i + 2 <= n){
        v_0 = _mm_and_pd(_mm_loadu_pd(v_base), mask_ABS);
        v_1 = _mm_and_pd(_mm_loadu_pd(v_base + (incv * 2)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 2, v_base += (incv * 4);
      }
      if(i + 1 <= n){
        v_0 = _mm_and_pd(_mm_loadu_pd(v_base), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        i += 1, v_base += (incv * 2);
      }
    }
    _mm_store_pd(tmp_max, m_0);
    (&max)[0] = ((double complex*)tmp_max)[0];
    return max;
    //[[[end]]]
  }
#else
  double complex zamax(int n, double complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.DoubleComplex, vectorizations.SISD), args, params))
    ]]]*/
    int i;
    double complex max;

    double* v_base = (double*) v;
    double v_0, v_1;
    double m_0, m_1;
    m_0 = 0;
    m_1 = 0;

    if(incv == 1){

      for(i = 0; i + 1 <= n; i += 1, v_base += 2){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
      }
    }else{

      for(i = 0; i + 1 <= n; i += 1, v_base += (incv * 2)){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
      }
    }
    ((double*)(&max))[0] = m_0;
    ((double*)(&max))[1] = m_1;
    return max;
    //[[[end]]]
  }
#endif
