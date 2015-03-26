#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "../config.h"
#include "../Common/Common.h"
#include "indexedBLAS.h"
#include <immintrin.h>
#include <emmintrin.h>

/*[[[cog
import cog
import sys, os
import generate
import dataTypes
import vectorizations
import amax
]]]*/
//[[[end]]]

#if defined( __AVX__ )
  float complex camax(int n, float complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.FloatComplex, vectorizations.AVX), args, params))
    ]]]*/
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    float tmp_max[8] __attribute__((aligned(32)));
    int i;
    float complex max;

    float* v_base = (float*) v;
    __m256 v_0, v_1;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 8 <= n; i += 8, v_base += 16){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
      if(i + 4 <= n){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        i += 4, v_base += 8;
      }
      if(i < n){
        v_0 = _mm256_and_ps((__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[2]:0, (n - i)>1?((double*)v_base)[1]:0, ((double*)v_base)[0]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
      }
    }else{

      for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
        v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
      if(i + 4 <= n){
        v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        i += 4, v_base += (incv * 8);
      }
      if(i < n){
        v_0 = _mm256_and_ps((__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[(incv * 2)]:0, (n - i)>1?((double*)v_base)[incv]:0, ((double*)v_base)[0]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
      }
    }
    _mm256_store_ps(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[1] = (tmp_max[1] > tmp_max[3] ? tmp_max[1]: tmp_max[3]);
    tmp_max[0] = (tmp_max[0] > tmp_max[4] ? tmp_max[0]: tmp_max[4]);
    tmp_max[1] = (tmp_max[1] > tmp_max[5] ? tmp_max[1]: tmp_max[5]);
    tmp_max[0] = (tmp_max[0] > tmp_max[6] ? tmp_max[0]: tmp_max[6]);
    tmp_max[1] = (tmp_max[1] > tmp_max[7] ? tmp_max[1]: tmp_max[7]);
    (&max)[0] = ((float complex*)tmp_max)[0];
    return max;
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  float complex camax(int n, float complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.FloatComplex, vectorizations.SSE), args, params))
    ]]]*/
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    float tmp_max[4] __attribute__((aligned(16)));
    int i;
    float complex max;

    float* v_base = (float*) v;
    __m128 v_0, v_1;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 4 <= n; i += 4, v_base += 8){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
      if(i + 2 <= n){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        i += 2, v_base += 4;
      }
      if(i < n){
        v_0 = _mm_and_ps(_mm_set_ps(0, 0, v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
      }
    }else{

      for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
      if(i + 2 <= n){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        i += 2, v_base += (incv * 4);
      }
      if(i < n){
        v_0 = _mm_and_ps(_mm_set_ps(0, 0, v_base[1], v_base[0]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
      }
    }
    _mm_store_ps(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[1] = (tmp_max[1] > tmp_max[3] ? tmp_max[1]: tmp_max[3]);
    (&max)[0] = ((float complex*)tmp_max)[0];
    return max;
    //[[[end]]]
  }
#else
  float complex camax(int n, float complex* v, int incv){
    /*[[[cog
    cog.out(generate.generate(amax.AMax(dataTypes.FloatComplex, vectorizations.SISD), args, params))
    ]]]*/
    int i;
    float complex max;

    float* v_base = (float*) v;
    float v_0, v_1, v_2, v_3;
    float m_0, m_1;
    m_0 = 0;
    m_1 = 0;

    if(incv == 1){

      for(i = 0; i + 2 <= n; i += 2, v_base += 4){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[2]);
        v_3 = fabs(v_base[3]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
      }
      if(i + 1 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        i += 1, v_base += 2;
      }
    }else{

      for(i = 0; i + 2 <= n; i += 2, v_base += (incv * 4)){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[(incv * 2)]);
        v_3 = fabs(v_base[((incv * 2) + 1)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
      }
      if(i + 1 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        i += 1, v_base += (incv * 2);
      }
    }
    ((float*)(&max))[0] = m_0;
    ((float*)(&max))[1] = m_1;
    return max;
    //[[[end]]]
  }
#endif
