#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../config.h"
#include "../Common/Common.h"
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
  float samaxm(int n, float* v, int incv, float* y, int incy){
    /*[[[cog
    cog.out(generate.generate(amax.AMaxM(dataTypes.Float, vectorizations.AVX), args, params))
    ]]]*/
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    float tmp_max[8] __attribute__((aligned(32)));
    int i;
    float max;

    __m256 v_0, v_1;
    __m256 y_0, y_1;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 16 <= n; i += 16, v += 16, y += 16){
        v_0 = _mm256_loadu_ps(v);
        v_1 = _mm256_loadu_ps(v + 8);
        y_0 = _mm256_loadu_ps(y);
        y_1 = _mm256_loadu_ps(y + 8);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
      if(i + 8 <= n){
        v_0 = _mm256_loadu_ps(v);
        y_0 = _mm256_loadu_ps(y);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        i += 8, v += 8, y += 8;
      }
      if(i < n){
        v_0 = _mm256_set_ps(0, (n - i)>6?v[6]:0, (n - i)>5?v[5]:0, (n - i)>4?v[4]:0, (n - i)>3?v[3]:0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
        y_0 = _mm256_set_ps(0, (n - i)>6?y[6]:0, (n - i)>5?y[5]:0, (n - i)>4?y[4]:0, (n - i)>3?y[3]:0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
      }
    }else{

      for(i = 0; i + 16 <= n; i += 16, v += (incv * 16), y += (incy * 16)){
        v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm256_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)], y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
      if(i + 8 <= n){
        v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        i += 8, v += (incv * 8), y += (incy * 8);
      }
      if(i < n){
        v_0 = _mm256_set_ps(0, (n - i)>6?v[(incv * 6)]:0, (n - i)>5?v[(incv * 5)]:0, (n - i)>4?v[(incv * 4)]:0, (n - i)>3?v[(incv * 3)]:0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
        y_0 = _mm256_set_ps(0, (n - i)>6?y[(incy * 6)]:0, (n - i)>5?y[(incy * 5)]:0, (n - i)>4?y[(incy * 4)]:0, (n - i)>3?y[(incy * 3)]:0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
      }
    }
    _mm256_store_ps(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[1] ? tmp_max[0]: tmp_max[1]);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[0] = (tmp_max[0] > tmp_max[3] ? tmp_max[0]: tmp_max[3]);
    tmp_max[0] = (tmp_max[0] > tmp_max[4] ? tmp_max[0]: tmp_max[4]);
    tmp_max[0] = (tmp_max[0] > tmp_max[5] ? tmp_max[0]: tmp_max[5]);
    tmp_max[0] = (tmp_max[0] > tmp_max[6] ? tmp_max[0]: tmp_max[6]);
    tmp_max[0] = (tmp_max[0] > tmp_max[7] ? tmp_max[0]: tmp_max[7]);
    (&max)[0] = ((float*)tmp_max)[0];
    return max;
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  float samaxm(int n, float* v, int incv, float* y, int incy){
    /*[[[cog
    cog.out(generate.generate(amax.AMaxM(dataTypes.Float, vectorizations.SSE), args, params))
    ]]]*/
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    float tmp_max[4] __attribute__((aligned(16)));
    int i;
    float max;

    __m128 v_0, v_1;
    __m128 y_0, y_1;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
        v_0 = _mm_loadu_ps(v);
        v_1 = _mm_loadu_ps(v + 4);
        y_0 = _mm_loadu_ps(y);
        y_1 = _mm_loadu_ps(y + 4);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
      if(i + 4 <= n){
        v_0 = _mm_loadu_ps(v);
        y_0 = _mm_loadu_ps(y);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        i += 4, v += 4, y += 4;
      }
      if(i < n){
        v_0 = _mm_set_ps(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
        y_0 = _mm_set_ps(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
      }
    }else{

      for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
        v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
      if(i + 4 <= n){
        v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        i += 4, v += (incv * 4), y += (incy * 4);
      }
      if(i < n){
        v_0 = _mm_set_ps(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
        y_0 = _mm_set_ps(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
      }
    }
    _mm_store_ps(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[1] ? tmp_max[0]: tmp_max[1]);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[0] = (tmp_max[0] > tmp_max[3] ? tmp_max[0]: tmp_max[3]);
    (&max)[0] = ((float*)tmp_max)[0];
    return max;
    //[[[end]]]
  }
#else
  float samaxm(int n, float* v, int incv, float* y, int incy){
    /*[[[cog
    cog.out(generate.generate(amax.AMaxM(dataTypes.Float, vectorizations.SISD), args, params))
    ]]]*/
    int i;
    float max;

    float v_0, v_1;
    float y_0, y_1;
    float m_0;
    m_0 = 0;

    if(incv == 1 && incy == 1){

      for(i = 0; i + 2 <= n; i += 2, v += 2, y += 2){
        v_0 = v[0];
        v_1 = v[1];
        y_0 = y[0];
        y_1 = y[1];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
      }
      if(i + 1 <= n){
        v_0 = v[0];
        y_0 = y[0];
        v_0 = fabs(v_0 * y_0);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += 1, y += 1;
      }
    }else{

      for(i = 0; i + 2 <= n; i += 2, v += (incv * 2), y += (incy * 2)){
        v_0 = v[0];
        v_1 = v[incv];
        y_0 = y[0];
        y_1 = y[incy];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
      }
      if(i + 1 <= n){
        v_0 = v[0];
        y_0 = y[0];
        v_0 = fabs(v_0 * y_0);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += incv, y += incy;
      }
    }
    (&max)[0] = m_0;
    return max;
    //[[[end]]]
  }
#endif
