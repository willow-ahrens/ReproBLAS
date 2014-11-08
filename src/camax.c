#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  float complex camax(int n, float complex* v, int incv){
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    float tmp_max[8] __attribute__((aligned(32)));
    int i;
    float complex max;

    float* v_base = (float*) v;
    __m256 v_0, v_1, v_2, v_3;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 16 <= n; i += 16, v_base += 32){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_loadu_ps(v_base + 16), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_loadu_ps(v_base + 24), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
      }
      if(i + 8 <= n){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        i += 8, v_base += 16;
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

      for(i = 0; i + 16 <= n; i += 16, v_base += (incv * 32)){
        v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)], v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)], v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
      }
      if(i + 8 <= n){
        v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        i += 8, v_base += (incv * 16);
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
  }
#elif defined( __SSE2__ )
  float complex camax(int n, float complex* v, int incv){
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    float tmp_max[4] __attribute__((aligned(16)));
    int i;
    float complex max;

    float* v_base = (float*) v;
    __m128 v_0, v_1, v_2, v_3;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 8 <= n; i += 8, v_base += 16){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
        v_2 = _mm_and_ps(_mm_loadu_ps(v_base + 8), mask_ABS);
        v_3 = _mm_and_ps(_mm_loadu_ps(v_base + 12), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
      }
      if(i + 4 <= n){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        i += 4, v_base += 8;
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

      for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        v_2 = _mm_and_ps(_mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_3 = _mm_and_ps(_mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
      }
      if(i + 4 <= n){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        i += 4, v_base += (incv * 8);
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
  }
#else
  float complex camax(int n, float complex* v, int incv){
    int i;
    float complex max;

    float* v_base = (float*) v;
    float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    float m_0, m_1;
    m_0 = 0;
    m_1 = 0;

    if(incv == 1){

      for(i = 0; i + 4 <= n; i += 4, v_base += 8){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[2]);
        v_3 = fabs(v_base[3]);
        v_4 = fabs(v_base[4]);
        v_5 = fabs(v_base[5]);
        v_6 = fabs(v_base[6]);
        v_7 = fabs(v_base[7]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
      }
      if(i + 2 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[2]);
        v_3 = fabs(v_base[3]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        i += 2, v_base += 4;
      }
      if(i + 1 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        i += 1, v_base += 2;
      }
    }else{

      for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[(incv * 2)]);
        v_3 = fabs(v_base[((incv * 2) + 1)]);
        v_4 = fabs(v_base[(incv * 4)]);
        v_5 = fabs(v_base[((incv * 4) + 1)]);
        v_6 = fabs(v_base[(incv * 6)]);
        v_7 = fabs(v_base[((incv * 6) + 1)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
      }
      if(i + 2 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[(incv * 2)]);
        v_3 = fabs(v_base[((incv * 2) + 1)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        i += 2, v_base += (incv * 4);
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
  }
#endif