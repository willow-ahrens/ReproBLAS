#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  double damax(int n, double* v, int incv){
    __m256d mask_ABS; AVX_ABS_MASKD(mask_ABS);
    double tmp_max[4] __attribute__((aligned(32)));
    int i;
    double max;

    __m256d v_0, v_1, v_2, v_3;
    __m256d m_0;
    m_0 = _mm256_setzero_pd();

    if(incv == 1){

      for(i = 0; i + 16 <= n; i += 16, v += 16){
        v_0 = _mm256_and_pd(_mm256_loadu_pd(v), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_loadu_pd(v + 4), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_loadu_pd(v + 8), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_loadu_pd(v + 12), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
      }
      if(i + 8 <= n){
        v_0 = _mm256_and_pd(_mm256_loadu_pd(v), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_loadu_pd(v + 4), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 8, v += 8;
      }
      if(i + 4 <= n){
        v_0 = _mm256_and_pd(_mm256_loadu_pd(v), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        i += 4, v += 4;
      }
      if(i < n){
        v_0 = _mm256_and_pd(_mm256_set_pd(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
      }
    }else{

      for(i = 0; i + 16 <= n; i += 16, v += (incv * 16)){
        v_0 = _mm256_and_pd(_mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_set_pd(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_set_pd(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
      }
      if(i + 8 <= n){
        v_0 = _mm256_and_pd(_mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 8, v += (incv * 8);
      }
      if(i + 4 <= n){
        v_0 = _mm256_and_pd(_mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        i += 4, v += (incv * 4);
      }
      if(i < n){
        v_0 = _mm256_and_pd(_mm256_set_pd(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
      }
    }
    _mm256_store_pd(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[1] ? tmp_max[0]: tmp_max[1]);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[0] = (tmp_max[0] > tmp_max[3] ? tmp_max[0]: tmp_max[3]);
    (&max)[0] = ((double*)tmp_max)[0];
    return max;
  }
#elif defined( __SSE2__ )
  double damax(int n, double* v, int incv){
    __m128d mask_ABS; SSE_ABS_MASKD(mask_ABS);
    double tmp_max[2] __attribute__((aligned(16)));
    int i;
    double max;

    __m128d v_0, v_1, v_2, v_3;
    __m128d m_0;
    m_0 = _mm_setzero_pd();

    if(incv == 1){

      for(i = 0; i + 8 <= n; i += 8, v += 8){
        v_0 = _mm_and_pd(_mm_loadu_pd(v), mask_ABS);
        v_1 = _mm_and_pd(_mm_loadu_pd(v + 2), mask_ABS);
        v_2 = _mm_and_pd(_mm_loadu_pd(v + 4), mask_ABS);
        v_3 = _mm_and_pd(_mm_loadu_pd(v + 6), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
      }
      if(i + 4 <= n){
        v_0 = _mm_and_pd(_mm_loadu_pd(v), mask_ABS);
        v_1 = _mm_and_pd(_mm_loadu_pd(v + 2), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 4, v += 4;
      }
      if(i + 2 <= n){
        v_0 = _mm_and_pd(_mm_loadu_pd(v), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        i += 2, v += 2;
      }
      if(i < n){
        v_0 = _mm_and_pd(_mm_set_pd(0, v[0]), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
      }
    }else{

      for(i = 0; i + 8 <= n; i += 8, v += (incv * 8)){
        v_0 = _mm_and_pd(_mm_set_pd(v[incv], v[0]), mask_ABS);
        v_1 = _mm_and_pd(_mm_set_pd(v[(incv * 3)], v[(incv * 2)]), mask_ABS);
        v_2 = _mm_and_pd(_mm_set_pd(v[(incv * 5)], v[(incv * 4)]), mask_ABS);
        v_3 = _mm_and_pd(_mm_set_pd(v[(incv * 7)], v[(incv * 6)]), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
      }
      if(i + 4 <= n){
        v_0 = _mm_and_pd(_mm_set_pd(v[incv], v[0]), mask_ABS);
        v_1 = _mm_and_pd(_mm_set_pd(v[(incv * 3)], v[(incv * 2)]), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 4, v += (incv * 4);
      }
      if(i + 2 <= n){
        v_0 = _mm_and_pd(_mm_set_pd(v[incv], v[0]), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        i += 2, v += (incv * 2);
      }
      if(i < n){
        v_0 = _mm_and_pd(_mm_set_pd(0, v[0]), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
      }
    }
    _mm_store_pd(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[1] ? tmp_max[0]: tmp_max[1]);
    (&max)[0] = ((double*)tmp_max)[0];
    return max;
  }
#else
  double damax(int n, double* v, int incv){
    int i;
    double max;

    double v_0, v_1, v_2, v_3;
    double m_0;
    m_0 = 0;

    if(incv == 1){

      for(i = 0; i + 4 <= n; i += 4, v += 4){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        v_2 = fabs(v[2]);
        v_3 = fabs(v[3]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
      }
      if(i + 2 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        i += 2, v += 2;
      }
      if(i + 1 <= n){
        v_0 = fabs(v[0]);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += 1;
      }
    }else{

      for(i = 0; i + 4 <= n; i += 4, v += (incv * 4)){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        v_2 = fabs(v[(incv * 2)]);
        v_3 = fabs(v[(incv * 3)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
      }
      if(i + 2 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        i += 2, v += (incv * 2);
      }
      if(i + 1 <= n){
        v_0 = fabs(v[0]);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += incv;
      }
    }
    (&max)[0] = m_0;
    return max;
  }
#endif