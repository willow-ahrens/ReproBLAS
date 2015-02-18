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
  float complex camaxm(int n, float complex* v, int incv, float complex* y, int incy){
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    float tmp_max[8] __attribute__((aligned(32)));
    int i;
    float complex max;

    float* v_base = (float*) v;
    float* y_base = (float*) y;
    __m256 v_0, v_1;
    __m256 y_0;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 4 <= n; i += 4, v_base += 8, y_base += 8){
        v_0 = _mm256_loadu_ps(v_base);
        y_0 = _mm256_loadu_ps(y_base);
        v_1 = _mm256_and_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_ABS);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000)), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
      if(i < n){
        v_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[2]:0, (n - i)>1?((double*)v_base)[1]:0, ((double*)v_base)[0]);
        y_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)y_base)[2]:0, (n - i)>1?((double*)y_base)[1]:0, ((double*)y_base)[0]);
        v_1 = _mm256_and_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_ABS);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000)), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
    }else{

      for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8), y_base += (incy * 8)){
        v_0 = _mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
        y_0 = _mm256_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)], y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
        v_1 = _mm256_and_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_ABS);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000)), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
      }
      if(i < n){
        v_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[(incv * 2)]:0, (n - i)>1?((double*)v_base)[incv]:0, ((double*)v_base)[0]);
        y_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)y_base)[(incy * 2)]:0, (n - i)>1?((double*)y_base)[incy]:0, ((double*)y_base)[0]);
        v_1 = _mm256_and_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_ABS);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000)), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
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
  float complex camaxm(int n, float complex* v, int incv, float complex* y, int incy){
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    float tmp_max[4] __attribute__((aligned(16)));
    int i;
    float complex max;

    float* v_base = (float*) v;
    float* y_base = (float*) y;
    __m128 v_0, v_1;
    __m128 y_0;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 2 <= n; i += 2, v_base += 4, y_base += 4){
        v_0 = _mm_loadu_ps(v_base);
        y_0 = _mm_loadu_ps(y_base);
        v_1 = _mm_and_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_ABS);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000)), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
      if(i < n){
        v_0 = _mm_set_ps(0, 0, v_base[1], v_base[0]);
        y_0 = _mm_set_ps(0, 0, y_base[1], y_base[0]);
        v_1 = _mm_and_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_ABS);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000)), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
    }else{

      for(i = 0; i + 2 <= n; i += 2, v_base += (incv * 4), y_base += (incy * 4)){
        v_0 = _mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
        y_0 = _mm_set_ps(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
        v_1 = _mm_and_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_ABS);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000)), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
      if(i < n){
        v_0 = _mm_set_ps(0, 0, v_base[1], v_base[0]);
        y_0 = _mm_set_ps(0, 0, y_base[1], y_base[0]);
        v_1 = _mm_and_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_ABS);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000)), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
      }
    }
    _mm_store_ps(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[1] = (tmp_max[1] > tmp_max[3] ? tmp_max[1]: tmp_max[3]);
    (&max)[0] = ((float complex*)tmp_max)[0];
    return max;
  }
#else
  float complex camaxm(int n, float complex* v, int incv, float complex* y, int incy){
    int i;
    float complex max;

    float* v_base = (float*) v;
    float* y_base = (float*) y;
    float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    float y_0, y_1, y_2, y_3;
    float m_0, m_1;
    m_0 = 0;
    m_1 = 0;

    if(incv == 1 && incy == 1){

      for(i = 0; i + 2 <= n; i += 2, v_base += 4, y_base += 4){
        v_0 = v_base[0];
        v_1 = v_base[1];
        v_2 = v_base[2];
        v_3 = v_base[3];
        y_0 = y_base[0];
        y_1 = y_base[1];
        y_2 = y_base[2];
        y_3 = y_base[3];
        v_4 = fabs(v_1 * y_1);
        v_5 = fabs(v_0 * y_1);
        v_6 = fabs(v_3 * y_3);
        v_7 = fabs(v_2 * y_3);
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_0);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_2);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
      }
      if(i + 1 <= n){
        v_0 = v_base[0];
        v_1 = v_base[1];
        y_0 = y_base[0];
        y_1 = y_base[1];
        v_2 = fabs(v_1 * y_1);
        v_3 = fabs(v_0 * y_1);
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_0);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        i += 1, v_base += 2, y_base += 2;
      }
    }else{

      for(i = 0; i + 2 <= n; i += 2, v_base += (incv * 4), y_base += (incy * 4)){
        v_0 = v_base[0];
        v_1 = v_base[1];
        v_2 = v_base[(incv * 2)];
        v_3 = v_base[((incv * 2) + 1)];
        y_0 = y_base[0];
        y_1 = y_base[1];
        y_2 = y_base[(incy * 2)];
        y_3 = y_base[((incy * 2) + 1)];
        v_4 = fabs(v_1 * y_1);
        v_5 = fabs(v_0 * y_1);
        v_6 = fabs(v_3 * y_3);
        v_7 = fabs(v_2 * y_3);
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_0);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_2);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
      }
      if(i + 1 <= n){
        v_0 = v_base[0];
        v_1 = v_base[1];
        y_0 = y_base[0];
        y_1 = y_base[1];
        v_2 = fabs(v_1 * y_1);
        v_3 = fabs(v_0 * y_1);
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_0);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        i += 1, v_base += (incv * 2), y_base += (incy * 2);
      }
    }
    ((float*)(&max))[0] = m_0;
    ((float*)(&max))[1] = m_1;
    return max;
  }
#endif