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
    __m256 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 32 <= n; i += 32, v_base += 64){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_loadu_ps(v_base + 16), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_loadu_ps(v_base + 24), mask_ABS);
        v_4 = _mm256_and_ps(_mm256_loadu_ps(v_base + 32), mask_ABS);
        v_5 = _mm256_and_ps(_mm256_loadu_ps(v_base + 40), mask_ABS);
        v_6 = _mm256_and_ps(_mm256_loadu_ps(v_base + 48), mask_ABS);
        v_7 = _mm256_and_ps(_mm256_loadu_ps(v_base + 56), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
        m_0 = _mm256_max_ps(m_0, v_4);
        m_0 = _mm256_max_ps(m_0, v_5);
        m_0 = _mm256_max_ps(m_0, v_6);
        m_0 = _mm256_max_ps(m_0, v_7);
      }
      if(i + 16 <= n){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_loadu_ps(v_base + 16), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_loadu_ps(v_base + 24), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
        i += 16, v_base += 32;
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

      for(i = 0; i + 32 <= n; i += 32, v_base += (incv * 64)){
        v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)], v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)], v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]), mask_ABS);
        v_4 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 38) + 1)], v_base[(incv * 38)], v_base[((incv * 36) + 1)], v_base[(incv * 36)], v_base[((incv * 34) + 1)], v_base[(incv * 34)], v_base[((incv * 32) + 1)], v_base[(incv * 32)]), mask_ABS);
        v_5 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 46) + 1)], v_base[(incv * 46)], v_base[((incv * 44) + 1)], v_base[(incv * 44)], v_base[((incv * 42) + 1)], v_base[(incv * 42)], v_base[((incv * 40) + 1)], v_base[(incv * 40)]), mask_ABS);
        v_6 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 54) + 1)], v_base[(incv * 54)], v_base[((incv * 52) + 1)], v_base[(incv * 52)], v_base[((incv * 50) + 1)], v_base[(incv * 50)], v_base[((incv * 48) + 1)], v_base[(incv * 48)]), mask_ABS);
        v_7 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 62) + 1)], v_base[(incv * 62)], v_base[((incv * 60) + 1)], v_base[(incv * 60)], v_base[((incv * 58) + 1)], v_base[(incv * 58)], v_base[((incv * 56) + 1)], v_base[(incv * 56)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
        m_0 = _mm256_max_ps(m_0, v_4);
        m_0 = _mm256_max_ps(m_0, v_5);
        m_0 = _mm256_max_ps(m_0, v_6);
        m_0 = _mm256_max_ps(m_0, v_7);
      }
      if(i + 16 <= n){
        v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)], v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)], v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
        i += 16, v_base += (incv * 32);
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
    __m128 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11, v_12, v_13, v_14, v_15;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 32 <= n; i += 32, v_base += 64){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
        v_2 = _mm_and_ps(_mm_loadu_ps(v_base + 8), mask_ABS);
        v_3 = _mm_and_ps(_mm_loadu_ps(v_base + 12), mask_ABS);
        v_4 = _mm_and_ps(_mm_loadu_ps(v_base + 16), mask_ABS);
        v_5 = _mm_and_ps(_mm_loadu_ps(v_base + 20), mask_ABS);
        v_6 = _mm_and_ps(_mm_loadu_ps(v_base + 24), mask_ABS);
        v_7 = _mm_and_ps(_mm_loadu_ps(v_base + 28), mask_ABS);
        v_8 = _mm_and_ps(_mm_loadu_ps(v_base + 32), mask_ABS);
        v_9 = _mm_and_ps(_mm_loadu_ps(v_base + 36), mask_ABS);
        v_10 = _mm_and_ps(_mm_loadu_ps(v_base + 40), mask_ABS);
        v_11 = _mm_and_ps(_mm_loadu_ps(v_base + 44), mask_ABS);
        v_12 = _mm_and_ps(_mm_loadu_ps(v_base + 48), mask_ABS);
        v_13 = _mm_and_ps(_mm_loadu_ps(v_base + 52), mask_ABS);
        v_14 = _mm_and_ps(_mm_loadu_ps(v_base + 56), mask_ABS);
        v_15 = _mm_and_ps(_mm_loadu_ps(v_base + 60), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        m_0 = _mm_max_ps(m_0, v_4);
        m_0 = _mm_max_ps(m_0, v_5);
        m_0 = _mm_max_ps(m_0, v_6);
        m_0 = _mm_max_ps(m_0, v_7);
        m_0 = _mm_max_ps(m_0, v_8);
        m_0 = _mm_max_ps(m_0, v_9);
        m_0 = _mm_max_ps(m_0, v_10);
        m_0 = _mm_max_ps(m_0, v_11);
        m_0 = _mm_max_ps(m_0, v_12);
        m_0 = _mm_max_ps(m_0, v_13);
        m_0 = _mm_max_ps(m_0, v_14);
        m_0 = _mm_max_ps(m_0, v_15);
      }
      if(i + 16 <= n){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
        v_2 = _mm_and_ps(_mm_loadu_ps(v_base + 8), mask_ABS);
        v_3 = _mm_and_ps(_mm_loadu_ps(v_base + 12), mask_ABS);
        v_4 = _mm_and_ps(_mm_loadu_ps(v_base + 16), mask_ABS);
        v_5 = _mm_and_ps(_mm_loadu_ps(v_base + 20), mask_ABS);
        v_6 = _mm_and_ps(_mm_loadu_ps(v_base + 24), mask_ABS);
        v_7 = _mm_and_ps(_mm_loadu_ps(v_base + 28), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        m_0 = _mm_max_ps(m_0, v_4);
        m_0 = _mm_max_ps(m_0, v_5);
        m_0 = _mm_max_ps(m_0, v_6);
        m_0 = _mm_max_ps(m_0, v_7);
        i += 16, v_base += 32;
      }
      if(i + 8 <= n){
        v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
        v_2 = _mm_and_ps(_mm_loadu_ps(v_base + 8), mask_ABS);
        v_3 = _mm_and_ps(_mm_loadu_ps(v_base + 12), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        i += 8, v_base += 16;
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

      for(i = 0; i + 32 <= n; i += 32, v_base += (incv * 64)){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        v_2 = _mm_and_ps(_mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_3 = _mm_and_ps(_mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
        v_4 = _mm_and_ps(_mm_set_ps(v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]), mask_ABS);
        v_5 = _mm_and_ps(_mm_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)]), mask_ABS);
        v_6 = _mm_and_ps(_mm_set_ps(v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]), mask_ABS);
        v_7 = _mm_and_ps(_mm_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)]), mask_ABS);
        v_8 = _mm_and_ps(_mm_set_ps(v_base[((incv * 34) + 1)], v_base[(incv * 34)], v_base[((incv * 32) + 1)], v_base[(incv * 32)]), mask_ABS);
        v_9 = _mm_and_ps(_mm_set_ps(v_base[((incv * 38) + 1)], v_base[(incv * 38)], v_base[((incv * 36) + 1)], v_base[(incv * 36)]), mask_ABS);
        v_10 = _mm_and_ps(_mm_set_ps(v_base[((incv * 42) + 1)], v_base[(incv * 42)], v_base[((incv * 40) + 1)], v_base[(incv * 40)]), mask_ABS);
        v_11 = _mm_and_ps(_mm_set_ps(v_base[((incv * 46) + 1)], v_base[(incv * 46)], v_base[((incv * 44) + 1)], v_base[(incv * 44)]), mask_ABS);
        v_12 = _mm_and_ps(_mm_set_ps(v_base[((incv * 50) + 1)], v_base[(incv * 50)], v_base[((incv * 48) + 1)], v_base[(incv * 48)]), mask_ABS);
        v_13 = _mm_and_ps(_mm_set_ps(v_base[((incv * 54) + 1)], v_base[(incv * 54)], v_base[((incv * 52) + 1)], v_base[(incv * 52)]), mask_ABS);
        v_14 = _mm_and_ps(_mm_set_ps(v_base[((incv * 58) + 1)], v_base[(incv * 58)], v_base[((incv * 56) + 1)], v_base[(incv * 56)]), mask_ABS);
        v_15 = _mm_and_ps(_mm_set_ps(v_base[((incv * 62) + 1)], v_base[(incv * 62)], v_base[((incv * 60) + 1)], v_base[(incv * 60)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        m_0 = _mm_max_ps(m_0, v_4);
        m_0 = _mm_max_ps(m_0, v_5);
        m_0 = _mm_max_ps(m_0, v_6);
        m_0 = _mm_max_ps(m_0, v_7);
        m_0 = _mm_max_ps(m_0, v_8);
        m_0 = _mm_max_ps(m_0, v_9);
        m_0 = _mm_max_ps(m_0, v_10);
        m_0 = _mm_max_ps(m_0, v_11);
        m_0 = _mm_max_ps(m_0, v_12);
        m_0 = _mm_max_ps(m_0, v_13);
        m_0 = _mm_max_ps(m_0, v_14);
        m_0 = _mm_max_ps(m_0, v_15);
      }
      if(i + 16 <= n){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        v_2 = _mm_and_ps(_mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_3 = _mm_and_ps(_mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
        v_4 = _mm_and_ps(_mm_set_ps(v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]), mask_ABS);
        v_5 = _mm_and_ps(_mm_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)]), mask_ABS);
        v_6 = _mm_and_ps(_mm_set_ps(v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]), mask_ABS);
        v_7 = _mm_and_ps(_mm_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        m_0 = _mm_max_ps(m_0, v_4);
        m_0 = _mm_max_ps(m_0, v_5);
        m_0 = _mm_max_ps(m_0, v_6);
        m_0 = _mm_max_ps(m_0, v_7);
        i += 16, v_base += (incv * 32);
      }
      if(i + 8 <= n){
        v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
        v_2 = _mm_and_ps(_mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
        v_3 = _mm_and_ps(_mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        i += 8, v_base += (incv * 16);
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
    float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11, v_12, v_13, v_14, v_15, v_16, v_17, v_18, v_19, v_20, v_21, v_22, v_23, v_24, v_25, v_26, v_27, v_28, v_29, v_30, v_31, v_32, v_33, v_34, v_35, v_36, v_37, v_38, v_39, v_40, v_41, v_42, v_43, v_44, v_45, v_46, v_47, v_48, v_49, v_50, v_51, v_52, v_53, v_54, v_55, v_56, v_57, v_58, v_59, v_60, v_61, v_62, v_63;
    float m_0, m_1;
    m_0 = 0;
    m_1 = 0;

    if(incv == 1){

      for(i = 0; i + 32 <= n; i += 32, v_base += 64){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[2]);
        v_3 = fabs(v_base[3]);
        v_4 = fabs(v_base[4]);
        v_5 = fabs(v_base[5]);
        v_6 = fabs(v_base[6]);
        v_7 = fabs(v_base[7]);
        v_8 = fabs(v_base[8]);
        v_9 = fabs(v_base[9]);
        v_10 = fabs(v_base[10]);
        v_11 = fabs(v_base[11]);
        v_12 = fabs(v_base[12]);
        v_13 = fabs(v_base[13]);
        v_14 = fabs(v_base[14]);
        v_15 = fabs(v_base[15]);
        v_16 = fabs(v_base[16]);
        v_17 = fabs(v_base[17]);
        v_18 = fabs(v_base[18]);
        v_19 = fabs(v_base[19]);
        v_20 = fabs(v_base[20]);
        v_21 = fabs(v_base[21]);
        v_22 = fabs(v_base[22]);
        v_23 = fabs(v_base[23]);
        v_24 = fabs(v_base[24]);
        v_25 = fabs(v_base[25]);
        v_26 = fabs(v_base[26]);
        v_27 = fabs(v_base[27]);
        v_28 = fabs(v_base[28]);
        v_29 = fabs(v_base[29]);
        v_30 = fabs(v_base[30]);
        v_31 = fabs(v_base[31]);
        v_32 = fabs(v_base[32]);
        v_33 = fabs(v_base[33]);
        v_34 = fabs(v_base[34]);
        v_35 = fabs(v_base[35]);
        v_36 = fabs(v_base[36]);
        v_37 = fabs(v_base[37]);
        v_38 = fabs(v_base[38]);
        v_39 = fabs(v_base[39]);
        v_40 = fabs(v_base[40]);
        v_41 = fabs(v_base[41]);
        v_42 = fabs(v_base[42]);
        v_43 = fabs(v_base[43]);
        v_44 = fabs(v_base[44]);
        v_45 = fabs(v_base[45]);
        v_46 = fabs(v_base[46]);
        v_47 = fabs(v_base[47]);
        v_48 = fabs(v_base[48]);
        v_49 = fabs(v_base[49]);
        v_50 = fabs(v_base[50]);
        v_51 = fabs(v_base[51]);
        v_52 = fabs(v_base[52]);
        v_53 = fabs(v_base[53]);
        v_54 = fabs(v_base[54]);
        v_55 = fabs(v_base[55]);
        v_56 = fabs(v_base[56]);
        v_57 = fabs(v_base[57]);
        v_58 = fabs(v_base[58]);
        v_59 = fabs(v_base[59]);
        v_60 = fabs(v_base[60]);
        v_61 = fabs(v_base[61]);
        v_62 = fabs(v_base[62]);
        v_63 = fabs(v_base[63]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_1 = (m_1 > v_9? m_1: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_1 = (m_1 > v_11? m_1: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_1 = (m_1 > v_13? m_1: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_1 = (m_1 > v_15? m_1: v_15);
        m_0 = (m_0 > v_16? m_0: v_16);
        m_1 = (m_1 > v_17? m_1: v_17);
        m_0 = (m_0 > v_18? m_0: v_18);
        m_1 = (m_1 > v_19? m_1: v_19);
        m_0 = (m_0 > v_20? m_0: v_20);
        m_1 = (m_1 > v_21? m_1: v_21);
        m_0 = (m_0 > v_22? m_0: v_22);
        m_1 = (m_1 > v_23? m_1: v_23);
        m_0 = (m_0 > v_24? m_0: v_24);
        m_1 = (m_1 > v_25? m_1: v_25);
        m_0 = (m_0 > v_26? m_0: v_26);
        m_1 = (m_1 > v_27? m_1: v_27);
        m_0 = (m_0 > v_28? m_0: v_28);
        m_1 = (m_1 > v_29? m_1: v_29);
        m_0 = (m_0 > v_30? m_0: v_30);
        m_1 = (m_1 > v_31? m_1: v_31);
        m_0 = (m_0 > v_32? m_0: v_32);
        m_1 = (m_1 > v_33? m_1: v_33);
        m_0 = (m_0 > v_34? m_0: v_34);
        m_1 = (m_1 > v_35? m_1: v_35);
        m_0 = (m_0 > v_36? m_0: v_36);
        m_1 = (m_1 > v_37? m_1: v_37);
        m_0 = (m_0 > v_38? m_0: v_38);
        m_1 = (m_1 > v_39? m_1: v_39);
        m_0 = (m_0 > v_40? m_0: v_40);
        m_1 = (m_1 > v_41? m_1: v_41);
        m_0 = (m_0 > v_42? m_0: v_42);
        m_1 = (m_1 > v_43? m_1: v_43);
        m_0 = (m_0 > v_44? m_0: v_44);
        m_1 = (m_1 > v_45? m_1: v_45);
        m_0 = (m_0 > v_46? m_0: v_46);
        m_1 = (m_1 > v_47? m_1: v_47);
        m_0 = (m_0 > v_48? m_0: v_48);
        m_1 = (m_1 > v_49? m_1: v_49);
        m_0 = (m_0 > v_50? m_0: v_50);
        m_1 = (m_1 > v_51? m_1: v_51);
        m_0 = (m_0 > v_52? m_0: v_52);
        m_1 = (m_1 > v_53? m_1: v_53);
        m_0 = (m_0 > v_54? m_0: v_54);
        m_1 = (m_1 > v_55? m_1: v_55);
        m_0 = (m_0 > v_56? m_0: v_56);
        m_1 = (m_1 > v_57? m_1: v_57);
        m_0 = (m_0 > v_58? m_0: v_58);
        m_1 = (m_1 > v_59? m_1: v_59);
        m_0 = (m_0 > v_60? m_0: v_60);
        m_1 = (m_1 > v_61? m_1: v_61);
        m_0 = (m_0 > v_62? m_0: v_62);
        m_1 = (m_1 > v_63? m_1: v_63);
      }
      if(i + 16 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[2]);
        v_3 = fabs(v_base[3]);
        v_4 = fabs(v_base[4]);
        v_5 = fabs(v_base[5]);
        v_6 = fabs(v_base[6]);
        v_7 = fabs(v_base[7]);
        v_8 = fabs(v_base[8]);
        v_9 = fabs(v_base[9]);
        v_10 = fabs(v_base[10]);
        v_11 = fabs(v_base[11]);
        v_12 = fabs(v_base[12]);
        v_13 = fabs(v_base[13]);
        v_14 = fabs(v_base[14]);
        v_15 = fabs(v_base[15]);
        v_16 = fabs(v_base[16]);
        v_17 = fabs(v_base[17]);
        v_18 = fabs(v_base[18]);
        v_19 = fabs(v_base[19]);
        v_20 = fabs(v_base[20]);
        v_21 = fabs(v_base[21]);
        v_22 = fabs(v_base[22]);
        v_23 = fabs(v_base[23]);
        v_24 = fabs(v_base[24]);
        v_25 = fabs(v_base[25]);
        v_26 = fabs(v_base[26]);
        v_27 = fabs(v_base[27]);
        v_28 = fabs(v_base[28]);
        v_29 = fabs(v_base[29]);
        v_30 = fabs(v_base[30]);
        v_31 = fabs(v_base[31]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_1 = (m_1 > v_9? m_1: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_1 = (m_1 > v_11? m_1: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_1 = (m_1 > v_13? m_1: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_1 = (m_1 > v_15? m_1: v_15);
        m_0 = (m_0 > v_16? m_0: v_16);
        m_1 = (m_1 > v_17? m_1: v_17);
        m_0 = (m_0 > v_18? m_0: v_18);
        m_1 = (m_1 > v_19? m_1: v_19);
        m_0 = (m_0 > v_20? m_0: v_20);
        m_1 = (m_1 > v_21? m_1: v_21);
        m_0 = (m_0 > v_22? m_0: v_22);
        m_1 = (m_1 > v_23? m_1: v_23);
        m_0 = (m_0 > v_24? m_0: v_24);
        m_1 = (m_1 > v_25? m_1: v_25);
        m_0 = (m_0 > v_26? m_0: v_26);
        m_1 = (m_1 > v_27? m_1: v_27);
        m_0 = (m_0 > v_28? m_0: v_28);
        m_1 = (m_1 > v_29? m_1: v_29);
        m_0 = (m_0 > v_30? m_0: v_30);
        m_1 = (m_1 > v_31? m_1: v_31);
        i += 16, v_base += 32;
      }
      if(i + 8 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[2]);
        v_3 = fabs(v_base[3]);
        v_4 = fabs(v_base[4]);
        v_5 = fabs(v_base[5]);
        v_6 = fabs(v_base[6]);
        v_7 = fabs(v_base[7]);
        v_8 = fabs(v_base[8]);
        v_9 = fabs(v_base[9]);
        v_10 = fabs(v_base[10]);
        v_11 = fabs(v_base[11]);
        v_12 = fabs(v_base[12]);
        v_13 = fabs(v_base[13]);
        v_14 = fabs(v_base[14]);
        v_15 = fabs(v_base[15]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_1 = (m_1 > v_9? m_1: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_1 = (m_1 > v_11? m_1: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_1 = (m_1 > v_13? m_1: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_1 = (m_1 > v_15? m_1: v_15);
        i += 8, v_base += 16;
      }
      if(i + 4 <= n){
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
        i += 4, v_base += 8;
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

      for(i = 0; i + 32 <= n; i += 32, v_base += (incv * 64)){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[(incv * 2)]);
        v_3 = fabs(v_base[((incv * 2) + 1)]);
        v_4 = fabs(v_base[(incv * 4)]);
        v_5 = fabs(v_base[((incv * 4) + 1)]);
        v_6 = fabs(v_base[(incv * 6)]);
        v_7 = fabs(v_base[((incv * 6) + 1)]);
        v_8 = fabs(v_base[(incv * 8)]);
        v_9 = fabs(v_base[((incv * 8) + 1)]);
        v_10 = fabs(v_base[(incv * 10)]);
        v_11 = fabs(v_base[((incv * 10) + 1)]);
        v_12 = fabs(v_base[(incv * 12)]);
        v_13 = fabs(v_base[((incv * 12) + 1)]);
        v_14 = fabs(v_base[(incv * 14)]);
        v_15 = fabs(v_base[((incv * 14) + 1)]);
        v_16 = fabs(v_base[(incv * 16)]);
        v_17 = fabs(v_base[((incv * 16) + 1)]);
        v_18 = fabs(v_base[(incv * 18)]);
        v_19 = fabs(v_base[((incv * 18) + 1)]);
        v_20 = fabs(v_base[(incv * 20)]);
        v_21 = fabs(v_base[((incv * 20) + 1)]);
        v_22 = fabs(v_base[(incv * 22)]);
        v_23 = fabs(v_base[((incv * 22) + 1)]);
        v_24 = fabs(v_base[(incv * 24)]);
        v_25 = fabs(v_base[((incv * 24) + 1)]);
        v_26 = fabs(v_base[(incv * 26)]);
        v_27 = fabs(v_base[((incv * 26) + 1)]);
        v_28 = fabs(v_base[(incv * 28)]);
        v_29 = fabs(v_base[((incv * 28) + 1)]);
        v_30 = fabs(v_base[(incv * 30)]);
        v_31 = fabs(v_base[((incv * 30) + 1)]);
        v_32 = fabs(v_base[(incv * 32)]);
        v_33 = fabs(v_base[((incv * 32) + 1)]);
        v_34 = fabs(v_base[(incv * 34)]);
        v_35 = fabs(v_base[((incv * 34) + 1)]);
        v_36 = fabs(v_base[(incv * 36)]);
        v_37 = fabs(v_base[((incv * 36) + 1)]);
        v_38 = fabs(v_base[(incv * 38)]);
        v_39 = fabs(v_base[((incv * 38) + 1)]);
        v_40 = fabs(v_base[(incv * 40)]);
        v_41 = fabs(v_base[((incv * 40) + 1)]);
        v_42 = fabs(v_base[(incv * 42)]);
        v_43 = fabs(v_base[((incv * 42) + 1)]);
        v_44 = fabs(v_base[(incv * 44)]);
        v_45 = fabs(v_base[((incv * 44) + 1)]);
        v_46 = fabs(v_base[(incv * 46)]);
        v_47 = fabs(v_base[((incv * 46) + 1)]);
        v_48 = fabs(v_base[(incv * 48)]);
        v_49 = fabs(v_base[((incv * 48) + 1)]);
        v_50 = fabs(v_base[(incv * 50)]);
        v_51 = fabs(v_base[((incv * 50) + 1)]);
        v_52 = fabs(v_base[(incv * 52)]);
        v_53 = fabs(v_base[((incv * 52) + 1)]);
        v_54 = fabs(v_base[(incv * 54)]);
        v_55 = fabs(v_base[((incv * 54) + 1)]);
        v_56 = fabs(v_base[(incv * 56)]);
        v_57 = fabs(v_base[((incv * 56) + 1)]);
        v_58 = fabs(v_base[(incv * 58)]);
        v_59 = fabs(v_base[((incv * 58) + 1)]);
        v_60 = fabs(v_base[(incv * 60)]);
        v_61 = fabs(v_base[((incv * 60) + 1)]);
        v_62 = fabs(v_base[(incv * 62)]);
        v_63 = fabs(v_base[((incv * 62) + 1)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_1 = (m_1 > v_9? m_1: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_1 = (m_1 > v_11? m_1: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_1 = (m_1 > v_13? m_1: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_1 = (m_1 > v_15? m_1: v_15);
        m_0 = (m_0 > v_16? m_0: v_16);
        m_1 = (m_1 > v_17? m_1: v_17);
        m_0 = (m_0 > v_18? m_0: v_18);
        m_1 = (m_1 > v_19? m_1: v_19);
        m_0 = (m_0 > v_20? m_0: v_20);
        m_1 = (m_1 > v_21? m_1: v_21);
        m_0 = (m_0 > v_22? m_0: v_22);
        m_1 = (m_1 > v_23? m_1: v_23);
        m_0 = (m_0 > v_24? m_0: v_24);
        m_1 = (m_1 > v_25? m_1: v_25);
        m_0 = (m_0 > v_26? m_0: v_26);
        m_1 = (m_1 > v_27? m_1: v_27);
        m_0 = (m_0 > v_28? m_0: v_28);
        m_1 = (m_1 > v_29? m_1: v_29);
        m_0 = (m_0 > v_30? m_0: v_30);
        m_1 = (m_1 > v_31? m_1: v_31);
        m_0 = (m_0 > v_32? m_0: v_32);
        m_1 = (m_1 > v_33? m_1: v_33);
        m_0 = (m_0 > v_34? m_0: v_34);
        m_1 = (m_1 > v_35? m_1: v_35);
        m_0 = (m_0 > v_36? m_0: v_36);
        m_1 = (m_1 > v_37? m_1: v_37);
        m_0 = (m_0 > v_38? m_0: v_38);
        m_1 = (m_1 > v_39? m_1: v_39);
        m_0 = (m_0 > v_40? m_0: v_40);
        m_1 = (m_1 > v_41? m_1: v_41);
        m_0 = (m_0 > v_42? m_0: v_42);
        m_1 = (m_1 > v_43? m_1: v_43);
        m_0 = (m_0 > v_44? m_0: v_44);
        m_1 = (m_1 > v_45? m_1: v_45);
        m_0 = (m_0 > v_46? m_0: v_46);
        m_1 = (m_1 > v_47? m_1: v_47);
        m_0 = (m_0 > v_48? m_0: v_48);
        m_1 = (m_1 > v_49? m_1: v_49);
        m_0 = (m_0 > v_50? m_0: v_50);
        m_1 = (m_1 > v_51? m_1: v_51);
        m_0 = (m_0 > v_52? m_0: v_52);
        m_1 = (m_1 > v_53? m_1: v_53);
        m_0 = (m_0 > v_54? m_0: v_54);
        m_1 = (m_1 > v_55? m_1: v_55);
        m_0 = (m_0 > v_56? m_0: v_56);
        m_1 = (m_1 > v_57? m_1: v_57);
        m_0 = (m_0 > v_58? m_0: v_58);
        m_1 = (m_1 > v_59? m_1: v_59);
        m_0 = (m_0 > v_60? m_0: v_60);
        m_1 = (m_1 > v_61? m_1: v_61);
        m_0 = (m_0 > v_62? m_0: v_62);
        m_1 = (m_1 > v_63? m_1: v_63);
      }
      if(i + 16 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[(incv * 2)]);
        v_3 = fabs(v_base[((incv * 2) + 1)]);
        v_4 = fabs(v_base[(incv * 4)]);
        v_5 = fabs(v_base[((incv * 4) + 1)]);
        v_6 = fabs(v_base[(incv * 6)]);
        v_7 = fabs(v_base[((incv * 6) + 1)]);
        v_8 = fabs(v_base[(incv * 8)]);
        v_9 = fabs(v_base[((incv * 8) + 1)]);
        v_10 = fabs(v_base[(incv * 10)]);
        v_11 = fabs(v_base[((incv * 10) + 1)]);
        v_12 = fabs(v_base[(incv * 12)]);
        v_13 = fabs(v_base[((incv * 12) + 1)]);
        v_14 = fabs(v_base[(incv * 14)]);
        v_15 = fabs(v_base[((incv * 14) + 1)]);
        v_16 = fabs(v_base[(incv * 16)]);
        v_17 = fabs(v_base[((incv * 16) + 1)]);
        v_18 = fabs(v_base[(incv * 18)]);
        v_19 = fabs(v_base[((incv * 18) + 1)]);
        v_20 = fabs(v_base[(incv * 20)]);
        v_21 = fabs(v_base[((incv * 20) + 1)]);
        v_22 = fabs(v_base[(incv * 22)]);
        v_23 = fabs(v_base[((incv * 22) + 1)]);
        v_24 = fabs(v_base[(incv * 24)]);
        v_25 = fabs(v_base[((incv * 24) + 1)]);
        v_26 = fabs(v_base[(incv * 26)]);
        v_27 = fabs(v_base[((incv * 26) + 1)]);
        v_28 = fabs(v_base[(incv * 28)]);
        v_29 = fabs(v_base[((incv * 28) + 1)]);
        v_30 = fabs(v_base[(incv * 30)]);
        v_31 = fabs(v_base[((incv * 30) + 1)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_1 = (m_1 > v_9? m_1: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_1 = (m_1 > v_11? m_1: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_1 = (m_1 > v_13? m_1: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_1 = (m_1 > v_15? m_1: v_15);
        m_0 = (m_0 > v_16? m_0: v_16);
        m_1 = (m_1 > v_17? m_1: v_17);
        m_0 = (m_0 > v_18? m_0: v_18);
        m_1 = (m_1 > v_19? m_1: v_19);
        m_0 = (m_0 > v_20? m_0: v_20);
        m_1 = (m_1 > v_21? m_1: v_21);
        m_0 = (m_0 > v_22? m_0: v_22);
        m_1 = (m_1 > v_23? m_1: v_23);
        m_0 = (m_0 > v_24? m_0: v_24);
        m_1 = (m_1 > v_25? m_1: v_25);
        m_0 = (m_0 > v_26? m_0: v_26);
        m_1 = (m_1 > v_27? m_1: v_27);
        m_0 = (m_0 > v_28? m_0: v_28);
        m_1 = (m_1 > v_29? m_1: v_29);
        m_0 = (m_0 > v_30? m_0: v_30);
        m_1 = (m_1 > v_31? m_1: v_31);
        i += 16, v_base += (incv * 32);
      }
      if(i + 8 <= n){
        v_0 = fabs(v_base[0]);
        v_1 = fabs(v_base[1]);
        v_2 = fabs(v_base[(incv * 2)]);
        v_3 = fabs(v_base[((incv * 2) + 1)]);
        v_4 = fabs(v_base[(incv * 4)]);
        v_5 = fabs(v_base[((incv * 4) + 1)]);
        v_6 = fabs(v_base[(incv * 6)]);
        v_7 = fabs(v_base[((incv * 6) + 1)]);
        v_8 = fabs(v_base[(incv * 8)]);
        v_9 = fabs(v_base[((incv * 8) + 1)]);
        v_10 = fabs(v_base[(incv * 10)]);
        v_11 = fabs(v_base[((incv * 10) + 1)]);
        v_12 = fabs(v_base[(incv * 12)]);
        v_13 = fabs(v_base[((incv * 12) + 1)]);
        v_14 = fabs(v_base[(incv * 14)]);
        v_15 = fabs(v_base[((incv * 14) + 1)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_1 = (m_1 > v_1? m_1: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_1 = (m_1 > v_3? m_1: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_1 = (m_1 > v_5? m_1: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_1 = (m_1 > v_7? m_1: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_1 = (m_1 > v_9? m_1: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_1 = (m_1 > v_11? m_1: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_1 = (m_1 > v_13? m_1: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_1 = (m_1 > v_15? m_1: v_15);
        i += 8, v_base += (incv * 16);
      }
      if(i + 4 <= n){
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
        i += 4, v_base += (incv * 8);
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