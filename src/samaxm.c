#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  float samaxm(int n, float* v, int incv, float* y, int incy){
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    float tmp_max[8] __attribute__((aligned(32)));
    int i;
    float max;
    float* max_ptr = (float*) &max;

    __m256 v_0, v_1, v_2, v_3;
    __m256 y_0, y_1, y_2, y_3;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32, y += 32){
        v_0 = _mm256_loadu_ps(v);
        v_1 = _mm256_loadu_ps(v + 8);
        v_2 = _mm256_loadu_ps(v + 16);
        v_3 = _mm256_loadu_ps(v + 24);
        y_0 = _mm256_loadu_ps(y);
        y_1 = _mm256_loadu_ps(y + 8);
        y_2 = _mm256_loadu_ps(y + 16);
        y_3 = _mm256_loadu_ps(y + 24);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_mul_ps(v_1, y_1), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_mul_ps(v_2, y_2), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_mul_ps(v_3, y_3), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
      }
      if(i + 16 <= n){
        v_0 = _mm256_loadu_ps(v);
        v_1 = _mm256_loadu_ps(v + 8);
        y_0 = _mm256_loadu_ps(y);
        y_1 = _mm256_loadu_ps(y + 8);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        i += 16, v += 16, y += 16;
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

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32), y += (incy * 32)){
        v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        v_2 = _mm256_set_ps(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)], v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]);
        v_3 = _mm256_set_ps(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)], v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]);
        y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm256_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)], y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        y_2 = _mm256_set_ps(y[(incy * 23)], y[(incy * 22)], y[(incy * 21)], y[(incy * 20)], y[(incy * 19)], y[(incy * 18)], y[(incy * 17)], y[(incy * 16)]);
        y_3 = _mm256_set_ps(y[(incy * 31)], y[(incy * 30)], y[(incy * 29)], y[(incy * 28)], y[(incy * 27)], y[(incy * 26)], y[(incy * 25)], y[(incy * 24)]);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_mul_ps(v_1, y_1), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_mul_ps(v_2, y_2), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_mul_ps(v_3, y_3), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
      }
      if(i + 16 <= n){
        v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm256_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)], y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        v_0 = _mm256_and_ps(_mm256_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        i += 16, v += (incv * 16), y += (incy * 16);
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
    max_ptr[0] = ((float*)tmp_max)[0];
    return max;
  }
#elif defined( __SSE2__ )
  float samaxm(int n, float* v, int incv, float* y, int incy){
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    float tmp_max[4] __attribute__((aligned(16)));
    int i;
    float max;
    float* max_ptr = (float*) &max;

    __m128 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m128 y_0, y_1, y_2, y_3, y_4, y_5, y_6, y_7;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32, y += 32){
        v_0 = _mm_loadu_ps(v);
        v_1 = _mm_loadu_ps(v + 4);
        v_2 = _mm_loadu_ps(v + 8);
        v_3 = _mm_loadu_ps(v + 12);
        v_4 = _mm_loadu_ps(v + 16);
        v_5 = _mm_loadu_ps(v + 20);
        v_6 = _mm_loadu_ps(v + 24);
        v_7 = _mm_loadu_ps(v + 28);
        y_0 = _mm_loadu_ps(y);
        y_1 = _mm_loadu_ps(y + 4);
        y_2 = _mm_loadu_ps(y + 8);
        y_3 = _mm_loadu_ps(y + 12);
        y_4 = _mm_loadu_ps(y + 16);
        y_5 = _mm_loadu_ps(y + 20);
        y_6 = _mm_loadu_ps(y + 24);
        y_7 = _mm_loadu_ps(y + 28);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        v_2 = _mm_and_ps(_mm_mul_ps(v_2, y_2), mask_ABS);
        v_3 = _mm_and_ps(_mm_mul_ps(v_3, y_3), mask_ABS);
        v_4 = _mm_and_ps(_mm_mul_ps(v_4, y_4), mask_ABS);
        v_5 = _mm_and_ps(_mm_mul_ps(v_5, y_5), mask_ABS);
        v_6 = _mm_and_ps(_mm_mul_ps(v_6, y_6), mask_ABS);
        v_7 = _mm_and_ps(_mm_mul_ps(v_7, y_7), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        m_0 = _mm_max_ps(m_0, v_4);
        m_0 = _mm_max_ps(m_0, v_5);
        m_0 = _mm_max_ps(m_0, v_6);
        m_0 = _mm_max_ps(m_0, v_7);
      }
      if(i + 16 <= n){
        v_0 = _mm_loadu_ps(v);
        v_1 = _mm_loadu_ps(v + 4);
        v_2 = _mm_loadu_ps(v + 8);
        v_3 = _mm_loadu_ps(v + 12);
        y_0 = _mm_loadu_ps(y);
        y_1 = _mm_loadu_ps(y + 4);
        y_2 = _mm_loadu_ps(y + 8);
        y_3 = _mm_loadu_ps(y + 12);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        v_2 = _mm_and_ps(_mm_mul_ps(v_2, y_2), mask_ABS);
        v_3 = _mm_and_ps(_mm_mul_ps(v_3, y_3), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        i += 16, v += 16, y += 16;
      }
      if(i + 8 <= n){
        v_0 = _mm_loadu_ps(v);
        v_1 = _mm_loadu_ps(v + 4);
        y_0 = _mm_loadu_ps(y);
        y_1 = _mm_loadu_ps(y + 4);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        i += 8, v += 8, y += 8;
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

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32), y += (incy * 32)){
        v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        v_2 = _mm_set_ps(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        v_3 = _mm_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]);
        v_4 = _mm_set_ps(v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]);
        v_5 = _mm_set_ps(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)]);
        v_6 = _mm_set_ps(v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]);
        v_7 = _mm_set_ps(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)]);
        y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        y_2 = _mm_set_ps(y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        y_3 = _mm_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)]);
        y_4 = _mm_set_ps(y[(incy * 19)], y[(incy * 18)], y[(incy * 17)], y[(incy * 16)]);
        y_5 = _mm_set_ps(y[(incy * 23)], y[(incy * 22)], y[(incy * 21)], y[(incy * 20)]);
        y_6 = _mm_set_ps(y[(incy * 27)], y[(incy * 26)], y[(incy * 25)], y[(incy * 24)]);
        y_7 = _mm_set_ps(y[(incy * 31)], y[(incy * 30)], y[(incy * 29)], y[(incy * 28)]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        v_2 = _mm_and_ps(_mm_mul_ps(v_2, y_2), mask_ABS);
        v_3 = _mm_and_ps(_mm_mul_ps(v_3, y_3), mask_ABS);
        v_4 = _mm_and_ps(_mm_mul_ps(v_4, y_4), mask_ABS);
        v_5 = _mm_and_ps(_mm_mul_ps(v_5, y_5), mask_ABS);
        v_6 = _mm_and_ps(_mm_mul_ps(v_6, y_6), mask_ABS);
        v_7 = _mm_and_ps(_mm_mul_ps(v_7, y_7), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        m_0 = _mm_max_ps(m_0, v_4);
        m_0 = _mm_max_ps(m_0, v_5);
        m_0 = _mm_max_ps(m_0, v_6);
        m_0 = _mm_max_ps(m_0, v_7);
      }
      if(i + 16 <= n){
        v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        v_2 = _mm_set_ps(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        v_3 = _mm_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]);
        y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        y_2 = _mm_set_ps(y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        y_3 = _mm_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        v_2 = _mm_and_ps(_mm_mul_ps(v_2, y_2), mask_ABS);
        v_3 = _mm_and_ps(_mm_mul_ps(v_3, y_3), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        i += 16, v += (incv * 16), y += (incy * 16);
      }
      if(i + 8 <= n){
        v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        v_0 = _mm_and_ps(_mm_mul_ps(v_0, y_0), mask_ABS);
        v_1 = _mm_and_ps(_mm_mul_ps(v_1, y_1), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        i += 8, v += (incv * 8), y += (incy * 8);
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
    max_ptr[0] = ((float*)tmp_max)[0];
    return max;
  }
#else
  float samaxm(int n, float* v, int incv, float* y, int incy){
    int i;
    float max;
    float* max_ptr = (float*) &max;

    float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11, v_12, v_13, v_14, v_15, v_16, v_17, v_18, v_19, v_20, v_21, v_22, v_23, v_24, v_25, v_26, v_27, v_28, v_29, v_30, v_31;
    float y_0, y_1, y_2, y_3, y_4, y_5, y_6, y_7, y_8, y_9, y_10, y_11, y_12, y_13, y_14, y_15, y_16, y_17, y_18, y_19, y_20, y_21, y_22, y_23, y_24, y_25, y_26, y_27, y_28, y_29, y_30, y_31;
    float m_0;
    m_0 = 0;

    if(incv == 1 && incy == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32, y += 32){
        v_0 = v[0];
        v_1 = v[1];
        v_2 = v[2];
        v_3 = v[3];
        v_4 = v[4];
        v_5 = v[5];
        v_6 = v[6];
        v_7 = v[7];
        v_8 = v[8];
        v_9 = v[9];
        v_10 = v[10];
        v_11 = v[11];
        v_12 = v[12];
        v_13 = v[13];
        v_14 = v[14];
        v_15 = v[15];
        v_16 = v[16];
        v_17 = v[17];
        v_18 = v[18];
        v_19 = v[19];
        v_20 = v[20];
        v_21 = v[21];
        v_22 = v[22];
        v_23 = v[23];
        v_24 = v[24];
        v_25 = v[25];
        v_26 = v[26];
        v_27 = v[27];
        v_28 = v[28];
        v_29 = v[29];
        v_30 = v[30];
        v_31 = v[31];
        y_0 = y[0];
        y_1 = y[1];
        y_2 = y[2];
        y_3 = y[3];
        y_4 = y[4];
        y_5 = y[5];
        y_6 = y[6];
        y_7 = y[7];
        y_8 = y[8];
        y_9 = y[9];
        y_10 = y[10];
        y_11 = y[11];
        y_12 = y[12];
        y_13 = y[13];
        y_14 = y[14];
        y_15 = y[15];
        y_16 = y[16];
        y_17 = y[17];
        y_18 = y[18];
        y_19 = y[19];
        y_20 = y[20];
        y_21 = y[21];
        y_22 = y[22];
        y_23 = y[23];
        y_24 = y[24];
        y_25 = y[25];
        y_26 = y[26];
        y_27 = y[27];
        y_28 = y[28];
        y_29 = y[29];
        y_30 = y[30];
        y_31 = y[31];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        v_4 = fabs(v_4 * y_4);
        v_5 = fabs(v_5 * y_5);
        v_6 = fabs(v_6 * y_6);
        v_7 = fabs(v_7 * y_7);
        v_8 = fabs(v_8 * y_8);
        v_9 = fabs(v_9 * y_9);
        v_10 = fabs(v_10 * y_10);
        v_11 = fabs(v_11 * y_11);
        v_12 = fabs(v_12 * y_12);
        v_13 = fabs(v_13 * y_13);
        v_14 = fabs(v_14 * y_14);
        v_15 = fabs(v_15 * y_15);
        v_16 = fabs(v_16 * y_16);
        v_17 = fabs(v_17 * y_17);
        v_18 = fabs(v_18 * y_18);
        v_19 = fabs(v_19 * y_19);
        v_20 = fabs(v_20 * y_20);
        v_21 = fabs(v_21 * y_21);
        v_22 = fabs(v_22 * y_22);
        v_23 = fabs(v_23 * y_23);
        v_24 = fabs(v_24 * y_24);
        v_25 = fabs(v_25 * y_25);
        v_26 = fabs(v_26 * y_26);
        v_27 = fabs(v_27 * y_27);
        v_28 = fabs(v_28 * y_28);
        v_29 = fabs(v_29 * y_29);
        v_30 = fabs(v_30 * y_30);
        v_31 = fabs(v_31 * y_31);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_0 = (m_0 > v_9? m_0: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_0 = (m_0 > v_11? m_0: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_0 = (m_0 > v_13? m_0: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_0 = (m_0 > v_15? m_0: v_15);
        m_0 = (m_0 > v_16? m_0: v_16);
        m_0 = (m_0 > v_17? m_0: v_17);
        m_0 = (m_0 > v_18? m_0: v_18);
        m_0 = (m_0 > v_19? m_0: v_19);
        m_0 = (m_0 > v_20? m_0: v_20);
        m_0 = (m_0 > v_21? m_0: v_21);
        m_0 = (m_0 > v_22? m_0: v_22);
        m_0 = (m_0 > v_23? m_0: v_23);
        m_0 = (m_0 > v_24? m_0: v_24);
        m_0 = (m_0 > v_25? m_0: v_25);
        m_0 = (m_0 > v_26? m_0: v_26);
        m_0 = (m_0 > v_27? m_0: v_27);
        m_0 = (m_0 > v_28? m_0: v_28);
        m_0 = (m_0 > v_29? m_0: v_29);
        m_0 = (m_0 > v_30? m_0: v_30);
        m_0 = (m_0 > v_31? m_0: v_31);
      }
      if(i + 16 <= n){
        v_0 = v[0];
        v_1 = v[1];
        v_2 = v[2];
        v_3 = v[3];
        v_4 = v[4];
        v_5 = v[5];
        v_6 = v[6];
        v_7 = v[7];
        v_8 = v[8];
        v_9 = v[9];
        v_10 = v[10];
        v_11 = v[11];
        v_12 = v[12];
        v_13 = v[13];
        v_14 = v[14];
        v_15 = v[15];
        y_0 = y[0];
        y_1 = y[1];
        y_2 = y[2];
        y_3 = y[3];
        y_4 = y[4];
        y_5 = y[5];
        y_6 = y[6];
        y_7 = y[7];
        y_8 = y[8];
        y_9 = y[9];
        y_10 = y[10];
        y_11 = y[11];
        y_12 = y[12];
        y_13 = y[13];
        y_14 = y[14];
        y_15 = y[15];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        v_4 = fabs(v_4 * y_4);
        v_5 = fabs(v_5 * y_5);
        v_6 = fabs(v_6 * y_6);
        v_7 = fabs(v_7 * y_7);
        v_8 = fabs(v_8 * y_8);
        v_9 = fabs(v_9 * y_9);
        v_10 = fabs(v_10 * y_10);
        v_11 = fabs(v_11 * y_11);
        v_12 = fabs(v_12 * y_12);
        v_13 = fabs(v_13 * y_13);
        v_14 = fabs(v_14 * y_14);
        v_15 = fabs(v_15 * y_15);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_0 = (m_0 > v_9? m_0: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_0 = (m_0 > v_11? m_0: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_0 = (m_0 > v_13? m_0: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_0 = (m_0 > v_15? m_0: v_15);
        i += 16, v += 16, y += 16;
      }
      if(i + 8 <= n){
        v_0 = v[0];
        v_1 = v[1];
        v_2 = v[2];
        v_3 = v[3];
        v_4 = v[4];
        v_5 = v[5];
        v_6 = v[6];
        v_7 = v[7];
        y_0 = y[0];
        y_1 = y[1];
        y_2 = y[2];
        y_3 = y[3];
        y_4 = y[4];
        y_5 = y[5];
        y_6 = y[6];
        y_7 = y[7];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        v_4 = fabs(v_4 * y_4);
        v_5 = fabs(v_5 * y_5);
        v_6 = fabs(v_6 * y_6);
        v_7 = fabs(v_7 * y_7);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        i += 8, v += 8, y += 8;
      }
      if(i + 4 <= n){
        v_0 = v[0];
        v_1 = v[1];
        v_2 = v[2];
        v_3 = v[3];
        y_0 = y[0];
        y_1 = y[1];
        y_2 = y[2];
        y_3 = y[3];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        i += 4, v += 4, y += 4;
      }
      if(i + 2 <= n){
        v_0 = v[0];
        v_1 = v[1];
        y_0 = y[0];
        y_1 = y[1];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        i += 2, v += 2, y += 2;
      }
      if(i + 1 <= n){
        v_0 = v[0];
        y_0 = y[0];
        v_0 = fabs(v_0 * y_0);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += 1, y += 1;
      }
    }else{

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32), y += (incy * 32)){
        v_0 = v[0];
        v_1 = v[incv];
        v_2 = v[(incv * 2)];
        v_3 = v[(incv * 3)];
        v_4 = v[(incv * 4)];
        v_5 = v[(incv * 5)];
        v_6 = v[(incv * 6)];
        v_7 = v[(incv * 7)];
        v_8 = v[(incv * 8)];
        v_9 = v[(incv * 9)];
        v_10 = v[(incv * 10)];
        v_11 = v[(incv * 11)];
        v_12 = v[(incv * 12)];
        v_13 = v[(incv * 13)];
        v_14 = v[(incv * 14)];
        v_15 = v[(incv * 15)];
        v_16 = v[(incv * 16)];
        v_17 = v[(incv * 17)];
        v_18 = v[(incv * 18)];
        v_19 = v[(incv * 19)];
        v_20 = v[(incv * 20)];
        v_21 = v[(incv * 21)];
        v_22 = v[(incv * 22)];
        v_23 = v[(incv * 23)];
        v_24 = v[(incv * 24)];
        v_25 = v[(incv * 25)];
        v_26 = v[(incv * 26)];
        v_27 = v[(incv * 27)];
        v_28 = v[(incv * 28)];
        v_29 = v[(incv * 29)];
        v_30 = v[(incv * 30)];
        v_31 = v[(incv * 31)];
        y_0 = y[0];
        y_1 = y[incy];
        y_2 = y[(incy * 2)];
        y_3 = y[(incy * 3)];
        y_4 = y[(incy * 4)];
        y_5 = y[(incy * 5)];
        y_6 = y[(incy * 6)];
        y_7 = y[(incy * 7)];
        y_8 = y[(incy * 8)];
        y_9 = y[(incy * 9)];
        y_10 = y[(incy * 10)];
        y_11 = y[(incy * 11)];
        y_12 = y[(incy * 12)];
        y_13 = y[(incy * 13)];
        y_14 = y[(incy * 14)];
        y_15 = y[(incy * 15)];
        y_16 = y[(incy * 16)];
        y_17 = y[(incy * 17)];
        y_18 = y[(incy * 18)];
        y_19 = y[(incy * 19)];
        y_20 = y[(incy * 20)];
        y_21 = y[(incy * 21)];
        y_22 = y[(incy * 22)];
        y_23 = y[(incy * 23)];
        y_24 = y[(incy * 24)];
        y_25 = y[(incy * 25)];
        y_26 = y[(incy * 26)];
        y_27 = y[(incy * 27)];
        y_28 = y[(incy * 28)];
        y_29 = y[(incy * 29)];
        y_30 = y[(incy * 30)];
        y_31 = y[(incy * 31)];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        v_4 = fabs(v_4 * y_4);
        v_5 = fabs(v_5 * y_5);
        v_6 = fabs(v_6 * y_6);
        v_7 = fabs(v_7 * y_7);
        v_8 = fabs(v_8 * y_8);
        v_9 = fabs(v_9 * y_9);
        v_10 = fabs(v_10 * y_10);
        v_11 = fabs(v_11 * y_11);
        v_12 = fabs(v_12 * y_12);
        v_13 = fabs(v_13 * y_13);
        v_14 = fabs(v_14 * y_14);
        v_15 = fabs(v_15 * y_15);
        v_16 = fabs(v_16 * y_16);
        v_17 = fabs(v_17 * y_17);
        v_18 = fabs(v_18 * y_18);
        v_19 = fabs(v_19 * y_19);
        v_20 = fabs(v_20 * y_20);
        v_21 = fabs(v_21 * y_21);
        v_22 = fabs(v_22 * y_22);
        v_23 = fabs(v_23 * y_23);
        v_24 = fabs(v_24 * y_24);
        v_25 = fabs(v_25 * y_25);
        v_26 = fabs(v_26 * y_26);
        v_27 = fabs(v_27 * y_27);
        v_28 = fabs(v_28 * y_28);
        v_29 = fabs(v_29 * y_29);
        v_30 = fabs(v_30 * y_30);
        v_31 = fabs(v_31 * y_31);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_0 = (m_0 > v_9? m_0: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_0 = (m_0 > v_11? m_0: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_0 = (m_0 > v_13? m_0: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_0 = (m_0 > v_15? m_0: v_15);
        m_0 = (m_0 > v_16? m_0: v_16);
        m_0 = (m_0 > v_17? m_0: v_17);
        m_0 = (m_0 > v_18? m_0: v_18);
        m_0 = (m_0 > v_19? m_0: v_19);
        m_0 = (m_0 > v_20? m_0: v_20);
        m_0 = (m_0 > v_21? m_0: v_21);
        m_0 = (m_0 > v_22? m_0: v_22);
        m_0 = (m_0 > v_23? m_0: v_23);
        m_0 = (m_0 > v_24? m_0: v_24);
        m_0 = (m_0 > v_25? m_0: v_25);
        m_0 = (m_0 > v_26? m_0: v_26);
        m_0 = (m_0 > v_27? m_0: v_27);
        m_0 = (m_0 > v_28? m_0: v_28);
        m_0 = (m_0 > v_29? m_0: v_29);
        m_0 = (m_0 > v_30? m_0: v_30);
        m_0 = (m_0 > v_31? m_0: v_31);
      }
      if(i + 16 <= n){
        v_0 = v[0];
        v_1 = v[incv];
        v_2 = v[(incv * 2)];
        v_3 = v[(incv * 3)];
        v_4 = v[(incv * 4)];
        v_5 = v[(incv * 5)];
        v_6 = v[(incv * 6)];
        v_7 = v[(incv * 7)];
        v_8 = v[(incv * 8)];
        v_9 = v[(incv * 9)];
        v_10 = v[(incv * 10)];
        v_11 = v[(incv * 11)];
        v_12 = v[(incv * 12)];
        v_13 = v[(incv * 13)];
        v_14 = v[(incv * 14)];
        v_15 = v[(incv * 15)];
        y_0 = y[0];
        y_1 = y[incy];
        y_2 = y[(incy * 2)];
        y_3 = y[(incy * 3)];
        y_4 = y[(incy * 4)];
        y_5 = y[(incy * 5)];
        y_6 = y[(incy * 6)];
        y_7 = y[(incy * 7)];
        y_8 = y[(incy * 8)];
        y_9 = y[(incy * 9)];
        y_10 = y[(incy * 10)];
        y_11 = y[(incy * 11)];
        y_12 = y[(incy * 12)];
        y_13 = y[(incy * 13)];
        y_14 = y[(incy * 14)];
        y_15 = y[(incy * 15)];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        v_4 = fabs(v_4 * y_4);
        v_5 = fabs(v_5 * y_5);
        v_6 = fabs(v_6 * y_6);
        v_7 = fabs(v_7 * y_7);
        v_8 = fabs(v_8 * y_8);
        v_9 = fabs(v_9 * y_9);
        v_10 = fabs(v_10 * y_10);
        v_11 = fabs(v_11 * y_11);
        v_12 = fabs(v_12 * y_12);
        v_13 = fabs(v_13 * y_13);
        v_14 = fabs(v_14 * y_14);
        v_15 = fabs(v_15 * y_15);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        m_0 = (m_0 > v_8? m_0: v_8);
        m_0 = (m_0 > v_9? m_0: v_9);
        m_0 = (m_0 > v_10? m_0: v_10);
        m_0 = (m_0 > v_11? m_0: v_11);
        m_0 = (m_0 > v_12? m_0: v_12);
        m_0 = (m_0 > v_13? m_0: v_13);
        m_0 = (m_0 > v_14? m_0: v_14);
        m_0 = (m_0 > v_15? m_0: v_15);
        i += 16, v += (incv * 16), y += (incy * 16);
      }
      if(i + 8 <= n){
        v_0 = v[0];
        v_1 = v[incv];
        v_2 = v[(incv * 2)];
        v_3 = v[(incv * 3)];
        v_4 = v[(incv * 4)];
        v_5 = v[(incv * 5)];
        v_6 = v[(incv * 6)];
        v_7 = v[(incv * 7)];
        y_0 = y[0];
        y_1 = y[incy];
        y_2 = y[(incy * 2)];
        y_3 = y[(incy * 3)];
        y_4 = y[(incy * 4)];
        y_5 = y[(incy * 5)];
        y_6 = y[(incy * 6)];
        y_7 = y[(incy * 7)];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        v_4 = fabs(v_4 * y_4);
        v_5 = fabs(v_5 * y_5);
        v_6 = fabs(v_6 * y_6);
        v_7 = fabs(v_7 * y_7);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        i += 8, v += (incv * 8), y += (incy * 8);
      }
      if(i + 4 <= n){
        v_0 = v[0];
        v_1 = v[incv];
        v_2 = v[(incv * 2)];
        v_3 = v[(incv * 3)];
        y_0 = y[0];
        y_1 = y[incy];
        y_2 = y[(incy * 2)];
        y_3 = y[(incy * 3)];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        v_2 = fabs(v_2 * y_2);
        v_3 = fabs(v_3 * y_3);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        i += 4, v += (incv * 4), y += (incy * 4);
      }
      if(i + 2 <= n){
        v_0 = v[0];
        v_1 = v[incv];
        y_0 = y[0];
        y_1 = y[incy];
        v_0 = fabs(v_0 * y_0);
        v_1 = fabs(v_1 * y_1);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        i += 2, v += (incv * 2), y += (incy * 2);
      }
      if(i + 1 <= n){
        v_0 = v[0];
        y_0 = y[0];
        v_0 = fabs(v_0 * y_0);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += incv, y += incy;
      }
    }
    max_ptr[0] = m_0;
    return max;
  }
#endif