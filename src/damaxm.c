#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  double damaxm(int n, double* v, int incv, double* y, int incy){
    __m256d mask_ABS; AVX_ABS_MASKD(mask_ABS);
    double tmp_max[4] __attribute__((aligned(32)));
    int i;
    double max;

    __m256d v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m256d y_0, y_1, y_2, y_3, y_4, y_5, y_6, y_7;
    __m256d m_0;
    m_0 = _mm256_setzero_pd();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32, y += 32){
        v_0 = _mm256_loadu_pd(v);
        v_1 = _mm256_loadu_pd(v + 4);
        v_2 = _mm256_loadu_pd(v + 8);
        v_3 = _mm256_loadu_pd(v + 12);
        v_4 = _mm256_loadu_pd(v + 16);
        v_5 = _mm256_loadu_pd(v + 20);
        v_6 = _mm256_loadu_pd(v + 24);
        v_7 = _mm256_loadu_pd(v + 28);
        y_0 = _mm256_loadu_pd(y);
        y_1 = _mm256_loadu_pd(y + 4);
        y_2 = _mm256_loadu_pd(y + 8);
        y_3 = _mm256_loadu_pd(y + 12);
        y_4 = _mm256_loadu_pd(y + 16);
        y_5 = _mm256_loadu_pd(y + 20);
        y_6 = _mm256_loadu_pd(y + 24);
        y_7 = _mm256_loadu_pd(y + 28);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(v_3, y_3), mask_ABS);
        v_4 = _mm256_and_pd(_mm256_mul_pd(v_4, y_4), mask_ABS);
        v_5 = _mm256_and_pd(_mm256_mul_pd(v_5, y_5), mask_ABS);
        v_6 = _mm256_and_pd(_mm256_mul_pd(v_6, y_6), mask_ABS);
        v_7 = _mm256_and_pd(_mm256_mul_pd(v_7, y_7), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        m_0 = _mm256_max_pd(m_0, v_4);
        m_0 = _mm256_max_pd(m_0, v_5);
        m_0 = _mm256_max_pd(m_0, v_6);
        m_0 = _mm256_max_pd(m_0, v_7);
      }
      if(i + 16 <= n){
        v_0 = _mm256_loadu_pd(v);
        v_1 = _mm256_loadu_pd(v + 4);
        v_2 = _mm256_loadu_pd(v + 8);
        v_3 = _mm256_loadu_pd(v + 12);
        y_0 = _mm256_loadu_pd(y);
        y_1 = _mm256_loadu_pd(y + 4);
        y_2 = _mm256_loadu_pd(y + 8);
        y_3 = _mm256_loadu_pd(y + 12);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(v_3, y_3), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        i += 16, v += 16, y += 16;
      }
      if(i + 8 <= n){
        v_0 = _mm256_loadu_pd(v);
        v_1 = _mm256_loadu_pd(v + 4);
        y_0 = _mm256_loadu_pd(y);
        y_1 = _mm256_loadu_pd(y + 4);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, y_1), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 8, v += 8, y += 8;
      }
      if(i + 4 <= n){
        v_0 = _mm256_loadu_pd(v);
        y_0 = _mm256_loadu_pd(y);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        i += 4, v += 4, y += 4;
      }
      if(i < n){
        v_0 = _mm256_set_pd(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
        y_0 = _mm256_set_pd(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
      }
    }else{

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32), y += (incy * 32)){
        v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        v_2 = _mm256_set_pd(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        v_3 = _mm256_set_pd(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]);
        v_4 = _mm256_set_pd(v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]);
        v_5 = _mm256_set_pd(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)]);
        v_6 = _mm256_set_pd(v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]);
        v_7 = _mm256_set_pd(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)]);
        y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        y_2 = _mm256_set_pd(y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        y_3 = _mm256_set_pd(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)]);
        y_4 = _mm256_set_pd(y[(incy * 19)], y[(incy * 18)], y[(incy * 17)], y[(incy * 16)]);
        y_5 = _mm256_set_pd(y[(incy * 23)], y[(incy * 22)], y[(incy * 21)], y[(incy * 20)]);
        y_6 = _mm256_set_pd(y[(incy * 27)], y[(incy * 26)], y[(incy * 25)], y[(incy * 24)]);
        y_7 = _mm256_set_pd(y[(incy * 31)], y[(incy * 30)], y[(incy * 29)], y[(incy * 28)]);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(v_3, y_3), mask_ABS);
        v_4 = _mm256_and_pd(_mm256_mul_pd(v_4, y_4), mask_ABS);
        v_5 = _mm256_and_pd(_mm256_mul_pd(v_5, y_5), mask_ABS);
        v_6 = _mm256_and_pd(_mm256_mul_pd(v_6, y_6), mask_ABS);
        v_7 = _mm256_and_pd(_mm256_mul_pd(v_7, y_7), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        m_0 = _mm256_max_pd(m_0, v_4);
        m_0 = _mm256_max_pd(m_0, v_5);
        m_0 = _mm256_max_pd(m_0, v_6);
        m_0 = _mm256_max_pd(m_0, v_7);
      }
      if(i + 16 <= n){
        v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        v_2 = _mm256_set_pd(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
        v_3 = _mm256_set_pd(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]);
        y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        y_2 = _mm256_set_pd(y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
        y_3 = _mm256_set_pd(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)]);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(v_3, y_3), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        i += 16, v += (incv * 16), y += (incy * 16);
      }
      if(i + 8 <= n){
        v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
        y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, y_1), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 8, v += (incv * 8), y += (incy * 8);
      }
      if(i + 4 <= n){
        v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
        y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        i += 4, v += (incv * 4), y += (incy * 4);
      }
      if(i < n){
        v_0 = _mm256_set_pd(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
        y_0 = _mm256_set_pd(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, y_0), mask_ABS);
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
  double damaxm(int n, double* v, int incv, double* y, int incy){
    __m128d mask_ABS; SSE_ABS_MASKD(mask_ABS);
    double tmp_max[2] __attribute__((aligned(16)));
    int i;
    double max;

    __m128d v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m128d y_0, y_1, y_2, y_3, y_4, y_5, y_6, y_7;
    __m128d m_0;
    m_0 = _mm_setzero_pd();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 16 <= n; i += 16, v += 16, y += 16){
        v_0 = _mm_loadu_pd(v);
        v_1 = _mm_loadu_pd(v + 2);
        v_2 = _mm_loadu_pd(v + 4);
        v_3 = _mm_loadu_pd(v + 6);
        v_4 = _mm_loadu_pd(v + 8);
        v_5 = _mm_loadu_pd(v + 10);
        v_6 = _mm_loadu_pd(v + 12);
        v_7 = _mm_loadu_pd(v + 14);
        y_0 = _mm_loadu_pd(y);
        y_1 = _mm_loadu_pd(y + 2);
        y_2 = _mm_loadu_pd(y + 4);
        y_3 = _mm_loadu_pd(y + 6);
        y_4 = _mm_loadu_pd(y + 8);
        y_5 = _mm_loadu_pd(y + 10);
        y_6 = _mm_loadu_pd(y + 12);
        y_7 = _mm_loadu_pd(y + 14);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm_and_pd(_mm_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(v_3, y_3), mask_ABS);
        v_4 = _mm_and_pd(_mm_mul_pd(v_4, y_4), mask_ABS);
        v_5 = _mm_and_pd(_mm_mul_pd(v_5, y_5), mask_ABS);
        v_6 = _mm_and_pd(_mm_mul_pd(v_6, y_6), mask_ABS);
        v_7 = _mm_and_pd(_mm_mul_pd(v_7, y_7), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        m_0 = _mm_max_pd(m_0, v_4);
        m_0 = _mm_max_pd(m_0, v_5);
        m_0 = _mm_max_pd(m_0, v_6);
        m_0 = _mm_max_pd(m_0, v_7);
      }
      if(i + 8 <= n){
        v_0 = _mm_loadu_pd(v);
        v_1 = _mm_loadu_pd(v + 2);
        v_2 = _mm_loadu_pd(v + 4);
        v_3 = _mm_loadu_pd(v + 6);
        y_0 = _mm_loadu_pd(y);
        y_1 = _mm_loadu_pd(y + 2);
        y_2 = _mm_loadu_pd(y + 4);
        y_3 = _mm_loadu_pd(y + 6);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm_and_pd(_mm_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(v_3, y_3), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        i += 8, v += 8, y += 8;
      }
      if(i + 4 <= n){
        v_0 = _mm_loadu_pd(v);
        v_1 = _mm_loadu_pd(v + 2);
        y_0 = _mm_loadu_pd(y);
        y_1 = _mm_loadu_pd(y + 2);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, y_1), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 4, v += 4, y += 4;
      }
      if(i + 2 <= n){
        v_0 = _mm_loadu_pd(v);
        y_0 = _mm_loadu_pd(y);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        i += 2, v += 2, y += 2;
      }
      if(i < n){
        v_0 = _mm_set_pd(0, v[0]);
        y_0 = _mm_set_pd(0, y[0]);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
      }
    }else{

      for(i = 0; i + 16 <= n; i += 16, v += (incv * 16), y += (incy * 16)){
        v_0 = _mm_set_pd(v[incv], v[0]);
        v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
        v_2 = _mm_set_pd(v[(incv * 5)], v[(incv * 4)]);
        v_3 = _mm_set_pd(v[(incv * 7)], v[(incv * 6)]);
        v_4 = _mm_set_pd(v[(incv * 9)], v[(incv * 8)]);
        v_5 = _mm_set_pd(v[(incv * 11)], v[(incv * 10)]);
        v_6 = _mm_set_pd(v[(incv * 13)], v[(incv * 12)]);
        v_7 = _mm_set_pd(v[(incv * 15)], v[(incv * 14)]);
        y_0 = _mm_set_pd(y[incy], y[0]);
        y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
        y_2 = _mm_set_pd(y[(incy * 5)], y[(incy * 4)]);
        y_3 = _mm_set_pd(y[(incy * 7)], y[(incy * 6)]);
        y_4 = _mm_set_pd(y[(incy * 9)], y[(incy * 8)]);
        y_5 = _mm_set_pd(y[(incy * 11)], y[(incy * 10)]);
        y_6 = _mm_set_pd(y[(incy * 13)], y[(incy * 12)]);
        y_7 = _mm_set_pd(y[(incy * 15)], y[(incy * 14)]);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm_and_pd(_mm_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(v_3, y_3), mask_ABS);
        v_4 = _mm_and_pd(_mm_mul_pd(v_4, y_4), mask_ABS);
        v_5 = _mm_and_pd(_mm_mul_pd(v_5, y_5), mask_ABS);
        v_6 = _mm_and_pd(_mm_mul_pd(v_6, y_6), mask_ABS);
        v_7 = _mm_and_pd(_mm_mul_pd(v_7, y_7), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        m_0 = _mm_max_pd(m_0, v_4);
        m_0 = _mm_max_pd(m_0, v_5);
        m_0 = _mm_max_pd(m_0, v_6);
        m_0 = _mm_max_pd(m_0, v_7);
      }
      if(i + 8 <= n){
        v_0 = _mm_set_pd(v[incv], v[0]);
        v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
        v_2 = _mm_set_pd(v[(incv * 5)], v[(incv * 4)]);
        v_3 = _mm_set_pd(v[(incv * 7)], v[(incv * 6)]);
        y_0 = _mm_set_pd(y[incy], y[0]);
        y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
        y_2 = _mm_set_pd(y[(incy * 5)], y[(incy * 4)]);
        y_3 = _mm_set_pd(y[(incy * 7)], y[(incy * 6)]);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, y_1), mask_ABS);
        v_2 = _mm_and_pd(_mm_mul_pd(v_2, y_2), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(v_3, y_3), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        i += 8, v += (incv * 8), y += (incy * 8);
      }
      if(i + 4 <= n){
        v_0 = _mm_set_pd(v[incv], v[0]);
        v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
        y_0 = _mm_set_pd(y[incy], y[0]);
        y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, y_1), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 4, v += (incv * 4), y += (incy * 4);
      }
      if(i + 2 <= n){
        v_0 = _mm_set_pd(v[incv], v[0]);
        y_0 = _mm_set_pd(y[incy], y[0]);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        i += 2, v += (incv * 2), y += (incy * 2);
      }
      if(i < n){
        v_0 = _mm_set_pd(0, v[0]);
        y_0 = _mm_set_pd(0, y[0]);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, y_0), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
      }
    }
    _mm_store_pd(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[1] ? tmp_max[0]: tmp_max[1]);
    (&max)[0] = ((double*)tmp_max)[0];
    return max;
  }
#else
  double damaxm(int n, double* v, int incv, double* y, int incy){
    int i;
    double max;

    double v_0, v_1;
    double y_0, y_1;
    double m_0;
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
  }
#endif