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
  double complex zamaxm(int n, double complex* v, int incv, double complex* y, int incy){
    __m256d mask_ABS; AVX_ABS_MASKD(mask_ABS);
    double tmp_max[4] __attribute__((aligned(32)));
    int i;
    double complex max;

    double* v_base = (double*) v;
    double* y_base = (double*) y;
    __m256d v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m256d y_0, y_1, y_2, y_3;
    __m256d m_0;
    m_0 = _mm256_setzero_pd();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 8 <= n; i += 8, v_base += 16, y_base += 16){
        v_0 = _mm256_loadu_pd(v_base);
        v_1 = _mm256_loadu_pd(v_base + 4);
        v_2 = _mm256_loadu_pd(v_base + 8);
        v_3 = _mm256_loadu_pd(v_base + 12);
        y_0 = _mm256_loadu_pd(y_base);
        y_1 = _mm256_loadu_pd(y_base + 4);
        y_2 = _mm256_loadu_pd(y_base + 8);
        y_3 = _mm256_loadu_pd(y_base + 12);
        v_4 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_5 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_1, 0b0101), _mm256_permute_pd(y_1, 0b1111)), mask_ABS);
        v_6 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_2, 0b0101), _mm256_permute_pd(y_2, 0b1111)), mask_ABS);
        v_7 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_3, 0b0101), _mm256_permute_pd(y_3, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, _mm256_permute_pd(y_1, 0b0000)), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_mul_pd(v_2, _mm256_permute_pd(y_2, 0b0000)), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(v_3, _mm256_permute_pd(y_3, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        m_0 = _mm256_max_pd(m_0, v_4);
        m_0 = _mm256_max_pd(m_0, v_5);
        m_0 = _mm256_max_pd(m_0, v_6);
        m_0 = _mm256_max_pd(m_0, v_7);
      }
      if(i + 4 <= n){
        v_0 = _mm256_loadu_pd(v_base);
        v_1 = _mm256_loadu_pd(v_base + 4);
        y_0 = _mm256_loadu_pd(y_base);
        y_1 = _mm256_loadu_pd(y_base + 4);
        v_2 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_1, 0b0101), _mm256_permute_pd(y_1, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, _mm256_permute_pd(y_1, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        i += 4, v_base += 8, y_base += 8;
      }
      if(i + 2 <= n){
        v_0 = _mm256_loadu_pd(v_base);
        y_0 = _mm256_loadu_pd(y_base);
        v_1 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 2, v_base += 4, y_base += 4;
      }
      if(i < n){
        v_0 = _mm256_set_pd(0, 0, v_base[1], v_base[0]);
        y_0 = _mm256_set_pd(0, 0, y_base[1], y_base[0]);
        v_1 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
      }
    }else{

      for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16), y_base += (incy * 16)){
        v_0 = _mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
        v_1 = _mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]);
        v_2 = _mm256_set_pd(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]);
        v_3 = _mm256_set_pd(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]);
        y_0 = _mm256_set_pd(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
        y_1 = _mm256_set_pd(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)]);
        y_2 = _mm256_set_pd(y_base[((incy * 10) + 1)], y_base[(incy * 10)], y_base[((incy * 8) + 1)], y_base[(incy * 8)]);
        y_3 = _mm256_set_pd(y_base[((incy * 14) + 1)], y_base[(incy * 14)], y_base[((incy * 12) + 1)], y_base[(incy * 12)]);
        v_4 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_5 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_1, 0b0101), _mm256_permute_pd(y_1, 0b1111)), mask_ABS);
        v_6 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_2, 0b0101), _mm256_permute_pd(y_2, 0b1111)), mask_ABS);
        v_7 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_3, 0b0101), _mm256_permute_pd(y_3, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, _mm256_permute_pd(y_1, 0b0000)), mask_ABS);
        v_2 = _mm256_and_pd(_mm256_mul_pd(v_2, _mm256_permute_pd(y_2, 0b0000)), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(v_3, _mm256_permute_pd(y_3, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        m_0 = _mm256_max_pd(m_0, v_4);
        m_0 = _mm256_max_pd(m_0, v_5);
        m_0 = _mm256_max_pd(m_0, v_6);
        m_0 = _mm256_max_pd(m_0, v_7);
      }
      if(i + 4 <= n){
        v_0 = _mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
        v_1 = _mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]);
        y_0 = _mm256_set_pd(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
        y_1 = _mm256_set_pd(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)]);
        v_2 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_3 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_1, 0b0101), _mm256_permute_pd(y_1, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        v_1 = _mm256_and_pd(_mm256_mul_pd(v_1, _mm256_permute_pd(y_1, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        m_0 = _mm256_max_pd(m_0, v_2);
        m_0 = _mm256_max_pd(m_0, v_3);
        i += 4, v_base += (incv * 8), y_base += (incy * 8);
      }
      if(i + 2 <= n){
        v_0 = _mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
        y_0 = _mm256_set_pd(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
        v_1 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
        i += 2, v_base += (incv * 4), y_base += (incy * 4);
      }
      if(i < n){
        v_0 = _mm256_set_pd(0, 0, v_base[1], v_base[0]);
        y_0 = _mm256_set_pd(0, 0, y_base[1], y_base[0]);
        v_1 = _mm256_and_pd(_mm256_mul_pd(_mm256_permute_pd(v_0, 0b0101), _mm256_permute_pd(y_0, 0b1111)), mask_ABS);
        v_0 = _mm256_and_pd(_mm256_mul_pd(v_0, _mm256_permute_pd(y_0, 0b0000)), mask_ABS);
        m_0 = _mm256_max_pd(m_0, v_0);
        m_0 = _mm256_max_pd(m_0, v_1);
      }
    }
    _mm256_store_pd(tmp_max, m_0);
    tmp_max[0] = (tmp_max[0] > tmp_max[2] ? tmp_max[0]: tmp_max[2]);
    tmp_max[1] = (tmp_max[1] > tmp_max[3] ? tmp_max[1]: tmp_max[3]);
    (&max)[0] = ((double complex*)tmp_max)[0];
    return max;
  }
#elif defined( __SSE2__ )
  double complex zamaxm(int n, double complex* v, int incv, double complex* y, int incy){
    __m128d mask_ABS; SSE_ABS_MASKD(mask_ABS);
    double tmp_max[2] __attribute__((aligned(16)));
    int i;
    double complex max;

    double* v_base = (double*) v;
    double* y_base = (double*) y;
    __m128d v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m128d y_0, y_1, y_2, y_3;
    __m128d m_0;
    m_0 = _mm_setzero_pd();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 4 <= n; i += 4, v_base += 8, y_base += 8){
        v_0 = _mm_loadu_pd(v_base);
        v_1 = _mm_loadu_pd(v_base + 2);
        v_2 = _mm_loadu_pd(v_base + 4);
        v_3 = _mm_loadu_pd(v_base + 6);
        y_0 = _mm_loadu_pd(y_base);
        y_1 = _mm_loadu_pd(y_base + 2);
        y_2 = _mm_loadu_pd(y_base + 4);
        y_3 = _mm_loadu_pd(y_base + 6);
        v_4 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_0, v_0, 0b01), _mm_shuffle_pd(y_0, y_0, 0b11)), mask_ABS);
        v_5 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_1, v_1, 0b01), _mm_shuffle_pd(y_1, y_1, 0b11)), mask_ABS);
        v_6 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_2, v_2, 0b01), _mm_shuffle_pd(y_2, y_2, 0b11)), mask_ABS);
        v_7 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_3, v_3, 0b01), _mm_shuffle_pd(y_3, y_3, 0b11)), mask_ABS);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, _mm_shuffle_pd(y_0, y_0, 0b00)), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, _mm_shuffle_pd(y_1, y_1, 0b00)), mask_ABS);
        v_2 = _mm_and_pd(_mm_mul_pd(v_2, _mm_shuffle_pd(y_2, y_2, 0b00)), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(v_3, _mm_shuffle_pd(y_3, y_3, 0b00)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        m_0 = _mm_max_pd(m_0, v_4);
        m_0 = _mm_max_pd(m_0, v_5);
        m_0 = _mm_max_pd(m_0, v_6);
        m_0 = _mm_max_pd(m_0, v_7);
      }
      if(i + 2 <= n){
        v_0 = _mm_loadu_pd(v_base);
        v_1 = _mm_loadu_pd(v_base + 2);
        y_0 = _mm_loadu_pd(y_base);
        y_1 = _mm_loadu_pd(y_base + 2);
        v_2 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_0, v_0, 0b01), _mm_shuffle_pd(y_0, y_0, 0b11)), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_1, v_1, 0b01), _mm_shuffle_pd(y_1, y_1, 0b11)), mask_ABS);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, _mm_shuffle_pd(y_0, y_0, 0b00)), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, _mm_shuffle_pd(y_1, y_1, 0b00)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        i += 2, v_base += 4, y_base += 4;
      }
      if(i + 1 <= n){
        v_0 = _mm_loadu_pd(v_base);
        y_0 = _mm_loadu_pd(y_base);
        v_1 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_0, v_0, 0b01), _mm_shuffle_pd(y_0, y_0, 0b11)), mask_ABS);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, _mm_shuffle_pd(y_0, y_0, 0b00)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 1, v_base += 2, y_base += 2;
      }
    }else{

      for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8), y_base += (incy * 8)){
        v_0 = _mm_loadu_pd(v_base);
        v_1 = _mm_loadu_pd(v_base + (incv * 2));
        v_2 = _mm_loadu_pd(v_base + (incv * 4));
        v_3 = _mm_loadu_pd(v_base + (incv * 6));
        y_0 = _mm_loadu_pd(y_base);
        y_1 = _mm_loadu_pd(y_base + (incy * 2));
        y_2 = _mm_loadu_pd(y_base + (incy * 4));
        y_3 = _mm_loadu_pd(y_base + (incy * 6));
        v_4 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_0, v_0, 0b01), _mm_shuffle_pd(y_0, y_0, 0b11)), mask_ABS);
        v_5 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_1, v_1, 0b01), _mm_shuffle_pd(y_1, y_1, 0b11)), mask_ABS);
        v_6 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_2, v_2, 0b01), _mm_shuffle_pd(y_2, y_2, 0b11)), mask_ABS);
        v_7 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_3, v_3, 0b01), _mm_shuffle_pd(y_3, y_3, 0b11)), mask_ABS);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, _mm_shuffle_pd(y_0, y_0, 0b00)), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, _mm_shuffle_pd(y_1, y_1, 0b00)), mask_ABS);
        v_2 = _mm_and_pd(_mm_mul_pd(v_2, _mm_shuffle_pd(y_2, y_2, 0b00)), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(v_3, _mm_shuffle_pd(y_3, y_3, 0b00)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        m_0 = _mm_max_pd(m_0, v_4);
        m_0 = _mm_max_pd(m_0, v_5);
        m_0 = _mm_max_pd(m_0, v_6);
        m_0 = _mm_max_pd(m_0, v_7);
      }
      if(i + 2 <= n){
        v_0 = _mm_loadu_pd(v_base);
        v_1 = _mm_loadu_pd(v_base + (incv * 2));
        y_0 = _mm_loadu_pd(y_base);
        y_1 = _mm_loadu_pd(y_base + (incy * 2));
        v_2 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_0, v_0, 0b01), _mm_shuffle_pd(y_0, y_0, 0b11)), mask_ABS);
        v_3 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_1, v_1, 0b01), _mm_shuffle_pd(y_1, y_1, 0b11)), mask_ABS);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, _mm_shuffle_pd(y_0, y_0, 0b00)), mask_ABS);
        v_1 = _mm_and_pd(_mm_mul_pd(v_1, _mm_shuffle_pd(y_1, y_1, 0b00)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        m_0 = _mm_max_pd(m_0, v_2);
        m_0 = _mm_max_pd(m_0, v_3);
        i += 2, v_base += (incv * 4), y_base += (incy * 4);
      }
      if(i + 1 <= n){
        v_0 = _mm_loadu_pd(v_base);
        y_0 = _mm_loadu_pd(y_base);
        v_1 = _mm_and_pd(_mm_mul_pd(_mm_shuffle_pd(v_0, v_0, 0b01), _mm_shuffle_pd(y_0, y_0, 0b11)), mask_ABS);
        v_0 = _mm_and_pd(_mm_mul_pd(v_0, _mm_shuffle_pd(y_0, y_0, 0b00)), mask_ABS);
        m_0 = _mm_max_pd(m_0, v_0);
        m_0 = _mm_max_pd(m_0, v_1);
        i += 1, v_base += (incv * 2), y_base += (incy * 2);
      }
    }
    _mm_store_pd(tmp_max, m_0);
    (&max)[0] = ((double complex*)tmp_max)[0];
    return max;
  }
#else
  double complex zamaxm(int n, double complex* v, int incv, double complex* y, int incy){
    int i;
    double complex max;

    double* v_base = (double*) v;
    double* y_base = (double*) y;
    double v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    double y_0, y_1, y_2, y_3;
    double m_0, m_1;
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
    ((double*)(&max))[0] = m_0;
    ((double*)(&max))[1] = m_1;
    return max;
  }
#endif