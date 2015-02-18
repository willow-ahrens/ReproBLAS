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
  double complex zamaxm(int n, double complex* v, int incv, double complex* y, int incy){
    /*[[[cog
    cog.out(generate.generate(amax.AMaxM(dataTypes.DoubleComplex, vectorizations.AVX), args, params))
    ]]]*/
    __m256d mask_ABS; AVX_ABS_MASKD(mask_ABS);
    double tmp_max[4] __attribute__((aligned(32)));
    int i;
    double complex max;

    double* v_base = (double*) v;
    double* y_base = (double*) y;
    __m256d v_0, v_1, v_2, v_3;
    __m256d y_0, y_1;
    __m256d m_0;
    m_0 = _mm256_setzero_pd();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 4 <= n; i += 4, v_base += 8, y_base += 8){
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

      for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8), y_base += (incy * 8)){
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
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  double complex zamaxm(int n, double complex* v, int incv, double complex* y, int incy){
    /*[[[cog
    cog.out(generate.generate(amax.AMaxM(dataTypes.DoubleComplex, vectorizations.SSE), args, params))
    ]]]*/
    __m128d mask_ABS; SSE_ABS_MASKD(mask_ABS);
    double tmp_max[2] __attribute__((aligned(16)));
    int i;
    double complex max;

    double* v_base = (double*) v;
    double* y_base = (double*) y;
    __m128d v_0, v_1, v_2, v_3;
    __m128d y_0, y_1;
    __m128d m_0;
    m_0 = _mm_setzero_pd();

    if(incv == 1 && incy == 1){

      for(i = 0; i + 2 <= n; i += 2, v_base += 4, y_base += 4){
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

      for(i = 0; i + 2 <= n; i += 2, v_base += (incv * 4), y_base += (incy * 4)){
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
    //[[[end]]]
  }
#else
  double complex zamaxm(int n, double complex* v, int incv, double complex* y, int incy){
    /*[[[cog
    cog.out(generate.generate(amax.AMaxM(dataTypes.DoubleComplex, vectorizations.SISD), args, params))
    ]]]*/
    int i;
    double complex max;

    double* v_base = (double*) v;
    double* y_base = (double*) y;
    double v_0, v_1, v_2, v_3;
    double y_0, y_1;
    double m_0, m_1;
    m_0 = 0;
    m_1 = 0;

    if(incv == 1 && incy == 1){

      for(i = 0; i + 1 <= n; i += 1, v_base += 2, y_base += 2){
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
      }
    }else{

      for(i = 0; i + 1 <= n; i += 1, v_base += (incv * 2), y_base += (incy * 2)){
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
      }
    }
    ((double*)(&max))[0] = m_0;
    ((double*)(&max))[1] = m_1;
    return max;
    //[[[end]]]
  }
#endif
