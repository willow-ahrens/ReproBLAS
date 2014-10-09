#include <complex.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <emmintrin.h>


#if defined( __AVX__ )
  void cdotuI2(int n, float complex* v, int incv, float complex* y, int incy, int fold, float complex* sum){
    __m256 mask_NCONJ; AVX_NCONJ_MASKS(mask_NCONJ);
    __m256 mask_BLP; AVX_BLP_MASKS(mask_BLP);
    float complex tmp[4] __attribute__((aligned(32)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float* y_base = (float*) y;
        __m256 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        __m256 y_0, y_1, y_2, y_3;
        __m256 q_0;
        __m256 s_0_0;
        __m256 s_1_0;
        __m256 s_2_0;

        s_0_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base));
        s_1_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base + 2));
        s_2_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base + 4));
        if(incv == 1 && incy == 1){

          for(i = 0; i + 16 <= n; i += 16, v_base += 32, y_base += 32){
            v_0 = _mm256_loadu_ps(v_base);
            v_1 = _mm256_loadu_ps(v_base + 8);
            v_2 = _mm256_loadu_ps(v_base + 16);
            v_3 = _mm256_loadu_ps(v_base + 24);
            y_0 = _mm256_loadu_ps(y_base);
            y_1 = _mm256_loadu_ps(y_base + 8);
            y_2 = _mm256_loadu_ps(y_base + 16);
            y_3 = _mm256_loadu_ps(y_base + 24);
            v_4 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_5 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_1, 0b10110001), _mm256_permute_ps(y_1, 0b11110101)), mask_NCONJ);
            v_6 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_2, 0b10110001), _mm256_permute_ps(y_2, 0b11110101)), mask_NCONJ);
            v_7 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_3, 0b10110001), _mm256_permute_ps(y_3, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            v_1 = _mm256_mul_ps(v_1, _mm256_permute_ps(y_1, 0b10100000));
            v_2 = _mm256_mul_ps(v_2, _mm256_permute_ps(y_2, 0b10100000));
            v_3 = _mm256_mul_ps(v_3, _mm256_permute_ps(y_3, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_4, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_4 = _mm256_add_ps(v_4, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_4, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_4 = _mm256_add_ps(v_4, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_4, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_5 = _mm256_add_ps(v_5, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_5 = _mm256_add_ps(v_5, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_5, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_6, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_6 = _mm256_add_ps(v_6, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_6, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_6 = _mm256_add_ps(v_6, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_6, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_7 = _mm256_add_ps(v_7, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_7 = _mm256_add_ps(v_7, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_7, mask_BLP));
          }
          if(i + 8 <= n){
            v_0 = _mm256_loadu_ps(v_base);
            v_1 = _mm256_loadu_ps(v_base + 8);
            y_0 = _mm256_loadu_ps(y_base);
            y_1 = _mm256_loadu_ps(y_base + 8);
            v_2 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_1, 0b10110001), _mm256_permute_ps(y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            v_1 = _mm256_mul_ps(v_1, _mm256_permute_ps(y_1, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_3, mask_BLP));
            i += 8, v_base += 16, y_base += 16;
          }
          if(i + 4 <= n){
            v_0 = _mm256_loadu_ps(v_base);
            y_0 = _mm256_loadu_ps(y_base);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
            i += 4, v_base += 8, y_base += 8;
          }
          if(i < n){
            v_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[2]:0, (n - i)>1?((double*)v_base)[1]:0, ((double*)v_base)[0]);
            y_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)y_base)[2]:0, (n - i)>1?((double*)y_base)[1]:0, ((double*)y_base)[0]);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
          }
        }else{

          for(i = 0; i + 16 <= n; i += 16, v_base += (incv * 32), y_base += (incy * 32)){
            v_0 = _mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            v_1 = _mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]);
            v_2 = _mm256_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)], v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]);
            v_3 = _mm256_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)], v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]);
            y_0 = _mm256_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)], y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            y_1 = _mm256_set_ps(y_base[((incy * 14) + 1)], y_base[(incy * 14)], y_base[((incy * 12) + 1)], y_base[(incy * 12)], y_base[((incy * 10) + 1)], y_base[(incy * 10)], y_base[((incy * 8) + 1)], y_base[(incy * 8)]);
            y_2 = _mm256_set_ps(y_base[((incy * 22) + 1)], y_base[(incy * 22)], y_base[((incy * 20) + 1)], y_base[(incy * 20)], y_base[((incy * 18) + 1)], y_base[(incy * 18)], y_base[((incy * 16) + 1)], y_base[(incy * 16)]);
            y_3 = _mm256_set_ps(y_base[((incy * 30) + 1)], y_base[(incy * 30)], y_base[((incy * 28) + 1)], y_base[(incy * 28)], y_base[((incy * 26) + 1)], y_base[(incy * 26)], y_base[((incy * 24) + 1)], y_base[(incy * 24)]);
            v_4 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_5 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_1, 0b10110001), _mm256_permute_ps(y_1, 0b11110101)), mask_NCONJ);
            v_6 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_2, 0b10110001), _mm256_permute_ps(y_2, 0b11110101)), mask_NCONJ);
            v_7 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_3, 0b10110001), _mm256_permute_ps(y_3, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            v_1 = _mm256_mul_ps(v_1, _mm256_permute_ps(y_1, 0b10100000));
            v_2 = _mm256_mul_ps(v_2, _mm256_permute_ps(y_2, 0b10100000));
            v_3 = _mm256_mul_ps(v_3, _mm256_permute_ps(y_3, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_4, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_4 = _mm256_add_ps(v_4, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_4, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_4 = _mm256_add_ps(v_4, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_4, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_5 = _mm256_add_ps(v_5, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_5 = _mm256_add_ps(v_5, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_5, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_6, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_6 = _mm256_add_ps(v_6, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_6, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_6 = _mm256_add_ps(v_6, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_6, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_7 = _mm256_add_ps(v_7, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_7 = _mm256_add_ps(v_7, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_7, mask_BLP));
          }
          if(i + 8 <= n){
            v_0 = _mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            v_1 = _mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]);
            y_0 = _mm256_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)], y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            y_1 = _mm256_set_ps(y_base[((incy * 14) + 1)], y_base[(incy * 14)], y_base[((incy * 12) + 1)], y_base[(incy * 12)], y_base[((incy * 10) + 1)], y_base[(incy * 10)], y_base[((incy * 8) + 1)], y_base[(incy * 8)]);
            v_2 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_1, 0b10110001), _mm256_permute_ps(y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            v_1 = _mm256_mul_ps(v_1, _mm256_permute_ps(y_1, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_2 = _mm256_add_ps(v_2, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_3 = _mm256_add_ps(v_3, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_3, mask_BLP));
            i += 8, v_base += (incv * 16), y_base += (incy * 16);
          }
          if(i + 4 <= n){
            v_0 = _mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            y_0 = _mm256_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)], y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
            i += 4, v_base += (incv * 8), y_base += (incy * 8);
          }
          if(i < n){
            v_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[(incv * 2)]:0, (n - i)>1?((double*)v_base)[incv]:0, ((double*)v_base)[0]);
            y_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)y_base)[(incy * 2)]:0, (n - i)>1?((double*)y_base)[incy]:0, ((double*)y_base)[0]);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_1 = _mm256_add_ps(v_1, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_1, mask_BLP));
          }
        }
        s_0_0 = _mm256_sub_ps(s_0_0, _mm256_set_ps(sum_base[1], sum_base[0], sum_base[1], sum_base[0], sum_base[1], sum_base[0], 0, 0));
        _mm256_store_ps((float*)tmp, s_0_0);
        sum[0] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        s_1_0 = _mm256_sub_ps(s_1_0, _mm256_set_ps(sum_base[3], sum_base[2], sum_base[3], sum_base[2], sum_base[3], sum_base[2], 0, 0));
        _mm256_store_ps((float*)tmp, s_1_0);
        sum[1] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        s_2_0 = _mm256_sub_ps(s_2_0, _mm256_set_ps(sum_base[5], sum_base[4], sum_base[5], sum_base[4], sum_base[5], sum_base[4], 0, 0));
        _mm256_store_ps((float*)tmp, s_2_0);
        sum[2] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float* y_base = (float*) y;
        __m256 v_0, v_1, v_2, v_3;
        __m256 y_0, y_1;
        __m256 q_0, q_1, q_2, q_3;
        __m256 s_0, s_1, s_2, s_3;
        __m256 s_buffer[(MAX_FOLD * 4)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 3)] = (__m256)_mm256_broadcast_sd((double *)(sum_base + (j * 2)));
        }
        if(incv == 1 && incy == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16, y_base += 16){
            v_0 = _mm256_loadu_ps(v_base);
            v_1 = _mm256_loadu_ps(v_base + 8);
            y_0 = _mm256_loadu_ps(y_base);
            y_1 = _mm256_loadu_ps(y_base + 8);
            v_2 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_1, 0b10110001), _mm256_permute_ps(y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            v_1 = _mm256_mul_ps(v_1, _mm256_permute_ps(y_1, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              q_2 = _mm256_add_ps(s_2, _mm256_or_ps(v_2, mask_BLP));
              q_3 = _mm256_add_ps(s_3, _mm256_or_ps(v_3, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              q_2 = _mm256_sub_ps(s_2, q_2);
              q_3 = _mm256_sub_ps(s_3, q_3);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
              v_2 = _mm256_add_ps(v_2, q_2);
              v_3 = _mm256_add_ps(v_3, q_3);
            }
            s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm256_add_ps(s_buffer[((j * 4) + 1)], _mm256_or_ps(v_1, mask_BLP));
            s_buffer[((j * 4) + 2)] = _mm256_add_ps(s_buffer[((j * 4) + 2)], _mm256_or_ps(v_2, mask_BLP));
            s_buffer[((j * 4) + 3)] = _mm256_add_ps(s_buffer[((j * 4) + 3)], _mm256_or_ps(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm256_loadu_ps(v_base);
            y_0 = _mm256_loadu_ps(y_base);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
            }
            s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm256_add_ps(s_buffer[((j * 4) + 1)], _mm256_or_ps(v_1, mask_BLP));
            i += 4, v_base += 8, y_base += 8;
          }
          if(i < n){
            v_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[2]:0, (n - i)>1?((double*)v_base)[1]:0, ((double*)v_base)[0]);
            y_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)y_base)[2]:0, (n - i)>1?((double*)y_base)[1]:0, ((double*)y_base)[0]);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
            }
            s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm256_add_ps(s_buffer[((j * 4) + 1)], _mm256_or_ps(v_1, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16), y_base += (incy * 16)){
            v_0 = _mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            v_1 = _mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]);
            y_0 = _mm256_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)], y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            y_1 = _mm256_set_ps(y_base[((incy * 14) + 1)], y_base[(incy * 14)], y_base[((incy * 12) + 1)], y_base[(incy * 12)], y_base[((incy * 10) + 1)], y_base[(incy * 10)], y_base[((incy * 8) + 1)], y_base[(incy * 8)]);
            v_2 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_1, 0b10110001), _mm256_permute_ps(y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            v_1 = _mm256_mul_ps(v_1, _mm256_permute_ps(y_1, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              q_2 = _mm256_add_ps(s_2, _mm256_or_ps(v_2, mask_BLP));
              q_3 = _mm256_add_ps(s_3, _mm256_or_ps(v_3, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              q_2 = _mm256_sub_ps(s_2, q_2);
              q_3 = _mm256_sub_ps(s_3, q_3);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
              v_2 = _mm256_add_ps(v_2, q_2);
              v_3 = _mm256_add_ps(v_3, q_3);
            }
            s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm256_add_ps(s_buffer[((j * 4) + 1)], _mm256_or_ps(v_1, mask_BLP));
            s_buffer[((j * 4) + 2)] = _mm256_add_ps(s_buffer[((j * 4) + 2)], _mm256_or_ps(v_2, mask_BLP));
            s_buffer[((j * 4) + 3)] = _mm256_add_ps(s_buffer[((j * 4) + 3)], _mm256_or_ps(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            y_0 = _mm256_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)], y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
            }
            s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm256_add_ps(s_buffer[((j * 4) + 1)], _mm256_or_ps(v_1, mask_BLP));
            i += 4, v_base += (incv * 8), y_base += (incy * 8);
          }
          if(i < n){
            v_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[(incv * 2)]:0, (n - i)>1?((double*)v_base)[incv]:0, ((double*)v_base)[0]);
            y_0 = (__m256)_mm256_set_pd(0, (n - i)>2?((double*)y_base)[(incy * 2)]:0, (n - i)>1?((double*)y_base)[incy]:0, ((double*)y_base)[0]);
            v_1 = _mm256_xor_ps(_mm256_mul_ps(_mm256_permute_ps(v_0, 0b10110001), _mm256_permute_ps(y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm256_mul_ps(v_0, _mm256_permute_ps(y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
            }
            s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm256_add_ps(s_buffer[((j * 4) + 1)], _mm256_or_ps(v_1, mask_BLP));
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = _mm256_sub_ps(s_buffer[(j * 4)], _mm256_set_ps(sum_base[((j * 2) + 1)], sum_base[(j * 2)], sum_base[((j * 2) + 1)], sum_base[(j * 2)], sum_base[((j * 2) + 1)], sum_base[(j * 2)], 0, 0));
          q_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base + (j * 2)));
          s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_sub_ps(s_buffer[((j * 4) + 1)], q_0));
          s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_sub_ps(s_buffer[((j * 4) + 2)], q_0));
          s_buffer[(j * 4)] = _mm256_add_ps(s_buffer[(j * 4)], _mm256_sub_ps(s_buffer[((j * 4) + 3)], q_0));
          _mm256_store_ps((float*)tmp, s_buffer[(j * 4)]);
          sum[j] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#elif defined( __SSE2__ )
  void cdotuI2(int n, float complex* v, int incv, float complex* y, int incy, int fold, float complex* sum){
    __m128 mask_NCONJ; SSE_NCONJ_MASKS(mask_NCONJ);
    __m128 mask_BLP; SSE_BLP_MASKS(mask_BLP);
    float complex tmp[2] __attribute__((aligned(16)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float* y_base = (float*) y;
        __m128 v_0, v_1, v_2, v_3;
        __m128 y_0, y_1;
        __m128 q_0;
        __m128 s_0_0;
        __m128 s_1_0;
        __m128 s_2_0;

        s_0_0 = (__m128)_mm_load1_pd((double *)(sum_base));
        s_1_0 = (__m128)_mm_load1_pd((double *)(sum_base + 2));
        s_2_0 = (__m128)_mm_load1_pd((double *)(sum_base + 4));
        if(incv == 1 && incy == 1){

          for(i = 0; i + 4 <= n; i += 4, v_base += 8, y_base += 8){
            v_0 = _mm_loadu_ps(v_base);
            v_1 = _mm_loadu_ps(v_base + 4);
            y_0 = _mm_loadu_ps(y_base);
            y_1 = _mm_loadu_ps(y_base + 4);
            v_2 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_1, v_1, 0b10110001), _mm_shuffle_ps(y_1, y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            v_1 = _mm_mul_ps(v_1, _mm_shuffle_ps(y_1, y_1, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_1 = _mm_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_1 = _mm_add_ps(v_1, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_2, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_2 = _mm_add_ps(v_2, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_2, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_2 = _mm_add_ps(v_2, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_2, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_3 = _mm_add_ps(v_3, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_3 = _mm_add_ps(v_3, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_3, mask_BLP));
          }
          if(i + 2 <= n){
            v_0 = _mm_loadu_ps(v_base);
            y_0 = _mm_loadu_ps(y_base);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_1 = _mm_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_1 = _mm_add_ps(v_1, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_1, mask_BLP));
            i += 2, v_base += 4, y_base += 4;
          }
          if(i < n){
            v_0 = _mm_set_ps(0, 0, v_base[1], v_base[0]);
            y_0 = _mm_set_ps(0, 0, y_base[1], y_base[0]);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_1 = _mm_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_1 = _mm_add_ps(v_1, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_1, mask_BLP));
          }
        }else{

          for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8), y_base += (incy * 8)){
            v_0 = _mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            v_1 = _mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]);
            y_0 = _mm_set_ps(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            y_1 = _mm_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)]);
            v_2 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_1, v_1, 0b10110001), _mm_shuffle_ps(y_1, y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            v_1 = _mm_mul_ps(v_1, _mm_shuffle_ps(y_1, y_1, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_1 = _mm_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_1 = _mm_add_ps(v_1, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_2, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_2 = _mm_add_ps(v_2, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_2, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_2 = _mm_add_ps(v_2, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_2, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_3 = _mm_add_ps(v_3, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_3 = _mm_add_ps(v_3, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_3, mask_BLP));
          }
          if(i + 2 <= n){
            v_0 = _mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            y_0 = _mm_set_ps(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_1 = _mm_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_1 = _mm_add_ps(v_1, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_1, mask_BLP));
            i += 2, v_base += (incv * 4), y_base += (incy * 4);
          }
          if(i < n){
            v_0 = _mm_set_ps(0, 0, v_base[1], v_base[0]);
            y_0 = _mm_set_ps(0, 0, y_base[1], y_base[0]);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_1 = _mm_add_ps(v_1, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_1 = _mm_add_ps(v_1, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_1, mask_BLP));
          }
        }
        s_0_0 = _mm_sub_ps(s_0_0, _mm_set_ps(sum_base[1], sum_base[0], 0, 0));
        _mm_store_ps((float*)tmp, s_0_0);
        sum[0] = tmp[0] + tmp[1];
        s_1_0 = _mm_sub_ps(s_1_0, _mm_set_ps(sum_base[3], sum_base[2], 0, 0));
        _mm_store_ps((float*)tmp, s_1_0);
        sum[1] = tmp[0] + tmp[1];
        s_2_0 = _mm_sub_ps(s_2_0, _mm_set_ps(sum_base[5], sum_base[4], 0, 0));
        _mm_store_ps((float*)tmp, s_2_0);
        sum[2] = tmp[0] + tmp[1];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float* y_base = (float*) y;
        __m128 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        __m128 y_0, y_1, y_2, y_3;
        __m128 q_0, q_1, q_2, q_3, q_4, q_5, q_6, q_7;
        __m128 s_0, s_1, s_2, s_3, s_4, s_5, s_6, s_7;
        __m128 s_buffer[(MAX_FOLD * 8)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 8)] = s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 2)] = s_buffer[((j * 8) + 3)] = s_buffer[((j * 8) + 4)] = s_buffer[((j * 8) + 5)] = s_buffer[((j * 8) + 6)] = s_buffer[((j * 8) + 7)] = (__m128)_mm_load1_pd((double *)(sum_base + (j * 2)));
        }
        if(incv == 1 && incy == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16, y_base += 16){
            v_0 = _mm_loadu_ps(v_base);
            v_1 = _mm_loadu_ps(v_base + 4);
            v_2 = _mm_loadu_ps(v_base + 8);
            v_3 = _mm_loadu_ps(v_base + 12);
            y_0 = _mm_loadu_ps(y_base);
            y_1 = _mm_loadu_ps(y_base + 4);
            y_2 = _mm_loadu_ps(y_base + 8);
            y_3 = _mm_loadu_ps(y_base + 12);
            v_4 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_5 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_1, v_1, 0b10110001), _mm_shuffle_ps(y_1, y_1, 0b11110101)), mask_NCONJ);
            v_6 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_2, v_2, 0b10110001), _mm_shuffle_ps(y_2, y_2, 0b11110101)), mask_NCONJ);
            v_7 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_3, v_3, 0b10110001), _mm_shuffle_ps(y_3, y_3, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            v_1 = _mm_mul_ps(v_1, _mm_shuffle_ps(y_1, y_1, 0b10100000));
            v_2 = _mm_mul_ps(v_2, _mm_shuffle_ps(y_2, y_2, 0b10100000));
            v_3 = _mm_mul_ps(v_3, _mm_shuffle_ps(y_3, y_3, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              s_4 = s_buffer[((j * 8) + 4)];
              s_5 = s_buffer[((j * 8) + 5)];
              s_6 = s_buffer[((j * 8) + 6)];
              s_7 = s_buffer[((j * 8) + 7)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              q_2 = _mm_add_ps(s_2, _mm_or_ps(v_2, mask_BLP));
              q_3 = _mm_add_ps(s_3, _mm_or_ps(v_3, mask_BLP));
              q_4 = _mm_add_ps(s_4, _mm_or_ps(v_4, mask_BLP));
              q_5 = _mm_add_ps(s_5, _mm_or_ps(v_5, mask_BLP));
              q_6 = _mm_add_ps(s_6, _mm_or_ps(v_6, mask_BLP));
              q_7 = _mm_add_ps(s_7, _mm_or_ps(v_7, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              s_buffer[((j * 8) + 4)] = q_4;
              s_buffer[((j * 8) + 5)] = q_5;
              s_buffer[((j * 8) + 6)] = q_6;
              s_buffer[((j * 8) + 7)] = q_7;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              q_2 = _mm_sub_ps(s_2, q_2);
              q_3 = _mm_sub_ps(s_3, q_3);
              q_4 = _mm_sub_ps(s_4, q_4);
              q_5 = _mm_sub_ps(s_5, q_5);
              q_6 = _mm_sub_ps(s_6, q_6);
              q_7 = _mm_sub_ps(s_7, q_7);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
              v_2 = _mm_add_ps(v_2, q_2);
              v_3 = _mm_add_ps(v_3, q_3);
              v_4 = _mm_add_ps(v_4, q_4);
              v_5 = _mm_add_ps(v_5, q_5);
              v_6 = _mm_add_ps(v_6, q_6);
              v_7 = _mm_add_ps(v_7, q_7);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
            s_buffer[((j * 8) + 2)] = _mm_add_ps(s_buffer[((j * 8) + 2)], _mm_or_ps(v_2, mask_BLP));
            s_buffer[((j * 8) + 3)] = _mm_add_ps(s_buffer[((j * 8) + 3)], _mm_or_ps(v_3, mask_BLP));
            s_buffer[((j * 8) + 4)] = _mm_add_ps(s_buffer[((j * 8) + 4)], _mm_or_ps(v_4, mask_BLP));
            s_buffer[((j * 8) + 5)] = _mm_add_ps(s_buffer[((j * 8) + 5)], _mm_or_ps(v_5, mask_BLP));
            s_buffer[((j * 8) + 6)] = _mm_add_ps(s_buffer[((j * 8) + 6)], _mm_or_ps(v_6, mask_BLP));
            s_buffer[((j * 8) + 7)] = _mm_add_ps(s_buffer[((j * 8) + 7)], _mm_or_ps(v_7, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_loadu_ps(v_base);
            v_1 = _mm_loadu_ps(v_base + 4);
            y_0 = _mm_loadu_ps(y_base);
            y_1 = _mm_loadu_ps(y_base + 4);
            v_2 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_1, v_1, 0b10110001), _mm_shuffle_ps(y_1, y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            v_1 = _mm_mul_ps(v_1, _mm_shuffle_ps(y_1, y_1, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              q_2 = _mm_add_ps(s_2, _mm_or_ps(v_2, mask_BLP));
              q_3 = _mm_add_ps(s_3, _mm_or_ps(v_3, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              q_2 = _mm_sub_ps(s_2, q_2);
              q_3 = _mm_sub_ps(s_3, q_3);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
              v_2 = _mm_add_ps(v_2, q_2);
              v_3 = _mm_add_ps(v_3, q_3);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
            s_buffer[((j * 8) + 2)] = _mm_add_ps(s_buffer[((j * 8) + 2)], _mm_or_ps(v_2, mask_BLP));
            s_buffer[((j * 8) + 3)] = _mm_add_ps(s_buffer[((j * 8) + 3)], _mm_or_ps(v_3, mask_BLP));
            i += 4, v_base += 8, y_base += 8;
          }
          if(i + 2 <= n){
            v_0 = _mm_loadu_ps(v_base);
            y_0 = _mm_loadu_ps(y_base);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
            i += 2, v_base += 4, y_base += 4;
          }
          if(i < n){
            v_0 = _mm_set_ps(0, 0, v_base[1], v_base[0]);
            y_0 = _mm_set_ps(0, 0, y_base[1], y_base[0]);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16), y_base += (incy * 16)){
            v_0 = _mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            v_1 = _mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]);
            v_2 = _mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]);
            v_3 = _mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]);
            y_0 = _mm_set_ps(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            y_1 = _mm_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)]);
            y_2 = _mm_set_ps(y_base[((incy * 10) + 1)], y_base[(incy * 10)], y_base[((incy * 8) + 1)], y_base[(incy * 8)]);
            y_3 = _mm_set_ps(y_base[((incy * 14) + 1)], y_base[(incy * 14)], y_base[((incy * 12) + 1)], y_base[(incy * 12)]);
            v_4 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_5 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_1, v_1, 0b10110001), _mm_shuffle_ps(y_1, y_1, 0b11110101)), mask_NCONJ);
            v_6 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_2, v_2, 0b10110001), _mm_shuffle_ps(y_2, y_2, 0b11110101)), mask_NCONJ);
            v_7 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_3, v_3, 0b10110001), _mm_shuffle_ps(y_3, y_3, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            v_1 = _mm_mul_ps(v_1, _mm_shuffle_ps(y_1, y_1, 0b10100000));
            v_2 = _mm_mul_ps(v_2, _mm_shuffle_ps(y_2, y_2, 0b10100000));
            v_3 = _mm_mul_ps(v_3, _mm_shuffle_ps(y_3, y_3, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              s_4 = s_buffer[((j * 8) + 4)];
              s_5 = s_buffer[((j * 8) + 5)];
              s_6 = s_buffer[((j * 8) + 6)];
              s_7 = s_buffer[((j * 8) + 7)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              q_2 = _mm_add_ps(s_2, _mm_or_ps(v_2, mask_BLP));
              q_3 = _mm_add_ps(s_3, _mm_or_ps(v_3, mask_BLP));
              q_4 = _mm_add_ps(s_4, _mm_or_ps(v_4, mask_BLP));
              q_5 = _mm_add_ps(s_5, _mm_or_ps(v_5, mask_BLP));
              q_6 = _mm_add_ps(s_6, _mm_or_ps(v_6, mask_BLP));
              q_7 = _mm_add_ps(s_7, _mm_or_ps(v_7, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              s_buffer[((j * 8) + 4)] = q_4;
              s_buffer[((j * 8) + 5)] = q_5;
              s_buffer[((j * 8) + 6)] = q_6;
              s_buffer[((j * 8) + 7)] = q_7;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              q_2 = _mm_sub_ps(s_2, q_2);
              q_3 = _mm_sub_ps(s_3, q_3);
              q_4 = _mm_sub_ps(s_4, q_4);
              q_5 = _mm_sub_ps(s_5, q_5);
              q_6 = _mm_sub_ps(s_6, q_6);
              q_7 = _mm_sub_ps(s_7, q_7);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
              v_2 = _mm_add_ps(v_2, q_2);
              v_3 = _mm_add_ps(v_3, q_3);
              v_4 = _mm_add_ps(v_4, q_4);
              v_5 = _mm_add_ps(v_5, q_5);
              v_6 = _mm_add_ps(v_6, q_6);
              v_7 = _mm_add_ps(v_7, q_7);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
            s_buffer[((j * 8) + 2)] = _mm_add_ps(s_buffer[((j * 8) + 2)], _mm_or_ps(v_2, mask_BLP));
            s_buffer[((j * 8) + 3)] = _mm_add_ps(s_buffer[((j * 8) + 3)], _mm_or_ps(v_3, mask_BLP));
            s_buffer[((j * 8) + 4)] = _mm_add_ps(s_buffer[((j * 8) + 4)], _mm_or_ps(v_4, mask_BLP));
            s_buffer[((j * 8) + 5)] = _mm_add_ps(s_buffer[((j * 8) + 5)], _mm_or_ps(v_5, mask_BLP));
            s_buffer[((j * 8) + 6)] = _mm_add_ps(s_buffer[((j * 8) + 6)], _mm_or_ps(v_6, mask_BLP));
            s_buffer[((j * 8) + 7)] = _mm_add_ps(s_buffer[((j * 8) + 7)], _mm_or_ps(v_7, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            v_1 = _mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]);
            y_0 = _mm_set_ps(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            y_1 = _mm_set_ps(y_base[((incy * 6) + 1)], y_base[(incy * 6)], y_base[((incy * 4) + 1)], y_base[(incy * 4)]);
            v_2 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_3 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_1, v_1, 0b10110001), _mm_shuffle_ps(y_1, y_1, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            v_1 = _mm_mul_ps(v_1, _mm_shuffle_ps(y_1, y_1, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              q_2 = _mm_add_ps(s_2, _mm_or_ps(v_2, mask_BLP));
              q_3 = _mm_add_ps(s_3, _mm_or_ps(v_3, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              q_2 = _mm_sub_ps(s_2, q_2);
              q_3 = _mm_sub_ps(s_3, q_3);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
              v_2 = _mm_add_ps(v_2, q_2);
              v_3 = _mm_add_ps(v_3, q_3);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
            s_buffer[((j * 8) + 2)] = _mm_add_ps(s_buffer[((j * 8) + 2)], _mm_or_ps(v_2, mask_BLP));
            s_buffer[((j * 8) + 3)] = _mm_add_ps(s_buffer[((j * 8) + 3)], _mm_or_ps(v_3, mask_BLP));
            i += 4, v_base += (incv * 8), y_base += (incy * 8);
          }
          if(i + 2 <= n){
            v_0 = _mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]);
            y_0 = _mm_set_ps(y_base[((incy * 2) + 1)], y_base[(incy * 2)], y_base[1], y_base[0]);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
            i += 2, v_base += (incv * 4), y_base += (incy * 4);
          }
          if(i < n){
            v_0 = _mm_set_ps(0, 0, v_base[1], v_base[0]);
            y_0 = _mm_set_ps(0, 0, y_base[1], y_base[0]);
            v_1 = _mm_xor_ps(_mm_mul_ps(_mm_shuffle_ps(v_0, v_0, 0b10110001), _mm_shuffle_ps(y_0, y_0, 0b11110101)), mask_NCONJ);
            v_0 = _mm_mul_ps(v_0, _mm_shuffle_ps(y_0, y_0, 0b10100000));
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 8) + 1)] = _mm_add_ps(s_buffer[((j * 8) + 1)], _mm_or_ps(v_1, mask_BLP));
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 8)] = _mm_sub_ps(s_buffer[(j * 8)], _mm_set_ps(sum_base[((j * 2) + 1)], sum_base[(j * 2)], 0, 0));
          q_0 = (__m128)_mm_load1_pd((double *)(sum_base + (j * 2)));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 1)], q_0));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 2)], q_0));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 3)], q_0));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 4)], q_0));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 5)], q_0));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 6)], q_0));
          s_buffer[(j * 8)] = _mm_add_ps(s_buffer[(j * 8)], _mm_sub_ps(s_buffer[((j * 8) + 7)], q_0));
          _mm_store_ps((float*)tmp, s_buffer[(j * 8)]);
          sum[j] = tmp[0] + tmp[1];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#else
  void cdotuI2(int n, float complex* v, int incv, float complex* y, int incy, int fold, float complex* sum){
    i_float tmp_BLP;
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float* y_base = (float*) y;
        float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        float y_0, y_1, y_2, y_3;
        float q_0, q_1;
        float s_0_0, s_0_1;
        float s_1_0, s_1_1;
        float s_2_0, s_2_1;

        s_0_0 = sum_base[0];
        s_0_1 = sum_base[1];
        s_1_0 = sum_base[2];
        s_1_1 = sum_base[3];
        s_2_0 = sum_base[4];
        s_2_1 = sum_base[5];
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
            v_4 = v_1 * y_1 * -1;
            v_5 = v_0 * y_1;
            v_6 = v_3 * y_3 * -1;
            v_7 = v_2 * y_3;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            v_2 = v_2 * y_2;
            v_3 = v_3 * y_2;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_6 = v_6 + q_0;
            v_7 = v_7 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_6 = v_6 + q_0;
            v_7 = v_7 + q_1;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = v_base[0];
            v_1 = v_base[1];
            y_0 = y_base[0];
            y_1 = y_base[1];
            v_2 = v_1 * y_1 * -1;
            v_3 = v_0 * y_1;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
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
            v_4 = v_1 * y_1 * -1;
            v_5 = v_0 * y_1;
            v_6 = v_3 * y_3 * -1;
            v_7 = v_2 * y_3;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            v_2 = v_2 * y_2;
            v_3 = v_3 * y_2;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_6 = v_6 + q_0;
            v_7 = v_7 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_6 = v_6 + q_0;
            v_7 = v_7 + q_1;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = v_base[0];
            v_1 = v_base[1];
            y_0 = y_base[0];
            y_1 = y_base[1];
            v_2 = v_1 * y_1 * -1;
            v_3 = v_0 * y_1;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            q_0 = s_0_0;
            q_1 = s_0_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            q_0 = s_1_0;
            q_1 = s_1_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            v_2 = v_2 + q_0;
            v_3 = v_3 + q_1;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            i += 1, v_base += (incv * 2), y_base += (incy * 2);
          }
        }
        ((float*)sum)[0] = s_0_0;
        ((float*)sum)[1] = s_0_1;
        ((float*)sum)[2] = s_1_0;
        ((float*)sum)[3] = s_1_1;
        ((float*)sum)[4] = s_2_0;
        ((float*)sum)[5] = s_2_1;
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float* y_base = (float*) y;
        float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        float y_0, y_1, y_2, y_3;
        float q_0, q_1, q_2, q_3, q_4, q_5, q_6, q_7;
        float s_0, s_1, s_2, s_3, s_4, s_5, s_6, s_7;
        float s_buffer[(MAX_FOLD * 8)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 8)] = s_buffer[((j * 8) + 2)] = s_buffer[((j * 8) + 4)] = s_buffer[((j * 8) + 6)] = sum_base[(j * 2)];
          s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 3)] = s_buffer[((j * 8) + 5)] = s_buffer[((j * 8) + 7)] = sum_base[((j * 2) + 1)];
        }
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
            v_4 = v_1 * y_1 * -1;
            v_5 = v_0 * y_1;
            v_6 = v_3 * y_3 * -1;
            v_7 = v_2 * y_3;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            v_2 = v_2 * y_2;
            v_3 = v_3 * y_2;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              s_4 = s_buffer[((j * 8) + 4)];
              s_5 = s_buffer[((j * 8) + 5)];
              s_6 = s_buffer[((j * 8) + 6)];
              s_7 = s_buffer[((j * 8) + 7)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              tmp_BLP.f = v_2;
              tmp_BLP.i |= 1;
              q_2 = s_2 + tmp_BLP.f;
              tmp_BLP.f = v_3;
              tmp_BLP.i |= 1;
              q_3 = s_3 + tmp_BLP.f;
              tmp_BLP.f = v_4;
              tmp_BLP.i |= 1;
              q_4 = s_4 + tmp_BLP.f;
              tmp_BLP.f = v_5;
              tmp_BLP.i |= 1;
              q_5 = s_5 + tmp_BLP.f;
              tmp_BLP.f = v_6;
              tmp_BLP.i |= 1;
              q_6 = s_6 + tmp_BLP.f;
              tmp_BLP.f = v_7;
              tmp_BLP.i |= 1;
              q_7 = s_7 + tmp_BLP.f;
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              s_buffer[((j * 8) + 4)] = q_4;
              s_buffer[((j * 8) + 5)] = q_5;
              s_buffer[((j * 8) + 6)] = q_6;
              s_buffer[((j * 8) + 7)] = q_7;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              q_2 = s_2 - q_2;
              q_3 = s_3 - q_3;
              q_4 = s_4 - q_4;
              q_5 = s_5 - q_5;
              q_6 = s_6 - q_6;
              q_7 = s_7 - q_7;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              v_2 = v_2 + q_2;
              v_3 = v_3 + q_3;
              v_4 = v_4 + q_4;
              v_5 = v_5 + q_5;
              v_6 = v_6 + q_6;
              v_7 = v_7 + q_7;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 8)] = s_buffer[(j * 8)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 2)] = s_buffer[((j * 8) + 2)] + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 3)] = s_buffer[((j * 8) + 3)] + tmp_BLP.f;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 4)] = s_buffer[((j * 8) + 4)] + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 5)] = s_buffer[((j * 8) + 5)] + tmp_BLP.f;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 6)] = s_buffer[((j * 8) + 6)] + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 7)] = s_buffer[((j * 8) + 7)] + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = v_base[0];
            v_1 = v_base[1];
            y_0 = y_base[0];
            y_1 = y_base[1];
            v_2 = v_1 * y_1 * -1;
            v_3 = v_0 * y_1;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              tmp_BLP.f = v_2;
              tmp_BLP.i |= 1;
              q_2 = s_2 + tmp_BLP.f;
              tmp_BLP.f = v_3;
              tmp_BLP.i |= 1;
              q_3 = s_3 + tmp_BLP.f;
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              q_2 = s_2 - q_2;
              q_3 = s_3 - q_3;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              v_2 = v_2 + q_2;
              v_3 = v_3 + q_3;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 8)] = s_buffer[(j * 8)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 2)] = s_buffer[((j * 8) + 2)] + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 3)] = s_buffer[((j * 8) + 3)] + tmp_BLP.f;
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
            v_4 = v_1 * y_1 * -1;
            v_5 = v_0 * y_1;
            v_6 = v_3 * y_3 * -1;
            v_7 = v_2 * y_3;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            v_2 = v_2 * y_2;
            v_3 = v_3 * y_2;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              s_4 = s_buffer[((j * 8) + 4)];
              s_5 = s_buffer[((j * 8) + 5)];
              s_6 = s_buffer[((j * 8) + 6)];
              s_7 = s_buffer[((j * 8) + 7)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              tmp_BLP.f = v_2;
              tmp_BLP.i |= 1;
              q_2 = s_2 + tmp_BLP.f;
              tmp_BLP.f = v_3;
              tmp_BLP.i |= 1;
              q_3 = s_3 + tmp_BLP.f;
              tmp_BLP.f = v_4;
              tmp_BLP.i |= 1;
              q_4 = s_4 + tmp_BLP.f;
              tmp_BLP.f = v_5;
              tmp_BLP.i |= 1;
              q_5 = s_5 + tmp_BLP.f;
              tmp_BLP.f = v_6;
              tmp_BLP.i |= 1;
              q_6 = s_6 + tmp_BLP.f;
              tmp_BLP.f = v_7;
              tmp_BLP.i |= 1;
              q_7 = s_7 + tmp_BLP.f;
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              s_buffer[((j * 8) + 4)] = q_4;
              s_buffer[((j * 8) + 5)] = q_5;
              s_buffer[((j * 8) + 6)] = q_6;
              s_buffer[((j * 8) + 7)] = q_7;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              q_2 = s_2 - q_2;
              q_3 = s_3 - q_3;
              q_4 = s_4 - q_4;
              q_5 = s_5 - q_5;
              q_6 = s_6 - q_6;
              q_7 = s_7 - q_7;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              v_2 = v_2 + q_2;
              v_3 = v_3 + q_3;
              v_4 = v_4 + q_4;
              v_5 = v_5 + q_5;
              v_6 = v_6 + q_6;
              v_7 = v_7 + q_7;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 8)] = s_buffer[(j * 8)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 2)] = s_buffer[((j * 8) + 2)] + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 3)] = s_buffer[((j * 8) + 3)] + tmp_BLP.f;
            tmp_BLP.f = v_4;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 4)] = s_buffer[((j * 8) + 4)] + tmp_BLP.f;
            tmp_BLP.f = v_5;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 5)] = s_buffer[((j * 8) + 5)] + tmp_BLP.f;
            tmp_BLP.f = v_6;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 6)] = s_buffer[((j * 8) + 6)] + tmp_BLP.f;
            tmp_BLP.f = v_7;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 7)] = s_buffer[((j * 8) + 7)] + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = v_base[0];
            v_1 = v_base[1];
            y_0 = y_base[0];
            y_1 = y_base[1];
            v_2 = v_1 * y_1 * -1;
            v_3 = v_0 * y_1;
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_0;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 8)];
              s_1 = s_buffer[((j * 8) + 1)];
              s_2 = s_buffer[((j * 8) + 2)];
              s_3 = s_buffer[((j * 8) + 3)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              tmp_BLP.f = v_2;
              tmp_BLP.i |= 1;
              q_2 = s_2 + tmp_BLP.f;
              tmp_BLP.f = v_3;
              tmp_BLP.i |= 1;
              q_3 = s_3 + tmp_BLP.f;
              s_buffer[(j * 8)] = q_0;
              s_buffer[((j * 8) + 1)] = q_1;
              s_buffer[((j * 8) + 2)] = q_2;
              s_buffer[((j * 8) + 3)] = q_3;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              q_2 = s_2 - q_2;
              q_3 = s_3 - q_3;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              v_2 = v_2 + q_2;
              v_3 = v_3 + q_3;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 8)] = s_buffer[(j * 8)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 2)] = s_buffer[((j * 8) + 2)] + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_buffer[((j * 8) + 3)] = s_buffer[((j * 8) + 3)] + tmp_BLP.f;
            i += 1, v_base += (incv * 2), y_base += (incy * 2);
          }
        }
        for(j = 0; j < fold; j += 1){
          q_0 = ((float*)sum)[(j * 2)];
          s_buffer[(j * 8)] = s_buffer[(j * 8)] + (s_buffer[((j * 8) + 2)] - q_0);
          s_buffer[(j * 8)] = s_buffer[(j * 8)] + (s_buffer[((j * 8) + 4)] - q_0);
          s_buffer[(j * 8)] = s_buffer[(j * 8)] + (s_buffer[((j * 8) + 6)] - q_0);
          q_0 = ((float*)sum)[((j * 2) + 1)];
          s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + (s_buffer[((j * 8) + 3)] - q_0);
          s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + (s_buffer[((j * 8) + 5)] - q_0);
          s_buffer[((j * 8) + 1)] = s_buffer[((j * 8) + 1)] + (s_buffer[((j * 8) + 7)] - q_0);
          ((float*)sum)[(j * 2)] = s_buffer[(j * 8)];
          ((float*)sum)[((j * 2) + 1)] = s_buffer[((j * 8) + 1)];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#endif