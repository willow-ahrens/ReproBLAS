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
  void scasumI2(int n, float complex* v, int incv, int fold, float complex* sum){
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    __m256 mask_BLP; AVX_BLP_MASKS(mask_BLP);
    float complex tmp_cons[4] __attribute__((aligned(32)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        __m256 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        __m256 q_0, q_1;
        __m256 s_0_0, s_0_1;
        __m256 s_1_0, s_1_1;
        __m256 s_2_0, s_2_1;

        s_0_0 = s_0_1 = (__m256)_mm256_broadcast_sd((double *)(sum_base));
        s_1_0 = s_1_1 = (__m256)_mm256_broadcast_sd((double *)(sum_base + 2));
        s_2_0 = s_2_1 = (__m256)_mm256_broadcast_sd((double *)(sum_base + 4));
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
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_4, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_4 = _mm256_add_ps(v_4, q_0);
            v_5 = _mm256_add_ps(v_5, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_4, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_4 = _mm256_add_ps(v_4, q_0);
            v_5 = _mm256_add_ps(v_5, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_4, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_5, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_6, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_6 = _mm256_add_ps(v_6, q_0);
            v_7 = _mm256_add_ps(v_7, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_6, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_6 = _mm256_add_ps(v_6, q_0);
            v_7 = _mm256_add_ps(v_7, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_6, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_7, mask_BLP));
          }
          if(i + 16 <= n){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
            v_2 = _mm256_and_ps(_mm256_loadu_ps(v_base + 16), mask_ABS);
            v_3 = _mm256_and_ps(_mm256_loadu_ps(v_base + 24), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_3, mask_BLP));
            i += 16, v_base += 32;
          }
          if(i + 8 <= n){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_1, mask_BLP));
            i += 8, v_base += 16;
          }
          if(i + 4 <= n){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            i += 4, v_base += 8;
          }
          if(i < n){
            v_0 = _mm256_and_ps((__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[2]:0, (n - i)>1?((double*)v_base)[1]:0, ((double*)v_base)[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
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
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_4, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_4 = _mm256_add_ps(v_4, q_0);
            v_5 = _mm256_add_ps(v_5, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_4, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_5, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_4 = _mm256_add_ps(v_4, q_0);
            v_5 = _mm256_add_ps(v_5, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_4, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_5, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_6, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_6 = _mm256_add_ps(v_6, q_0);
            v_7 = _mm256_add_ps(v_7, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_6, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_7, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_6 = _mm256_add_ps(v_6, q_0);
            v_7 = _mm256_add_ps(v_7, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_6, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_7, mask_BLP));
          }
          if(i + 16 <= n){
            v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
            v_2 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 22) + 1)], v_base[(incv * 22)], v_base[((incv * 20) + 1)], v_base[(incv * 20)], v_base[((incv * 18) + 1)], v_base[(incv * 18)], v_base[((incv * 16) + 1)], v_base[(incv * 16)]), mask_ABS);
            v_3 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 30) + 1)], v_base[(incv * 30)], v_base[((incv * 28) + 1)], v_base[(incv * 28)], v_base[((incv * 26) + 1)], v_base[(incv * 26)], v_base[((incv * 24) + 1)], v_base[(incv * 24)]), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_2, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_2, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_3, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_2 = _mm256_add_ps(v_2, q_0);
            v_3 = _mm256_add_ps(v_3, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_2, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_3, mask_BLP));
            i += 16, v_base += (incv * 32);
          }
          if(i + 8 <= n){
            v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            s_0_1 = _mm256_add_ps(s_0_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            q_1 = _mm256_sub_ps(q_1, s_0_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            s_1_1 = _mm256_add_ps(s_1_1, _mm256_or_ps(v_1, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            q_1 = _mm256_sub_ps(q_1, s_1_1);
            v_0 = _mm256_add_ps(v_0, q_0);
            v_1 = _mm256_add_ps(v_1, q_1);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            s_2_1 = _mm256_add_ps(s_2_1, _mm256_or_ps(v_1, mask_BLP));
            i += 8, v_base += (incv * 16);
          }
          if(i + 4 <= n){
            v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            i += 4, v_base += (incv * 8);
          }
          if(i < n){
            v_0 = _mm256_and_ps((__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[(incv * 2)]:0, (n - i)>1?((double*)v_base)[incv]:0, ((double*)v_base)[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
          }
        }
        s_0_0 = _mm256_sub_ps(s_0_0, _mm256_set_ps(sum_base[1], sum_base[0], sum_base[1], sum_base[0], sum_base[1], sum_base[0], 0, 0));
        q_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base));
        s_0_0 = _mm256_add_ps(s_0_0, _mm256_sub_ps(s_0_1, q_0));
        _mm256_store_ps((float*)tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        s_1_0 = _mm256_sub_ps(s_1_0, _mm256_set_ps(sum_base[3], sum_base[2], sum_base[3], sum_base[2], sum_base[3], sum_base[2], 0, 0));
        q_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base + 2));
        s_1_0 = _mm256_add_ps(s_1_0, _mm256_sub_ps(s_1_1, q_0));
        _mm256_store_ps((float*)tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        s_2_0 = _mm256_sub_ps(s_2_0, _mm256_set_ps(sum_base[5], sum_base[4], sum_base[5], sum_base[4], sum_base[5], sum_base[4], 0, 0));
        q_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base + 4));
        s_2_0 = _mm256_add_ps(s_2_0, _mm256_sub_ps(s_2_1, q_0));
        _mm256_store_ps((float*)tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        __m256 v_0, v_1;
        __m256 q_0, q_1;
        __m256 s_0, s_1;
        __m256 s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = (__m256)_mm256_broadcast_sd((double *)(sum_base + (j * 2)));
        }
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_loadu_ps(v_base + 8), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
            }
            s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm256_add_ps(s_buffer[((j * 2) + 1)], _mm256_or_ps(v_1, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v_base), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_or_ps(v_0, mask_BLP));
            i += 4, v_base += 8;
          }
          if(i < n){
            v_0 = _mm256_and_ps((__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[2]:0, (n - i)>1?((double*)v_base)[1]:0, ((double*)v_base)[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_or_ps(v_0, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
            v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)], v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              q_1 = _mm256_add_ps(s_1, _mm256_or_ps(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm256_sub_ps(s_0, q_0);
              q_1 = _mm256_sub_ps(s_1, q_1);
              v_0 = _mm256_add_ps(v_0, q_0);
              v_1 = _mm256_add_ps(v_1, q_1);
            }
            s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_or_ps(v_0, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm256_add_ps(s_buffer[((j * 2) + 1)], _mm256_or_ps(v_1, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm256_and_ps(_mm256_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)], v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_or_ps(v_0, mask_BLP));
            i += 4, v_base += (incv * 8);
          }
          if(i < n){
            v_0 = _mm256_and_ps((__m256)_mm256_set_pd(0, (n - i)>2?((double*)v_base)[(incv * 2)]:0, (n - i)>1?((double*)v_base)[incv]:0, ((double*)v_base)[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_or_ps(v_0, mask_BLP));
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = _mm256_sub_ps(s_buffer[(j * 2)], _mm256_set_ps(sum_base[((j * 2) + 1)], sum_base[(j * 2)], sum_base[((j * 2) + 1)], sum_base[(j * 2)], sum_base[((j * 2) + 1)], sum_base[(j * 2)], 0, 0));
          q_0 = (__m256)_mm256_broadcast_sd((double *)(sum_base + (j * 2)));
          s_buffer[(j * 2)] = _mm256_add_ps(s_buffer[(j * 2)], _mm256_sub_ps(s_buffer[((j * 2) + 1)], q_0));
          _mm256_store_ps((float*)tmp_cons, s_buffer[(j * 2)]);
          sum[j] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#elif defined( __SSE2__ )
  void scasumI2(int n, float complex* v, int incv, int fold, float complex* sum){
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    __m128 mask_BLP; SSE_BLP_MASKS(mask_BLP);
    float complex tmp_cons[2] __attribute__((aligned(16)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        __m128 v_0, v_1, v_2, v_3;
        __m128 q_0, q_1, q_2, q_3;
        __m128 s_0_0, s_0_1, s_0_2, s_0_3;
        __m128 s_1_0, s_1_1, s_1_2, s_1_3;
        __m128 s_2_0, s_2_1, s_2_2, s_2_3;

        s_0_0 = s_0_1 = s_0_2 = s_0_3 = (__m128)_mm_load1_pd((double *)(sum_base));
        s_1_0 = s_1_1 = s_1_2 = s_1_3 = (__m128)_mm_load1_pd((double *)(sum_base + 2));
        s_2_0 = s_2_1 = s_2_2 = s_2_3 = (__m128)_mm_load1_pd((double *)(sum_base + 4));
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16){
            v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
            v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
            v_2 = _mm_and_ps(_mm_loadu_ps(v_base + 8), mask_ABS);
            v_3 = _mm_and_ps(_mm_loadu_ps(v_base + 12), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            s_0_1 = _mm_add_ps(s_0_1, _mm_or_ps(v_1, mask_BLP));
            s_0_2 = _mm_add_ps(s_0_2, _mm_or_ps(v_2, mask_BLP));
            s_0_3 = _mm_add_ps(s_0_3, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            q_1 = _mm_sub_ps(q_1, s_0_1);
            q_2 = _mm_sub_ps(q_2, s_0_2);
            q_3 = _mm_sub_ps(q_3, s_0_3);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            v_2 = _mm_add_ps(v_2, q_2);
            v_3 = _mm_add_ps(v_3, q_3);
            q_0 = s_1_0;
            q_1 = s_1_1;
            q_2 = s_1_2;
            q_3 = s_1_3;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            s_1_1 = _mm_add_ps(s_1_1, _mm_or_ps(v_1, mask_BLP));
            s_1_2 = _mm_add_ps(s_1_2, _mm_or_ps(v_2, mask_BLP));
            s_1_3 = _mm_add_ps(s_1_3, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            q_1 = _mm_sub_ps(q_1, s_1_1);
            q_2 = _mm_sub_ps(q_2, s_1_2);
            q_3 = _mm_sub_ps(q_3, s_1_3);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            v_2 = _mm_add_ps(v_2, q_2);
            v_3 = _mm_add_ps(v_3, q_3);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            s_2_1 = _mm_add_ps(s_2_1, _mm_or_ps(v_1, mask_BLP));
            s_2_2 = _mm_add_ps(s_2_2, _mm_or_ps(v_2, mask_BLP));
            s_2_3 = _mm_add_ps(s_2_3, _mm_or_ps(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
            v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            s_0_1 = _mm_add_ps(s_0_1, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            q_1 = _mm_sub_ps(q_1, s_0_1);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            s_1_1 = _mm_add_ps(s_1_1, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            q_1 = _mm_sub_ps(q_1, s_1_1);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            s_2_1 = _mm_add_ps(s_2_1, _mm_or_ps(v_1, mask_BLP));
            i += 4, v_base += 8;
          }
          if(i + 2 <= n){
            v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            i += 2, v_base += 4;
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, 0, v_base[1], v_base[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
            v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
            v_2 = _mm_and_ps(_mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
            v_3 = _mm_and_ps(_mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            s_0_1 = _mm_add_ps(s_0_1, _mm_or_ps(v_1, mask_BLP));
            s_0_2 = _mm_add_ps(s_0_2, _mm_or_ps(v_2, mask_BLP));
            s_0_3 = _mm_add_ps(s_0_3, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            q_1 = _mm_sub_ps(q_1, s_0_1);
            q_2 = _mm_sub_ps(q_2, s_0_2);
            q_3 = _mm_sub_ps(q_3, s_0_3);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            v_2 = _mm_add_ps(v_2, q_2);
            v_3 = _mm_add_ps(v_3, q_3);
            q_0 = s_1_0;
            q_1 = s_1_1;
            q_2 = s_1_2;
            q_3 = s_1_3;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            s_1_1 = _mm_add_ps(s_1_1, _mm_or_ps(v_1, mask_BLP));
            s_1_2 = _mm_add_ps(s_1_2, _mm_or_ps(v_2, mask_BLP));
            s_1_3 = _mm_add_ps(s_1_3, _mm_or_ps(v_3, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            q_1 = _mm_sub_ps(q_1, s_1_1);
            q_2 = _mm_sub_ps(q_2, s_1_2);
            q_3 = _mm_sub_ps(q_3, s_1_3);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            v_2 = _mm_add_ps(v_2, q_2);
            v_3 = _mm_add_ps(v_3, q_3);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            s_2_1 = _mm_add_ps(s_2_1, _mm_or_ps(v_1, mask_BLP));
            s_2_2 = _mm_add_ps(s_2_2, _mm_or_ps(v_2, mask_BLP));
            s_2_3 = _mm_add_ps(s_2_3, _mm_or_ps(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            s_0_1 = _mm_add_ps(s_0_1, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            q_1 = _mm_sub_ps(q_1, s_0_1);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            s_1_1 = _mm_add_ps(s_1_1, _mm_or_ps(v_1, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            q_1 = _mm_sub_ps(q_1, s_1_1);
            v_0 = _mm_add_ps(v_0, q_0);
            v_1 = _mm_add_ps(v_1, q_1);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            s_2_1 = _mm_add_ps(s_2_1, _mm_or_ps(v_1, mask_BLP));
            i += 4, v_base += (incv * 8);
          }
          if(i + 2 <= n){
            v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            i += 2, v_base += (incv * 4);
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, 0, v_base[1], v_base[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
          }
        }
        s_0_0 = _mm_sub_ps(s_0_0, _mm_set_ps(sum_base[1], sum_base[0], 0, 0));
        q_0 = (__m128)_mm_load1_pd((double *)(sum_base));
        s_0_0 = _mm_add_ps(s_0_0, _mm_sub_ps(s_0_1, q_0));
        s_0_0 = _mm_add_ps(s_0_0, _mm_sub_ps(s_0_2, q_0));
        s_0_0 = _mm_add_ps(s_0_0, _mm_sub_ps(s_0_3, q_0));
        _mm_store_ps((float*)tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1];
        s_1_0 = _mm_sub_ps(s_1_0, _mm_set_ps(sum_base[3], sum_base[2], 0, 0));
        q_0 = (__m128)_mm_load1_pd((double *)(sum_base + 2));
        s_1_0 = _mm_add_ps(s_1_0, _mm_sub_ps(s_1_1, q_0));
        s_1_0 = _mm_add_ps(s_1_0, _mm_sub_ps(s_1_2, q_0));
        s_1_0 = _mm_add_ps(s_1_0, _mm_sub_ps(s_1_3, q_0));
        _mm_store_ps((float*)tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1];
        s_2_0 = _mm_sub_ps(s_2_0, _mm_set_ps(sum_base[5], sum_base[4], 0, 0));
        q_0 = (__m128)_mm_load1_pd((double *)(sum_base + 4));
        s_2_0 = _mm_add_ps(s_2_0, _mm_sub_ps(s_2_1, q_0));
        s_2_0 = _mm_add_ps(s_2_0, _mm_sub_ps(s_2_2, q_0));
        s_2_0 = _mm_add_ps(s_2_0, _mm_sub_ps(s_2_3, q_0));
        _mm_store_ps((float*)tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        __m128 v_0, v_1, v_2, v_3;
        __m128 q_0, q_1, q_2, q_3;
        __m128 s_0, s_1, s_2, s_3;
        __m128 s_buffer[(MAX_FOLD * 4)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 3)] = (__m128)_mm_load1_pd((double *)(sum_base + (j * 2)));
        }
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16){
            v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
            v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
            v_2 = _mm_and_ps(_mm_loadu_ps(v_base + 8), mask_ABS);
            v_3 = _mm_and_ps(_mm_loadu_ps(v_base + 12), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              q_2 = _mm_add_ps(s_2, _mm_or_ps(v_2, mask_BLP));
              q_3 = _mm_add_ps(s_3, _mm_or_ps(v_3, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              q_2 = _mm_sub_ps(s_2, q_2);
              q_3 = _mm_sub_ps(s_3, q_3);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
              v_2 = _mm_add_ps(v_2, q_2);
              v_3 = _mm_add_ps(v_3, q_3);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm_add_ps(s_buffer[((j * 4) + 1)], _mm_or_ps(v_1, mask_BLP));
            s_buffer[((j * 4) + 2)] = _mm_add_ps(s_buffer[((j * 4) + 2)], _mm_or_ps(v_2, mask_BLP));
            s_buffer[((j * 4) + 3)] = _mm_add_ps(s_buffer[((j * 4) + 3)], _mm_or_ps(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
            v_1 = _mm_and_ps(_mm_loadu_ps(v_base + 4), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm_add_ps(s_buffer[((j * 4) + 1)], _mm_or_ps(v_1, mask_BLP));
            i += 4, v_base += 8;
          }
          if(i + 2 <= n){
            v_0 = _mm_and_ps(_mm_loadu_ps(v_base), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
            i += 2, v_base += 4;
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, 0, v_base[1], v_base[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
            v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
            v_2 = _mm_and_ps(_mm_set_ps(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), mask_ABS);
            v_3 = _mm_and_ps(_mm_set_ps(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              q_2 = _mm_add_ps(s_2, _mm_or_ps(v_2, mask_BLP));
              q_3 = _mm_add_ps(s_3, _mm_or_ps(v_3, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              q_2 = _mm_sub_ps(s_2, q_2);
              q_3 = _mm_sub_ps(s_3, q_3);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
              v_2 = _mm_add_ps(v_2, q_2);
              v_3 = _mm_add_ps(v_3, q_3);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm_add_ps(s_buffer[((j * 4) + 1)], _mm_or_ps(v_1, mask_BLP));
            s_buffer[((j * 4) + 2)] = _mm_add_ps(s_buffer[((j * 4) + 2)], _mm_or_ps(v_2, mask_BLP));
            s_buffer[((j * 4) + 3)] = _mm_add_ps(s_buffer[((j * 4) + 3)], _mm_or_ps(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            v_1 = _mm_and_ps(_mm_set_ps(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 4) + 1)] = _mm_add_ps(s_buffer[((j * 4) + 1)], _mm_or_ps(v_1, mask_BLP));
            i += 4, v_base += (incv * 8);
          }
          if(i + 2 <= n){
            v_0 = _mm_and_ps(_mm_set_ps(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
            i += 2, v_base += (incv * 4);
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, 0, v_base[1], v_base[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 4)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_or_ps(v_0, mask_BLP));
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = _mm_sub_ps(s_buffer[(j * 4)], _mm_set_ps(sum_base[((j * 2) + 1)], sum_base[(j * 2)], 0, 0));
          q_0 = (__m128)_mm_load1_pd((double *)(sum_base + (j * 2)));
          s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_sub_ps(s_buffer[((j * 4) + 1)], q_0));
          s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_sub_ps(s_buffer[((j * 4) + 2)], q_0));
          s_buffer[(j * 4)] = _mm_add_ps(s_buffer[(j * 4)], _mm_sub_ps(s_buffer[((j * 4) + 3)], q_0));
          _mm_store_ps((float*)tmp_cons, s_buffer[(j * 4)]);
          sum[j] = tmp_cons[0] + tmp_cons[1];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#else
  void scasumI2(int n, float complex* v, int incv, int fold, float complex* sum){
    int_float tmp_BLP;
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float v_0, v_1, v_2, v_3;
        float q_0, q_1, q_2, q_3;
        float s_0_0, s_0_1, s_0_2, s_0_3;
        float s_1_0, s_1_1, s_1_2, s_1_3;
        float s_2_0, s_2_1, s_2_2, s_2_3;

        s_0_0 = s_0_2 = sum_base[0];
        s_0_1 = s_0_3 = sum_base[1];
        s_1_0 = s_1_2 = sum_base[2];
        s_1_1 = s_1_3 = sum_base[3];
        s_2_0 = s_2_2 = sum_base[4];
        s_2_1 = s_2_3 = sum_base[5];
        if(incv == 1){

          for(i = 0; i + 2 <= n; i += 2, v_base += 4){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
            v_2 = fabs(v_base[2]);
            v_3 = fabs(v_base[3]);
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_0_2 = s_0_2 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_0_3 = s_0_3 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            q_2 = q_2 - s_0_2;
            q_3 = q_3 - s_0_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            q_0 = s_1_0;
            q_1 = s_1_1;
            q_2 = s_1_2;
            q_3 = s_1_3;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_1_2 = s_1_2 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_1_3 = s_1_3 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_2_2 = s_2_2 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_2_3 = s_2_3 + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
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
            i += 1, v_base += 2;
          }
        }else{

          for(i = 0; i + 2 <= n; i += 2, v_base += (incv * 4)){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
            v_2 = fabs(v_base[(incv * 2)]);
            v_3 = fabs(v_base[((incv * 2) + 1)]);
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_1 = s_0_1 + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_0_2 = s_0_2 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_0_3 = s_0_3 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            q_2 = q_2 - s_0_2;
            q_3 = q_3 - s_0_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            q_0 = s_1_0;
            q_1 = s_1_1;
            q_2 = s_1_2;
            q_3 = s_1_3;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_1 = s_1_1 + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_1_2 = s_1_2 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_1_3 = s_1_3 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_1 = s_2_1 + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_2_2 = s_2_2 + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_2_3 = s_2_3 + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
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
            i += 1, v_base += (incv * 2);
          }
        }
        q_0 = ((float*)sum)[0];
        s_0_0 = s_0_0 + (s_0_2 - q_0);
        q_0 = ((float*)sum)[1];
        s_0_1 = s_0_1 + (s_0_3 - q_0);
        ((float*)sum)[0] = s_0_0;
        ((float*)sum)[1] = s_0_1;
        q_0 = ((float*)sum)[2];
        s_1_0 = s_1_0 + (s_1_2 - q_0);
        q_0 = ((float*)sum)[3];
        s_1_1 = s_1_1 + (s_1_3 - q_0);
        ((float*)sum)[2] = s_1_0;
        ((float*)sum)[3] = s_1_1;
        q_0 = ((float*)sum)[4];
        s_2_0 = s_2_0 + (s_2_2 - q_0);
        q_0 = ((float*)sum)[5];
        s_2_1 = s_2_1 + (s_2_3 - q_0);
        ((float*)sum)[4] = s_2_0;
        ((float*)sum)[5] = s_2_1;
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float* sum_base = (float*) sum;
        float* v_base = (float*) v;
        float v_0, v_1, v_2, v_3;
        float q_0, q_1, q_2, q_3;
        float s_0, s_1, s_2, s_3;
        float s_buffer[(MAX_FOLD * 4)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = s_buffer[((j * 4) + 2)] = sum_base[(j * 2)];
          s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 3)] = sum_base[((j * 2) + 1)];
        }
        if(incv == 1){

          for(i = 0; i + 2 <= n; i += 2, v_base += 4){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
            v_2 = fabs(v_base[2]);
            v_3 = fabs(v_base[3]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
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
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
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
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 2)] + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 3)] = s_buffer[((j * 4) + 3)] + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.f;
            i += 1, v_base += 2;
          }
        }else{

          for(i = 0; i + 2 <= n; i += 2, v_base += (incv * 4)){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
            v_2 = fabs(v_base[(incv * 2)]);
            v_3 = fabs(v_base[((incv * 2) + 1)]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
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
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
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
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.f;
            tmp_BLP.f = v_2;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 2)] + tmp_BLP.f;
            tmp_BLP.f = v_3;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 3)] = s_buffer[((j * 4) + 3)] + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v_base[0]);
            v_1 = fabs(v_base[1]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.f;
            i += 1, v_base += (incv * 2);
          }
        }
        for(j = 0; j < fold; j += 1){
          q_0 = ((float*)sum)[(j * 2)];
          s_buffer[(j * 4)] = s_buffer[(j * 4)] + (s_buffer[((j * 4) + 2)] - q_0);
          q_0 = ((float*)sum)[((j * 2) + 1)];
          s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + (s_buffer[((j * 4) + 3)] - q_0);
          ((float*)sum)[(j * 2)] = s_buffer[(j * 4)];
          ((float*)sum)[((j * 2) + 1)] = s_buffer[((j * 4) + 1)];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#endif