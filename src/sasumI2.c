#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  void sasumI2(int n, float* v, int incv, int fold, float* sum){
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    __m256 mask_BLP; AVX_BLP_MASKS(mask_BLP);
    float tmp_cons[8] __attribute__((aligned(32)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        __m256 v_0, v_1, v_2, v_3;
        __m256 q_0;
        __m256 s_0_0;
        __m256 s_1_0;
        __m256 s_2_0;

        s_0_0 = _mm256_broadcast_ss(sum);
        s_1_0 = _mm256_broadcast_ss(sum + 1);
        s_2_0 = _mm256_broadcast_ss(sum + 2);
        if(incv == 1){

          for(i = 0; i + 32 <= n; i += 32, v += 32){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_loadu_ps(v + 8), mask_ABS);
            v_2 = _mm256_and_ps(_mm256_loadu_ps(v + 16), mask_ABS);
            v_3 = _mm256_and_ps(_mm256_loadu_ps(v + 24), mask_ABS);
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
          }
          if(i + 16 <= n){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_loadu_ps(v + 8), mask_ABS);
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
            i += 16, v += 16;
          }
          if(i + 8 <= n){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            i += 8, v += 8;
          }
          if(i < n){
            v_0 = _mm256_and_ps(_mm256_set_ps(0, (n - i)>6?v[6]:0, (n - i)>5?v[5]:0, (n - i)>4?v[4]:0, (n - i)>3?v[3]:0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
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

          for(i = 0; i + 32 <= n; i += 32, v += (incv * 32)){
            v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
            v_2 = _mm256_and_ps(_mm256_set_ps(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)], v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]), mask_ABS);
            v_3 = _mm256_and_ps(_mm256_set_ps(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)], v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]), mask_ABS);
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
          }
          if(i + 16 <= n){
            v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            v_1 = _mm256_and_ps(_mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
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
            i += 16, v += (incv * 16);
          }
          if(i + 8 <= n){
            v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            i += 8, v += (incv * 8);
          }
          if(i < n){
            v_0 = _mm256_and_ps(_mm256_set_ps(0, (n - i)>6?v[(incv * 6)]:0, (n - i)>5?v[(incv * 5)]:0, (n - i)>4?v[(incv * 4)]:0, (n - i)>3?v[(incv * 3)]:0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
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
        s_0_0 = _mm256_sub_ps(s_0_0, _mm256_set_ps(sum[0], sum[0], sum[0], sum[0], sum[0], sum[0], sum[0], 0));
        _mm256_store_ps(tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3] + tmp_cons[4] + tmp_cons[5] + tmp_cons[6] + tmp_cons[7];
        s_1_0 = _mm256_sub_ps(s_1_0, _mm256_set_ps(sum[1], sum[1], sum[1], sum[1], sum[1], sum[1], sum[1], 0));
        _mm256_store_ps(tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3] + tmp_cons[4] + tmp_cons[5] + tmp_cons[6] + tmp_cons[7];
        s_2_0 = _mm256_sub_ps(s_2_0, _mm256_set_ps(sum[2], sum[2], sum[2], sum[2], sum[2], sum[2], sum[2], 0));
        _mm256_store_ps(tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3] + tmp_cons[4] + tmp_cons[5] + tmp_cons[6] + tmp_cons[7];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        __m256 v_0;
        __m256 q_0;
        __m256 s_0;
        __m256 s_buffer[MAX_FOLD];

        for(j = 0; j < fold; j += 1){
          s_buffer[j] = _mm256_broadcast_ss(sum + j);
        }
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v += 8){
            v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[j];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[j] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[j] = _mm256_add_ps(s_buffer[j], _mm256_or_ps(v_0, mask_BLP));
          }
          if(i < n){
            v_0 = _mm256_and_ps(_mm256_set_ps(0, (n - i)>6?v[6]:0, (n - i)>5?v[5]:0, (n - i)>4?v[4]:0, (n - i)>3?v[3]:0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[j];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[j] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[j] = _mm256_add_ps(s_buffer[j], _mm256_or_ps(v_0, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v += (incv * 8)){
            v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[j];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[j] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[j] = _mm256_add_ps(s_buffer[j], _mm256_or_ps(v_0, mask_BLP));
          }
          if(i < n){
            v_0 = _mm256_and_ps(_mm256_set_ps(0, (n - i)>6?v[(incv * 6)]:0, (n - i)>5?v[(incv * 5)]:0, (n - i)>4?v[(incv * 4)]:0, (n - i)>3?v[(incv * 3)]:0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[j];
              q_0 = _mm256_add_ps(s_0, _mm256_or_ps(v_0, mask_BLP));
              s_buffer[j] = q_0;
              q_0 = _mm256_sub_ps(s_0, q_0);
              v_0 = _mm256_add_ps(v_0, q_0);
            }
            s_buffer[j] = _mm256_add_ps(s_buffer[j], _mm256_or_ps(v_0, mask_BLP));
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[j] = _mm256_sub_ps(s_buffer[j], _mm256_set_ps(sum[j], sum[j], sum[j], sum[j], sum[j], sum[j], sum[j], 0));
          _mm256_store_ps(tmp_cons, s_buffer[j]);
          sum[j] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3] + tmp_cons[4] + tmp_cons[5] + tmp_cons[6] + tmp_cons[7];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#elif defined( __SSE2__ )
  void sasumI2(int n, float* v, int incv, int fold, float* sum){
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    __m128 mask_BLP; SSE_BLP_MASKS(mask_BLP);
    float tmp_cons[4] __attribute__((aligned(16)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        __m128 v_0, v_1;
        __m128 q_0;
        __m128 s_0_0;
        __m128 s_1_0;
        __m128 s_2_0;

        s_0_0 = _mm_load1_ps(sum);
        s_1_0 = _mm_load1_ps(sum + 1);
        s_2_0 = _mm_load1_ps(sum + 2);
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v += 8){
            v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
            v_1 = _mm_and_ps(_mm_loadu_ps(v + 4), mask_ABS);
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
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            i += 4, v += 4;
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
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

          for(i = 0; i + 8 <= n; i += 8, v += (incv * 8)){
            v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            v_1 = _mm_and_ps(_mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
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
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            i += 4, v += (incv * 4);
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
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
        s_0_0 = _mm_sub_ps(s_0_0, _mm_set_ps(sum[0], sum[0], sum[0], 0));
        _mm_store_ps(tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        s_1_0 = _mm_sub_ps(s_1_0, _mm_set_ps(sum[1], sum[1], sum[1], 0));
        _mm_store_ps(tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        s_2_0 = _mm_sub_ps(s_2_0, _mm_set_ps(sum[2], sum[2], sum[2], 0));
        _mm_store_ps(tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        __m128 v_0, v_1;
        __m128 q_0, q_1;
        __m128 s_0, s_1;
        __m128 s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = _mm_load1_ps(sum + j);
        }
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v += 8){
            v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
            v_1 = _mm_and_ps(_mm_loadu_ps(v + 4), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm_add_ps(s_buffer[((j * 2) + 1)], _mm_or_ps(v_1, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
            i += 4, v += 4;
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
          }
        }else{

          for(i = 0; i + 8 <= n; i += 8, v += (incv * 8)){
            v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            v_1 = _mm_and_ps(_mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              q_1 = _mm_add_ps(s_1, _mm_or_ps(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_ps(s_0, q_0);
              q_1 = _mm_sub_ps(s_1, q_1);
              v_0 = _mm_add_ps(v_0, q_0);
              v_1 = _mm_add_ps(v_1, q_1);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm_add_ps(s_buffer[((j * 2) + 1)], _mm_or_ps(v_1, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
            i += 4, v += (incv * 4);
          }
          if(i < n){
            v_0 = _mm_and_ps(_mm_set_ps(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = _mm_sub_ps(s_buffer[(j * 2)], _mm_set_ps(sum[j], sum[j], sum[j], 0));
          q_0 = _mm_load1_ps(sum + j);
          s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_sub_ps(s_buffer[((j * 2) + 1)], q_0));
          _mm_store_ps(tmp_cons, s_buffer[(j * 2)]);
          sum[j] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#else
  void sasumI2(int n, float* v, int incv, int fold, float* sum){
    i_float tmp_BLP;
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float v_0, v_1;
        float q_0;
        float s_0_0;
        float s_1_0;
        float s_2_0;

        s_0_0 = sum[0];
        s_1_0 = sum[1];
        s_2_0 = sum[2];
        if(incv == 1){

          for(i = 0; i + 2 <= n; i += 2, v += 2){
            v_0 = fabs(v[0]);
            v_1 = fabs(v[1]);
            q_0 = s_0_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            v_0 = v_0 + q_0;
            q_0 = s_1_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            v_0 = v_0 + q_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            q_0 = s_0_0;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            v_1 = v_1 + q_0;
            q_0 = s_1_0;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            v_1 = v_1 + q_0;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v[0]);
            q_0 = s_0_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            v_0 = v_0 + q_0;
            q_0 = s_1_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            v_0 = v_0 + q_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            i += 1, v += 1;
          }
        }else{

          for(i = 0; i + 2 <= n; i += 2, v += (incv * 2)){
            v_0 = fabs(v[0]);
            v_1 = fabs(v[incv]);
            q_0 = s_0_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            v_0 = v_0 + q_0;
            q_0 = s_1_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            v_0 = v_0 + q_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            q_0 = s_0_0;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            v_1 = v_1 + q_0;
            q_0 = s_1_0;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            v_1 = v_1 + q_0;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v[0]);
            q_0 = s_0_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_0_0 = s_0_0 + tmp_BLP.f;
            q_0 = q_0 - s_0_0;
            v_0 = v_0 + q_0;
            q_0 = s_1_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_1_0 = s_1_0 + tmp_BLP.f;
            q_0 = q_0 - s_1_0;
            v_0 = v_0 + q_0;
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_2_0 = s_2_0 + tmp_BLP.f;
            i += 1, v += incv;
          }
        }
        sum[0] = s_0_0;
        sum[1] = s_1_0;
        sum[2] = s_2_0;
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        float v_0, v_1;
        float q_0, q_1;
        float s_0, s_1;
        float s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = sum[j];
        }
        if(incv == 1){

          for(i = 0; i + 2 <= n; i += 2, v += 2){
            v_0 = fabs(v[0]);
            v_1 = fabs(v[1]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 2) + 1)] = s_buffer[((j * 2) + 1)] + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v[0]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              s_buffer[(j * 2)] = q_0;
              q_0 = s_0 - q_0;
              v_0 = v_0 + q_0;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.f;
            i += 1, v += 1;
          }
        }else{

          for(i = 0; i + 2 <= n; i += 2, v += (incv * 2)){
            v_0 = fabs(v[0]);
            v_1 = fabs(v[incv]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              tmp_BLP.f = v_1;
              tmp_BLP.i |= 1;
              q_1 = s_1 + tmp_BLP.f;
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.f;
            tmp_BLP.f = v_1;
            tmp_BLP.i |= 1;
            s_buffer[((j * 2) + 1)] = s_buffer[((j * 2) + 1)] + tmp_BLP.f;
          }
          if(i + 1 <= n){
            v_0 = fabs(v[0]);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              tmp_BLP.f = v_0;
              tmp_BLP.i |= 1;
              q_0 = s_0 + tmp_BLP.f;
              s_buffer[(j * 2)] = q_0;
              q_0 = s_0 - q_0;
              v_0 = v_0 + q_0;
            }
            tmp_BLP.f = v_0;
            tmp_BLP.i |= 1;
            s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.f;
            i += 1, v += incv;
          }
        }
        for(j = 0; j < fold; j += 1){
          q_0 = sum[j];
          s_buffer[(j * 2)] = s_buffer[(j * 2)] + (s_buffer[((j * 2) + 1)] - q_0);
          sum[j] = s_buffer[(j * 2)];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#endif