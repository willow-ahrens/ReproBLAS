#include <float.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <emmintrin.h>


#if defined( __AVX__ )
  void sdotI2(int n, float* v, int incv, float* y, int incy, int fold, float* sum){
    __m256 mask_BLP; AVX_BLP_MASKS(mask_BLP);
    float tmp[8] __attribute__((aligned(32)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        __m256 v_0, v_1, v_2, v_3;
        __m256 y_0, y_1, y_2, y_3;
        __m256 q_0;
        __m256 s_0_0;
        __m256 s_1_0;
        __m256 s_2_0;

        s_0_0 = _mm256_broadcast_ss(sum);
        s_1_0 = _mm256_broadcast_ss(sum + 1);
        s_2_0 = _mm256_broadcast_ss(sum + 2);
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
            v_0 = _mm256_mul_ps(v_0, y_0);
            v_1 = _mm256_mul_ps(v_1, y_1);
            v_2 = _mm256_mul_ps(v_2, y_2);
            v_3 = _mm256_mul_ps(v_3, y_3);
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
            v_0 = _mm256_loadu_ps(v);
            v_1 = _mm256_loadu_ps(v + 8);
            y_0 = _mm256_loadu_ps(y);
            y_1 = _mm256_loadu_ps(y + 8);
            v_0 = _mm256_mul_ps(v_0, y_0);
            v_1 = _mm256_mul_ps(v_1, y_1);
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
            i += 16, v += 16, y += 16;
          }
          if(i + 8 <= n){
            v_0 = _mm256_loadu_ps(v);
            y_0 = _mm256_loadu_ps(y);
            v_0 = _mm256_mul_ps(v_0, y_0);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            i += 8, v += 8, y += 8;
          }
          if(i < n){
            v_0 = _mm256_set_ps(0, (n - i)>6?v[6]:0, (n - i)>5?v[5]:0, (n - i)>4?v[4]:0, (n - i)>3?v[3]:0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
            y_0 = _mm256_set_ps(0, (n - i)>6?y[6]:0, (n - i)>5?y[5]:0, (n - i)>4?y[4]:0, (n - i)>3?y[3]:0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
            v_0 = _mm256_mul_ps(v_0, y_0);
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

          for(i = 0; i + 32 <= n; i += 32, v += (incv * 32), y += (incy * 32)){
            v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            v_1 = _mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
            v_2 = _mm256_set_ps(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)], v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]);
            v_3 = _mm256_set_ps(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)], v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]);
            y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            y_1 = _mm256_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)], y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
            y_2 = _mm256_set_ps(y[(incy * 23)], y[(incy * 22)], y[(incy * 21)], y[(incy * 20)], y[(incy * 19)], y[(incy * 18)], y[(incy * 17)], y[(incy * 16)]);
            y_3 = _mm256_set_ps(y[(incy * 31)], y[(incy * 30)], y[(incy * 29)], y[(incy * 28)], y[(incy * 27)], y[(incy * 26)], y[(incy * 25)], y[(incy * 24)]);
            v_0 = _mm256_mul_ps(v_0, y_0);
            v_1 = _mm256_mul_ps(v_1, y_1);
            v_2 = _mm256_mul_ps(v_2, y_2);
            v_3 = _mm256_mul_ps(v_3, y_3);
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
            v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            v_1 = _mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
            y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            y_1 = _mm256_set_ps(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)], y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
            v_0 = _mm256_mul_ps(v_0, y_0);
            v_1 = _mm256_mul_ps(v_1, y_1);
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
            i += 16, v += (incv * 16), y += (incy * 16);
          }
          if(i + 8 <= n){
            v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            v_0 = _mm256_mul_ps(v_0, y_0);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_ps(s_0_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_0_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_ps(s_1_0, _mm256_or_ps(v_0, mask_BLP));
            q_0 = _mm256_sub_ps(q_0, s_1_0);
            v_0 = _mm256_add_ps(v_0, q_0);
            s_2_0 = _mm256_add_ps(s_2_0, _mm256_or_ps(v_0, mask_BLP));
            i += 8, v += (incv * 8), y += (incy * 8);
          }
          if(i < n){
            v_0 = _mm256_set_ps(0, (n - i)>6?v[(incv * 6)]:0, (n - i)>5?v[(incv * 5)]:0, (n - i)>4?v[(incv * 4)]:0, (n - i)>3?v[(incv * 3)]:0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
            y_0 = _mm256_set_ps(0, (n - i)>6?y[(incy * 6)]:0, (n - i)>5?y[(incy * 5)]:0, (n - i)>4?y[(incy * 4)]:0, (n - i)>3?y[(incy * 3)]:0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
            v_0 = _mm256_mul_ps(v_0, y_0);
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
        _mm256_store_ps(tmp, s_0_0);
        sum[0] = tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
        s_1_0 = _mm256_sub_ps(s_1_0, _mm256_set_ps(sum[1], sum[1], sum[1], sum[1], sum[1], sum[1], sum[1], 0));
        _mm256_store_ps(tmp, s_1_0);
        sum[1] = tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
        s_2_0 = _mm256_sub_ps(s_2_0, _mm256_set_ps(sum[2], sum[2], sum[2], sum[2], sum[2], sum[2], sum[2], 0));
        _mm256_store_ps(tmp, s_2_0);
        sum[2] = tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        __m256 v_0;
        __m256 y_0;
        __m256 q_0;
        __m256 s_0;
        __m256 s_buffer[MAX_FOLD];

        for(j = 0; j < fold; j += 1){
          s_buffer[j] = _mm256_broadcast_ss(sum + j);
        }
        if(incv == 1 && incy == 1){

          for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
            v_0 = _mm256_loadu_ps(v);
            y_0 = _mm256_loadu_ps(y);
            v_0 = _mm256_mul_ps(v_0, y_0);
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
            v_0 = _mm256_set_ps(0, (n - i)>6?v[6]:0, (n - i)>5?v[5]:0, (n - i)>4?v[4]:0, (n - i)>3?v[3]:0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
            y_0 = _mm256_set_ps(0, (n - i)>6?y[6]:0, (n - i)>5?y[5]:0, (n - i)>4?y[4]:0, (n - i)>3?y[3]:0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
            v_0 = _mm256_mul_ps(v_0, y_0);
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

          for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
            v_0 = _mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            y_0 = _mm256_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)], y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            v_0 = _mm256_mul_ps(v_0, y_0);
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
            v_0 = _mm256_set_ps(0, (n - i)>6?v[(incv * 6)]:0, (n - i)>5?v[(incv * 5)]:0, (n - i)>4?v[(incv * 4)]:0, (n - i)>3?v[(incv * 3)]:0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
            y_0 = _mm256_set_ps(0, (n - i)>6?y[(incy * 6)]:0, (n - i)>5?y[(incy * 5)]:0, (n - i)>4?y[(incy * 4)]:0, (n - i)>3?y[(incy * 3)]:0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
            v_0 = _mm256_mul_ps(v_0, y_0);
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
          _mm256_store_ps(tmp, s_buffer[j]);
          sum[j] = tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#elif defined( __SSE2__ )
  void sdotI2(int n, float* v, int incv, float* y, int incy, int fold, float* sum){
    __m128 mask_BLP; SSE_BLP_MASKS(mask_BLP);
    float tmp[4] __attribute__((aligned(16)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        __m128 v_0, v_1;
        __m128 y_0, y_1;
        __m128 q_0;
        __m128 s_0_0;
        __m128 s_1_0;
        __m128 s_2_0;

        s_0_0 = _mm_load1_ps(sum);
        s_1_0 = _mm_load1_ps(sum + 1);
        s_2_0 = _mm_load1_ps(sum + 2);
        if(incv == 1 && incy == 1){

          for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
            v_0 = _mm_loadu_ps(v);
            v_1 = _mm_loadu_ps(v + 4);
            y_0 = _mm_loadu_ps(y);
            y_1 = _mm_loadu_ps(y + 4);
            v_0 = _mm_mul_ps(v_0, y_0);
            v_1 = _mm_mul_ps(v_1, y_1);
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
            v_0 = _mm_loadu_ps(v);
            y_0 = _mm_loadu_ps(y);
            v_0 = _mm_mul_ps(v_0, y_0);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            i += 4, v += 4, y += 4;
          }
          if(i < n){
            v_0 = _mm_set_ps(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
            y_0 = _mm_set_ps(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
            v_0 = _mm_mul_ps(v_0, y_0);
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

          for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
            v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            v_1 = _mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
            y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            y_1 = _mm_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
            v_0 = _mm_mul_ps(v_0, y_0);
            v_1 = _mm_mul_ps(v_1, y_1);
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
            v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            v_0 = _mm_mul_ps(v_0, y_0);
            q_0 = s_0_0;
            s_0_0 = _mm_add_ps(s_0_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_0_0);
            v_0 = _mm_add_ps(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_ps(s_1_0, _mm_or_ps(v_0, mask_BLP));
            q_0 = _mm_sub_ps(q_0, s_1_0);
            v_0 = _mm_add_ps(v_0, q_0);
            s_2_0 = _mm_add_ps(s_2_0, _mm_or_ps(v_0, mask_BLP));
            i += 4, v += (incv * 4), y += (incy * 4);
          }
          if(i < n){
            v_0 = _mm_set_ps(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
            y_0 = _mm_set_ps(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
            v_0 = _mm_mul_ps(v_0, y_0);
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
        _mm_store_ps(tmp, s_0_0);
        sum[0] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        s_1_0 = _mm_sub_ps(s_1_0, _mm_set_ps(sum[1], sum[1], sum[1], 0));
        _mm_store_ps(tmp, s_1_0);
        sum[1] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        s_2_0 = _mm_sub_ps(s_2_0, _mm_set_ps(sum[2], sum[2], sum[2], 0));
        _mm_store_ps(tmp, s_2_0);
        sum[2] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        __m128 v_0, v_1;
        __m128 y_0, y_1;
        __m128 q_0, q_1;
        __m128 s_0, s_1;
        __m128 s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = _mm_load1_ps(sum + j);
        }
        if(incv == 1 && incy == 1){

          for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
            v_0 = _mm_loadu_ps(v);
            v_1 = _mm_loadu_ps(v + 4);
            y_0 = _mm_loadu_ps(y);
            y_1 = _mm_loadu_ps(y + 4);
            v_0 = _mm_mul_ps(v_0, y_0);
            v_1 = _mm_mul_ps(v_1, y_1);
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
            v_0 = _mm_loadu_ps(v);
            y_0 = _mm_loadu_ps(y);
            v_0 = _mm_mul_ps(v_0, y_0);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
            i += 4, v += 4, y += 4;
          }
          if(i < n){
            v_0 = _mm_set_ps(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
            y_0 = _mm_set_ps(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
            v_0 = _mm_mul_ps(v_0, y_0);
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

          for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
            v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            v_1 = _mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
            y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            y_1 = _mm_set_ps(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
            v_0 = _mm_mul_ps(v_0, y_0);
            v_1 = _mm_mul_ps(v_1, y_1);
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
            v_0 = _mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
            y_0 = _mm_set_ps(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
            v_0 = _mm_mul_ps(v_0, y_0);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_ps(s_0, _mm_or_ps(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_ps(s_0, q_0);
              v_0 = _mm_add_ps(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_ps(s_buffer[(j * 2)], _mm_or_ps(v_0, mask_BLP));
            i += 4, v += (incv * 4), y += (incy * 4);
          }
          if(i < n){
            v_0 = _mm_set_ps(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
            y_0 = _mm_set_ps(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
            v_0 = _mm_mul_ps(v_0, y_0);
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
          _mm_store_ps(tmp, s_buffer[(j * 2)]);
          sum[j] = tmp[0] + tmp[1] + tmp[2] + tmp[3];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#else
  void sdotI2(int n, float* v, int incv, float* y, int incy, int fold, float* sum){
    i_float tmp_BLP;
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        float v_0, v_1;
        float y_0, y_1;
        float q_0;
        float s_0_0;
        float s_1_0;
        float s_2_0;

        s_0_0 = sum[0];
        s_1_0 = sum[1];
        s_2_0 = sum[2];
        if(incv == 1 && incy == 1){

          for(i = 0; i + 2 <= n; i += 2, v += 2, y += 2){
            v_0 = v[0];
            v_1 = v[1];
            y_0 = y[0];
            y_1 = y[1];
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_1;
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
            v_0 = v[0];
            y_0 = y[0];
            v_0 = v_0 * y_0;
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
            i += 1, v += 1, y += 1;
          }
        }else{

          for(i = 0; i + 2 <= n; i += 2, v += (incv * 2), y += (incy * 2)){
            v_0 = v[0];
            v_1 = v[incv];
            y_0 = y[0];
            y_1 = y[incy];
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_1;
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
            v_0 = v[0];
            y_0 = y[0];
            v_0 = v_0 * y_0;
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
            i += 1, v += incv, y += incy;
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
        float y_0, y_1;
        float q_0, q_1;
        float s_0, s_1;
        float s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = sum[j];
        }
        if(incv == 1 && incy == 1){

          for(i = 0; i + 2 <= n; i += 2, v += 2, y += 2){
            v_0 = v[0];
            v_1 = v[1];
            y_0 = y[0];
            y_1 = y[1];
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_1;
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
            v_0 = v[0];
            y_0 = y[0];
            v_0 = v_0 * y_0;
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
            i += 1, v += 1, y += 1;
          }
        }else{

          for(i = 0; i + 2 <= n; i += 2, v += (incv * 2), y += (incy * 2)){
            v_0 = v[0];
            v_1 = v[incv];
            y_0 = y[0];
            y_1 = y[incy];
            v_0 = v_0 * y_0;
            v_1 = v_1 * y_1;
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
            v_0 = v[0];
            y_0 = y[0];
            v_0 = v_0 * y_0;
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
            i += 1, v += incv, y += incy;
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