#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  void ddotI2(int n, double* v, int incv, double* y, int incy, int fold, double* sum){
    __m256d mask_BLP; AVX_BLP_MASKD(mask_BLP);
    double tmp_cons[4] __attribute__((aligned(32)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        __m256d v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        __m256d y_0, y_1, y_2, y_3, y_4, y_5, y_6, y_7;
        __m256d q_0, q_1;
        __m256d s_0_0, s_0_1;
        __m256d s_1_0, s_1_1;
        __m256d s_2_0, s_2_1;

        s_0_0 = s_0_1 = _mm256_broadcast_sd(sum);
        s_1_0 = s_1_1 = _mm256_broadcast_sd(sum + 1);
        s_2_0 = s_2_1 = _mm256_broadcast_sd(sum + 2);
        if(incv == 1){
          if(incy == 1){

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
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              v_4 = _mm256_mul_pd(v_4, y_4);
              v_5 = _mm256_mul_pd(v_5, y_5);
              v_6 = _mm256_mul_pd(v_6, y_6);
              v_7 = _mm256_mul_pd(v_7, y_7);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_4, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_4, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_4, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_6, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_6, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_6, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_7, mask_BLP));
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
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              i += 16, v += 16, y += 16;
            }
            if(i + 8 <= n){
              v_0 = _mm256_loadu_pd(v);
              v_1 = _mm256_loadu_pd(v + 4);
              y_0 = _mm256_loadu_pd(y);
              y_1 = _mm256_loadu_pd(y + 4);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              i += 8, v += 8, y += 8;
            }
            if(i + 4 <= n){
              v_0 = _mm256_loadu_pd(v);
              y_0 = _mm256_loadu_pd(y);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += 4, y += 4;
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 32 <= n; i += 32, v += 32, y += (incy * 32)){
              v_0 = _mm256_loadu_pd(v);
              v_1 = _mm256_loadu_pd(v + 4);
              v_2 = _mm256_loadu_pd(v + 8);
              v_3 = _mm256_loadu_pd(v + 12);
              v_4 = _mm256_loadu_pd(v + 16);
              v_5 = _mm256_loadu_pd(v + 20);
              v_6 = _mm256_loadu_pd(v + 24);
              v_7 = _mm256_loadu_pd(v + 28);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
              y_2 = _mm256_set_pd(y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
              y_3 = _mm256_set_pd(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)]);
              y_4 = _mm256_set_pd(y[(incy * 19)], y[(incy * 18)], y[(incy * 17)], y[(incy * 16)]);
              y_5 = _mm256_set_pd(y[(incy * 23)], y[(incy * 22)], y[(incy * 21)], y[(incy * 20)]);
              y_6 = _mm256_set_pd(y[(incy * 27)], y[(incy * 26)], y[(incy * 25)], y[(incy * 24)]);
              y_7 = _mm256_set_pd(y[(incy * 31)], y[(incy * 30)], y[(incy * 29)], y[(incy * 28)]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              v_4 = _mm256_mul_pd(v_4, y_4);
              v_5 = _mm256_mul_pd(v_5, y_5);
              v_6 = _mm256_mul_pd(v_6, y_6);
              v_7 = _mm256_mul_pd(v_7, y_7);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_4, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_4, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_4, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_6, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_6, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_6, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_7, mask_BLP));
            }
            if(i + 16 <= n){
              v_0 = _mm256_loadu_pd(v);
              v_1 = _mm256_loadu_pd(v + 4);
              v_2 = _mm256_loadu_pd(v + 8);
              v_3 = _mm256_loadu_pd(v + 12);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
              y_2 = _mm256_set_pd(y[(incy * 11)], y[(incy * 10)], y[(incy * 9)], y[(incy * 8)]);
              y_3 = _mm256_set_pd(y[(incy * 15)], y[(incy * 14)], y[(incy * 13)], y[(incy * 12)]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              i += 16, v += 16, y += (incy * 16);
            }
            if(i + 8 <= n){
              v_0 = _mm256_loadu_pd(v);
              v_1 = _mm256_loadu_pd(v + 4);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              i += 8, v += 8, y += (incy * 8);
            }
            if(i + 4 <= n){
              v_0 = _mm256_loadu_pd(v);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += 4, y += (incy * 4);
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
            }
          }
        }else{
          if(incy == 1){

            for(i = 0; i + 32 <= n; i += 32, v += (incv * 32), y += 32){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
              v_2 = _mm256_set_pd(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
              v_3 = _mm256_set_pd(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]);
              v_4 = _mm256_set_pd(v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]);
              v_5 = _mm256_set_pd(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)]);
              v_6 = _mm256_set_pd(v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]);
              v_7 = _mm256_set_pd(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)]);
              y_0 = _mm256_loadu_pd(y);
              y_1 = _mm256_loadu_pd(y + 4);
              y_2 = _mm256_loadu_pd(y + 8);
              y_3 = _mm256_loadu_pd(y + 12);
              y_4 = _mm256_loadu_pd(y + 16);
              y_5 = _mm256_loadu_pd(y + 20);
              y_6 = _mm256_loadu_pd(y + 24);
              y_7 = _mm256_loadu_pd(y + 28);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              v_4 = _mm256_mul_pd(v_4, y_4);
              v_5 = _mm256_mul_pd(v_5, y_5);
              v_6 = _mm256_mul_pd(v_6, y_6);
              v_7 = _mm256_mul_pd(v_7, y_7);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_4, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_4, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_4, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_6, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_6, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_6, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_7, mask_BLP));
            }
            if(i + 16 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
              v_2 = _mm256_set_pd(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]);
              v_3 = _mm256_set_pd(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]);
              y_0 = _mm256_loadu_pd(y);
              y_1 = _mm256_loadu_pd(y + 4);
              y_2 = _mm256_loadu_pd(y + 8);
              y_3 = _mm256_loadu_pd(y + 12);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              i += 16, v += (incv * 16), y += 16;
            }
            if(i + 8 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
              y_0 = _mm256_loadu_pd(y);
              y_1 = _mm256_loadu_pd(y + 4);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              i += 8, v += (incv * 8), y += 8;
            }
            if(i + 4 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              y_0 = _mm256_loadu_pd(y);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += (incv * 4), y += 4;
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
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
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              v_4 = _mm256_mul_pd(v_4, y_4);
              v_5 = _mm256_mul_pd(v_5, y_5);
              v_6 = _mm256_mul_pd(v_6, y_6);
              v_7 = _mm256_mul_pd(v_7, y_7);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_4, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_4, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_4 = _mm256_add_pd(v_4, q_0);
              v_5 = _mm256_add_pd(v_5, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_4, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_5, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_6, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_6, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_7, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_6 = _mm256_add_pd(v_6, q_0);
              v_7 = _mm256_add_pd(v_7, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_6, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_7, mask_BLP));
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
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              v_2 = _mm256_mul_pd(v_2, y_2);
              v_3 = _mm256_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_2, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_2, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_3, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_2, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_3, mask_BLP));
              i += 16, v += (incv * 16), y += (incy * 16);
            }
            if(i + 8 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              s_0_1 = _mm256_add_pd(s_0_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              q_1 = _mm256_sub_pd(q_1, s_0_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              s_1_1 = _mm256_add_pd(s_1_1, _mm256_or_pd(v_1, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              q_1 = _mm256_sub_pd(q_1, s_1_1);
              v_0 = _mm256_add_pd(v_0, q_0);
              v_1 = _mm256_add_pd(v_1, q_1);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              s_2_1 = _mm256_add_pd(s_2_1, _mm256_or_pd(v_1, mask_BLP));
              i += 8, v += (incv * 8), y += (incy * 8);
            }
            if(i + 4 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += (incv * 4), y += (incy * 4);
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_0_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
              q_0 = _mm256_sub_pd(q_0, s_1_0);
              v_0 = _mm256_add_pd(v_0, q_0);
              s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
            }
          }
        }
        s_0_0 = _mm256_sub_pd(s_0_0, _mm256_set_pd(sum[0], sum[0], sum[0], 0));
        q_0 = _mm256_broadcast_sd(sum);
        s_0_0 = _mm256_add_pd(s_0_0, _mm256_sub_pd(s_0_1, q_0));
        _mm256_store_pd(tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        s_1_0 = _mm256_sub_pd(s_1_0, _mm256_set_pd(sum[1], sum[1], sum[1], 0));
        q_0 = _mm256_broadcast_sd(sum + 1);
        s_1_0 = _mm256_add_pd(s_1_0, _mm256_sub_pd(s_1_1, q_0));
        _mm256_store_pd(tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        s_2_0 = _mm256_sub_pd(s_2_0, _mm256_set_pd(sum[2], sum[2], sum[2], 0));
        q_0 = _mm256_broadcast_sd(sum + 2);
        s_2_0 = _mm256_add_pd(s_2_0, _mm256_sub_pd(s_2_1, q_0));
        _mm256_store_pd(tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        __m256d v_0, v_1;
        __m256d y_0, y_1;
        __m256d q_0, q_1;
        __m256d s_0, s_1;
        __m256d s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = _mm256_broadcast_sd(sum + j);
        }
        if(incv == 1){
          if(incy == 1){

            for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
              v_0 = _mm256_loadu_pd(v);
              v_1 = _mm256_loadu_pd(v + 4);
              y_0 = _mm256_loadu_pd(y);
              y_1 = _mm256_loadu_pd(y + 4);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                s_1 = s_buffer[((j * 2) + 1)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                q_1 = _mm256_add_pd(s_1, _mm256_or_pd(v_1, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = _mm256_sub_pd(s_0, q_0);
                q_1 = _mm256_sub_pd(s_1, q_1);
                v_0 = _mm256_add_pd(v_0, q_0);
                v_1 = _mm256_add_pd(v_1, q_1);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              s_buffer[((j * 2) + 1)] = _mm256_add_pd(s_buffer[((j * 2) + 1)], _mm256_or_pd(v_1, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm256_loadu_pd(v);
              y_0 = _mm256_loadu_pd(y);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += 4, y += 4;
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 8 <= n; i += 8, v += 8, y += (incy * 8)){
              v_0 = _mm256_loadu_pd(v);
              v_1 = _mm256_loadu_pd(v + 4);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                s_1 = s_buffer[((j * 2) + 1)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                q_1 = _mm256_add_pd(s_1, _mm256_or_pd(v_1, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = _mm256_sub_pd(s_0, q_0);
                q_1 = _mm256_sub_pd(s_1, q_1);
                v_0 = _mm256_add_pd(v_0, q_0);
                v_1 = _mm256_add_pd(v_1, q_1);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              s_buffer[((j * 2) + 1)] = _mm256_add_pd(s_buffer[((j * 2) + 1)], _mm256_or_pd(v_1, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm256_loadu_pd(v);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += 4, y += (incy * 4);
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
            }
          }
        }else{
          if(incy == 1){

            for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += 8){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
              y_0 = _mm256_loadu_pd(y);
              y_1 = _mm256_loadu_pd(y + 4);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                s_1 = s_buffer[((j * 2) + 1)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                q_1 = _mm256_add_pd(s_1, _mm256_or_pd(v_1, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = _mm256_sub_pd(s_0, q_0);
                q_1 = _mm256_sub_pd(s_1, q_1);
                v_0 = _mm256_add_pd(v_0, q_0);
                v_1 = _mm256_add_pd(v_1, q_1);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              s_buffer[((j * 2) + 1)] = _mm256_add_pd(s_buffer[((j * 2) + 1)], _mm256_or_pd(v_1, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              y_0 = _mm256_loadu_pd(y);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += (incv * 4), y += 4;
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[2]:0, (n - i)>1?y[1]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              v_1 = _mm256_set_pd(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              y_1 = _mm256_set_pd(y[(incy * 7)], y[(incy * 6)], y[(incy * 5)], y[(incy * 4)]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              v_1 = _mm256_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                s_1 = s_buffer[((j * 2) + 1)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                q_1 = _mm256_add_pd(s_1, _mm256_or_pd(v_1, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = _mm256_sub_pd(s_0, q_0);
                q_1 = _mm256_sub_pd(s_1, q_1);
                v_0 = _mm256_add_pd(v_0, q_0);
                v_1 = _mm256_add_pd(v_1, q_1);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              s_buffer[((j * 2) + 1)] = _mm256_add_pd(s_buffer[((j * 2) + 1)], _mm256_or_pd(v_1, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm256_set_pd(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]);
              y_0 = _mm256_set_pd(y[(incy * 3)], y[(incy * 2)], y[incy], y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
              i += 4, v += (incv * 4), y += (incy * 4);
            }
            if(i < n){
              v_0 = _mm256_set_pd(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]);
              y_0 = _mm256_set_pd(0, (n - i)>2?y[(incy * 2)]:0, (n - i)>1?y[incy]:0, y[0]);
              v_0 = _mm256_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
                s_buffer[(j * 2)] = q_0;
                q_0 = _mm256_sub_pd(s_0, q_0);
                v_0 = _mm256_add_pd(v_0, q_0);
              }
              s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
            }
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = _mm256_sub_pd(s_buffer[(j * 2)], _mm256_set_pd(sum[j], sum[j], sum[j], 0));
          q_0 = _mm256_broadcast_sd(sum + j);
          s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_sub_pd(s_buffer[((j * 2) + 1)], q_0));
          _mm256_store_pd(tmp_cons, s_buffer[(j * 2)]);
          sum[j] = tmp_cons[0] + tmp_cons[1] + tmp_cons[2] + tmp_cons[3];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#elif defined( __SSE2__ )
  void ddotI2(int n, double* v, int incv, double* y, int incy, int fold, double* sum){
    __m128d mask_BLP; SSE_BLP_MASKD(mask_BLP);
    double tmp_cons[2] __attribute__((aligned(16)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        __m128d v_0, v_1, v_2, v_3;
        __m128d y_0, y_1, y_2, y_3;
        __m128d q_0, q_1, q_2, q_3;
        __m128d s_0_0, s_0_1, s_0_2, s_0_3;
        __m128d s_1_0, s_1_1, s_1_2, s_1_3;
        __m128d s_2_0, s_2_1, s_2_2, s_2_3;

        s_0_0 = s_0_1 = s_0_2 = s_0_3 = _mm_load1_pd(sum);
        s_1_0 = s_1_1 = s_1_2 = s_1_3 = _mm_load1_pd(sum + 1);
        s_2_0 = s_2_1 = s_2_2 = s_2_3 = _mm_load1_pd(sum + 2);
        if(incv == 1){
          if(incy == 1){

            for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              v_2 = _mm_loadu_pd(v + 4);
              v_3 = _mm_loadu_pd(v + 6);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              y_2 = _mm_loadu_pd(y + 4);
              y_3 = _mm_loadu_pd(y + 6);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              q_2 = s_0_2;
              q_3 = s_0_3;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              s_0_2 = _mm_add_pd(s_0_2, _mm_or_pd(v_2, mask_BLP));
              s_0_3 = _mm_add_pd(s_0_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              q_2 = _mm_sub_pd(q_2, s_0_2);
              q_3 = _mm_sub_pd(q_3, s_0_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              q_0 = s_1_0;
              q_1 = s_1_1;
              q_2 = s_1_2;
              q_3 = s_1_3;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              s_1_2 = _mm_add_pd(s_1_2, _mm_or_pd(v_2, mask_BLP));
              s_1_3 = _mm_add_pd(s_1_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              q_2 = _mm_sub_pd(q_2, s_1_2);
              q_3 = _mm_sub_pd(q_3, s_1_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              s_2_2 = _mm_add_pd(s_2_2, _mm_or_pd(v_2, mask_BLP));
              s_2_3 = _mm_add_pd(s_2_3, _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              i += 4, v += 4, y += 4;
            }
            if(i + 2 <= n){
              v_0 = _mm_loadu_pd(v);
              y_0 = _mm_loadu_pd(y);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              i += 2, v += 2, y += 2;
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 8 <= n; i += 8, v += 8, y += (incy * 8)){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              v_2 = _mm_loadu_pd(v + 4);
              v_3 = _mm_loadu_pd(v + 6);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              y_2 = _mm_set_pd(y[(incy * 5)], y[(incy * 4)]);
              y_3 = _mm_set_pd(y[(incy * 7)], y[(incy * 6)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              q_2 = s_0_2;
              q_3 = s_0_3;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              s_0_2 = _mm_add_pd(s_0_2, _mm_or_pd(v_2, mask_BLP));
              s_0_3 = _mm_add_pd(s_0_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              q_2 = _mm_sub_pd(q_2, s_0_2);
              q_3 = _mm_sub_pd(q_3, s_0_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              q_0 = s_1_0;
              q_1 = s_1_1;
              q_2 = s_1_2;
              q_3 = s_1_3;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              s_1_2 = _mm_add_pd(s_1_2, _mm_or_pd(v_2, mask_BLP));
              s_1_3 = _mm_add_pd(s_1_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              q_2 = _mm_sub_pd(q_2, s_1_2);
              q_3 = _mm_sub_pd(q_3, s_1_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              s_2_2 = _mm_add_pd(s_2_2, _mm_or_pd(v_2, mask_BLP));
              s_2_3 = _mm_add_pd(s_2_3, _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              i += 4, v += 4, y += (incy * 4);
            }
            if(i + 2 <= n){
              v_0 = _mm_loadu_pd(v);
              y_0 = _mm_set_pd(y[incy], y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              i += 2, v += 2, y += (incy * 2);
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
            }
          }
        }else{
          if(incy == 1){

            for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += 8){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              v_2 = _mm_set_pd(v[(incv * 5)], v[(incv * 4)]);
              v_3 = _mm_set_pd(v[(incv * 7)], v[(incv * 6)]);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              y_2 = _mm_loadu_pd(y + 4);
              y_3 = _mm_loadu_pd(y + 6);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              q_2 = s_0_2;
              q_3 = s_0_3;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              s_0_2 = _mm_add_pd(s_0_2, _mm_or_pd(v_2, mask_BLP));
              s_0_3 = _mm_add_pd(s_0_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              q_2 = _mm_sub_pd(q_2, s_0_2);
              q_3 = _mm_sub_pd(q_3, s_0_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              q_0 = s_1_0;
              q_1 = s_1_1;
              q_2 = s_1_2;
              q_3 = s_1_3;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              s_1_2 = _mm_add_pd(s_1_2, _mm_or_pd(v_2, mask_BLP));
              s_1_3 = _mm_add_pd(s_1_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              q_2 = _mm_sub_pd(q_2, s_1_2);
              q_3 = _mm_sub_pd(q_3, s_1_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              s_2_2 = _mm_add_pd(s_2_2, _mm_or_pd(v_2, mask_BLP));
              s_2_3 = _mm_add_pd(s_2_3, _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              i += 4, v += (incv * 4), y += 4;
            }
            if(i + 2 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              y_0 = _mm_loadu_pd(y);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              i += 2, v += (incv * 2), y += 2;
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              v_2 = _mm_set_pd(v[(incv * 5)], v[(incv * 4)]);
              v_3 = _mm_set_pd(v[(incv * 7)], v[(incv * 6)]);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              y_2 = _mm_set_pd(y[(incy * 5)], y[(incy * 4)]);
              y_3 = _mm_set_pd(y[(incy * 7)], y[(incy * 6)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              q_0 = s_0_0;
              q_1 = s_0_1;
              q_2 = s_0_2;
              q_3 = s_0_3;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              s_0_2 = _mm_add_pd(s_0_2, _mm_or_pd(v_2, mask_BLP));
              s_0_3 = _mm_add_pd(s_0_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              q_2 = _mm_sub_pd(q_2, s_0_2);
              q_3 = _mm_sub_pd(q_3, s_0_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              q_0 = s_1_0;
              q_1 = s_1_1;
              q_2 = s_1_2;
              q_3 = s_1_3;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              s_1_2 = _mm_add_pd(s_1_2, _mm_or_pd(v_2, mask_BLP));
              s_1_3 = _mm_add_pd(s_1_3, _mm_or_pd(v_3, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              q_2 = _mm_sub_pd(q_2, s_1_2);
              q_3 = _mm_sub_pd(q_3, s_1_3);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              v_2 = _mm_add_pd(v_2, q_2);
              v_3 = _mm_add_pd(v_3, q_3);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              s_2_2 = _mm_add_pd(s_2_2, _mm_or_pd(v_2, mask_BLP));
              s_2_3 = _mm_add_pd(s_2_3, _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              q_0 = s_0_0;
              q_1 = s_0_1;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              q_1 = _mm_sub_pd(q_1, s_0_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              q_0 = s_1_0;
              q_1 = s_1_1;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_1, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              q_1 = _mm_sub_pd(q_1, s_1_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_1, mask_BLP));
              i += 4, v += (incv * 4), y += (incy * 4);
            }
            if(i + 2 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              y_0 = _mm_set_pd(y[incy], y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
              i += 2, v += (incv * 2), y += (incy * 2);
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              q_0 = s_0_0;
              s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_0_0);
              v_0 = _mm_add_pd(v_0, q_0);
              q_0 = s_1_0;
              s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
              q_0 = _mm_sub_pd(q_0, s_1_0);
              v_0 = _mm_add_pd(v_0, q_0);
              s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
            }
          }
        }
        s_0_0 = _mm_sub_pd(s_0_0, _mm_set_pd(sum[0], 0));
        q_0 = _mm_load1_pd(sum);
        s_0_0 = _mm_add_pd(s_0_0, _mm_sub_pd(s_0_1, q_0));
        s_0_0 = _mm_add_pd(s_0_0, _mm_sub_pd(s_0_2, q_0));
        s_0_0 = _mm_add_pd(s_0_0, _mm_sub_pd(s_0_3, q_0));
        _mm_store_pd(tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1];
        s_1_0 = _mm_sub_pd(s_1_0, _mm_set_pd(sum[1], 0));
        q_0 = _mm_load1_pd(sum + 1);
        s_1_0 = _mm_add_pd(s_1_0, _mm_sub_pd(s_1_1, q_0));
        s_1_0 = _mm_add_pd(s_1_0, _mm_sub_pd(s_1_2, q_0));
        s_1_0 = _mm_add_pd(s_1_0, _mm_sub_pd(s_1_3, q_0));
        _mm_store_pd(tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1];
        s_2_0 = _mm_sub_pd(s_2_0, _mm_set_pd(sum[2], 0));
        q_0 = _mm_load1_pd(sum + 2);
        s_2_0 = _mm_add_pd(s_2_0, _mm_sub_pd(s_2_1, q_0));
        s_2_0 = _mm_add_pd(s_2_0, _mm_sub_pd(s_2_2, q_0));
        s_2_0 = _mm_add_pd(s_2_0, _mm_sub_pd(s_2_3, q_0));
        _mm_store_pd(tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        __m128d v_0, v_1, v_2, v_3;
        __m128d y_0, y_1, y_2, y_3;
        __m128d q_0, q_1, q_2, q_3;
        __m128d s_0, s_1, s_2, s_3;
        __m128d s_buffer[(MAX_FOLD * 4)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 3)] = _mm_load1_pd(sum + j);
        }
        if(incv == 1){
          if(incy == 1){

            for(i = 0; i + 8 <= n; i += 8, v += 8, y += 8){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              v_2 = _mm_loadu_pd(v + 4);
              v_3 = _mm_loadu_pd(v + 6);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              y_2 = _mm_loadu_pd(y + 4);
              y_3 = _mm_loadu_pd(y + 6);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                s_2 = s_buffer[((j * 4) + 2)];
                s_3 = s_buffer[((j * 4) + 3)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                q_2 = _mm_add_pd(s_2, _mm_or_pd(v_2, mask_BLP));
                q_3 = _mm_add_pd(s_3, _mm_or_pd(v_3, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                s_buffer[((j * 4) + 2)] = q_2;
                s_buffer[((j * 4) + 3)] = q_3;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                q_2 = _mm_sub_pd(s_2, q_2);
                q_3 = _mm_sub_pd(s_3, q_3);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
                v_2 = _mm_add_pd(v_2, q_2);
                v_3 = _mm_add_pd(v_3, q_3);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              s_buffer[((j * 4) + 2)] = _mm_add_pd(s_buffer[((j * 4) + 2)], _mm_or_pd(v_2, mask_BLP));
              s_buffer[((j * 4) + 3)] = _mm_add_pd(s_buffer[((j * 4) + 3)], _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              i += 4, v += 4, y += 4;
            }
            if(i + 2 <= n){
              v_0 = _mm_loadu_pd(v);
              y_0 = _mm_loadu_pd(y);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              i += 2, v += 2, y += 2;
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 8 <= n; i += 8, v += 8, y += (incy * 8)){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              v_2 = _mm_loadu_pd(v + 4);
              v_3 = _mm_loadu_pd(v + 6);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              y_2 = _mm_set_pd(y[(incy * 5)], y[(incy * 4)]);
              y_3 = _mm_set_pd(y[(incy * 7)], y[(incy * 6)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                s_2 = s_buffer[((j * 4) + 2)];
                s_3 = s_buffer[((j * 4) + 3)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                q_2 = _mm_add_pd(s_2, _mm_or_pd(v_2, mask_BLP));
                q_3 = _mm_add_pd(s_3, _mm_or_pd(v_3, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                s_buffer[((j * 4) + 2)] = q_2;
                s_buffer[((j * 4) + 3)] = q_3;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                q_2 = _mm_sub_pd(s_2, q_2);
                q_3 = _mm_sub_pd(s_3, q_3);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
                v_2 = _mm_add_pd(v_2, q_2);
                v_3 = _mm_add_pd(v_3, q_3);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              s_buffer[((j * 4) + 2)] = _mm_add_pd(s_buffer[((j * 4) + 2)], _mm_or_pd(v_2, mask_BLP));
              s_buffer[((j * 4) + 3)] = _mm_add_pd(s_buffer[((j * 4) + 3)], _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_loadu_pd(v);
              v_1 = _mm_loadu_pd(v + 2);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              i += 4, v += 4, y += (incy * 4);
            }
            if(i + 2 <= n){
              v_0 = _mm_loadu_pd(v);
              y_0 = _mm_set_pd(y[incy], y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              i += 2, v += 2, y += (incy * 2);
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
            }
          }
        }else{
          if(incy == 1){

            for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += 8){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              v_2 = _mm_set_pd(v[(incv * 5)], v[(incv * 4)]);
              v_3 = _mm_set_pd(v[(incv * 7)], v[(incv * 6)]);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              y_2 = _mm_loadu_pd(y + 4);
              y_3 = _mm_loadu_pd(y + 6);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                s_2 = s_buffer[((j * 4) + 2)];
                s_3 = s_buffer[((j * 4) + 3)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                q_2 = _mm_add_pd(s_2, _mm_or_pd(v_2, mask_BLP));
                q_3 = _mm_add_pd(s_3, _mm_or_pd(v_3, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                s_buffer[((j * 4) + 2)] = q_2;
                s_buffer[((j * 4) + 3)] = q_3;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                q_2 = _mm_sub_pd(s_2, q_2);
                q_3 = _mm_sub_pd(s_3, q_3);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
                v_2 = _mm_add_pd(v_2, q_2);
                v_3 = _mm_add_pd(v_3, q_3);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              s_buffer[((j * 4) + 2)] = _mm_add_pd(s_buffer[((j * 4) + 2)], _mm_or_pd(v_2, mask_BLP));
              s_buffer[((j * 4) + 3)] = _mm_add_pd(s_buffer[((j * 4) + 3)], _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              y_0 = _mm_loadu_pd(y);
              y_1 = _mm_loadu_pd(y + 2);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              i += 4, v += (incv * 4), y += 4;
            }
            if(i + 2 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              y_0 = _mm_loadu_pd(y);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              i += 2, v += (incv * 2), y += 2;
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
            }
          }else{

            for(i = 0; i + 8 <= n; i += 8, v += (incv * 8), y += (incy * 8)){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              v_2 = _mm_set_pd(v[(incv * 5)], v[(incv * 4)]);
              v_3 = _mm_set_pd(v[(incv * 7)], v[(incv * 6)]);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              y_2 = _mm_set_pd(y[(incy * 5)], y[(incy * 4)]);
              y_3 = _mm_set_pd(y[(incy * 7)], y[(incy * 6)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              v_2 = _mm_mul_pd(v_2, y_2);
              v_3 = _mm_mul_pd(v_3, y_3);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                s_2 = s_buffer[((j * 4) + 2)];
                s_3 = s_buffer[((j * 4) + 3)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                q_2 = _mm_add_pd(s_2, _mm_or_pd(v_2, mask_BLP));
                q_3 = _mm_add_pd(s_3, _mm_or_pd(v_3, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                s_buffer[((j * 4) + 2)] = q_2;
                s_buffer[((j * 4) + 3)] = q_3;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                q_2 = _mm_sub_pd(s_2, q_2);
                q_3 = _mm_sub_pd(s_3, q_3);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
                v_2 = _mm_add_pd(v_2, q_2);
                v_3 = _mm_add_pd(v_3, q_3);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              s_buffer[((j * 4) + 2)] = _mm_add_pd(s_buffer[((j * 4) + 2)], _mm_or_pd(v_2, mask_BLP));
              s_buffer[((j * 4) + 3)] = _mm_add_pd(s_buffer[((j * 4) + 3)], _mm_or_pd(v_3, mask_BLP));
            }
            if(i + 4 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              v_1 = _mm_set_pd(v[(incv * 3)], v[(incv * 2)]);
              y_0 = _mm_set_pd(y[incy], y[0]);
              y_1 = _mm_set_pd(y[(incy * 3)], y[(incy * 2)]);
              v_0 = _mm_mul_pd(v_0, y_0);
              v_1 = _mm_mul_pd(v_1, y_1);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                s_1 = s_buffer[((j * 4) + 1)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                s_buffer[((j * 4) + 1)] = q_1;
                q_0 = _mm_sub_pd(s_0, q_0);
                q_1 = _mm_sub_pd(s_1, q_1);
                v_0 = _mm_add_pd(v_0, q_0);
                v_1 = _mm_add_pd(v_1, q_1);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              s_buffer[((j * 4) + 1)] = _mm_add_pd(s_buffer[((j * 4) + 1)], _mm_or_pd(v_1, mask_BLP));
              i += 4, v += (incv * 4), y += (incy * 4);
            }
            if(i + 2 <= n){
              v_0 = _mm_set_pd(v[incv], v[0]);
              y_0 = _mm_set_pd(y[incy], y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
              i += 2, v += (incv * 2), y += (incy * 2);
            }
            if(i < n){
              v_0 = _mm_set_pd(0, v[0]);
              y_0 = _mm_set_pd(0, y[0]);
              v_0 = _mm_mul_pd(v_0, y_0);
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 4)];
                q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
                s_buffer[(j * 4)] = q_0;
                q_0 = _mm_sub_pd(s_0, q_0);
                v_0 = _mm_add_pd(v_0, q_0);
              }
              s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_or_pd(v_0, mask_BLP));
            }
          }
        }
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = _mm_sub_pd(s_buffer[(j * 4)], _mm_set_pd(sum[j], 0));
          q_0 = _mm_load1_pd(sum + j);
          s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_sub_pd(s_buffer[((j * 4) + 1)], q_0));
          s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_sub_pd(s_buffer[((j * 4) + 2)], q_0));
          s_buffer[(j * 4)] = _mm_add_pd(s_buffer[(j * 4)], _mm_sub_pd(s_buffer[((j * 4) + 3)], q_0));
          _mm_store_pd(tmp_cons, s_buffer[(j * 4)]);
          sum[j] = tmp_cons[0] + tmp_cons[1];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
  }
#else
  void ddotI2(int n, double* v, int incv, double* y, int incy, int fold, double* sum){
    long_double tmp_BLP;
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        double v_0, v_1;
        double y_0, y_1;
        double q_0, q_1;
        double s_0_0, s_0_1;
        double s_1_0, s_1_1;
        double s_2_0, s_2_1;

        s_0_0 = s_0_1 = sum[0];
        s_1_0 = s_1_1 = sum[1];
        s_2_0 = s_2_1 = sum[2];
        if(incv == 1){
          if(incy == 1){

            for(i = 0; i + 2 <= n; i += 2, v += 2, y += 2){
              v_0 = v[0];
              v_1 = v[1];
              y_0 = y[0];
              y_1 = y[1];
              v_0 = v_0 * y_0;
              v_1 = v_1 * y_1;
              q_0 = s_0_0;
              q_1 = s_0_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_0_1 = s_0_1 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              q_1 = q_1 - s_0_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              q_0 = s_1_0;
              q_1 = s_1_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_1_1 = s_1_1 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              q_1 = q_1 - s_1_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_2_1 = s_2_1 + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              q_0 = s_0_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              v_0 = v_0 + q_0;
              q_0 = s_1_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              v_0 = v_0 + q_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              i += 1, v += 1, y += 1;
            }
          }else{

            for(i = 0; i + 2 <= n; i += 2, v += 2, y += (incy * 2)){
              v_0 = v[0];
              v_1 = v[1];
              y_0 = y[0];
              y_1 = y[incy];
              v_0 = v_0 * y_0;
              v_1 = v_1 * y_1;
              q_0 = s_0_0;
              q_1 = s_0_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_0_1 = s_0_1 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              q_1 = q_1 - s_0_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              q_0 = s_1_0;
              q_1 = s_1_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_1_1 = s_1_1 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              q_1 = q_1 - s_1_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_2_1 = s_2_1 + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              q_0 = s_0_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              v_0 = v_0 + q_0;
              q_0 = s_1_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              v_0 = v_0 + q_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              i += 1, v += 1, y += incy;
            }
          }
        }else{
          if(incy == 1){

            for(i = 0; i + 2 <= n; i += 2, v += (incv * 2), y += 2){
              v_0 = v[0];
              v_1 = v[incv];
              y_0 = y[0];
              y_1 = y[1];
              v_0 = v_0 * y_0;
              v_1 = v_1 * y_1;
              q_0 = s_0_0;
              q_1 = s_0_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_0_1 = s_0_1 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              q_1 = q_1 - s_0_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              q_0 = s_1_0;
              q_1 = s_1_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_1_1 = s_1_1 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              q_1 = q_1 - s_1_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_2_1 = s_2_1 + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              q_0 = s_0_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              v_0 = v_0 + q_0;
              q_0 = s_1_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              v_0 = v_0 + q_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              i += 1, v += incv, y += 1;
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
              q_1 = s_0_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_0_1 = s_0_1 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              q_1 = q_1 - s_0_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              q_0 = s_1_0;
              q_1 = s_1_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_1_1 = s_1_1 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              q_1 = q_1 - s_1_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_2_1 = s_2_1 + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              q_0 = s_0_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_0_0 = s_0_0 + tmp_BLP.d;
              q_0 = q_0 - s_0_0;
              v_0 = v_0 + q_0;
              q_0 = s_1_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_1_0 = s_1_0 + tmp_BLP.d;
              q_0 = q_0 - s_1_0;
              v_0 = v_0 + q_0;
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_2_0 = s_2_0 + tmp_BLP.d;
              i += 1, v += incv, y += incy;
            }
          }
        }
        q_0 = sum[0];
        s_0_0 = s_0_0 + (s_0_1 - q_0);
        sum[0] = s_0_0;
        q_0 = sum[1];
        s_1_0 = s_1_0 + (s_1_1 - q_0);
        sum[1] = s_1_0;
        q_0 = sum[2];
        s_2_0 = s_2_0 + (s_2_1 - q_0);
        sum[2] = s_2_0;
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        double v_0, v_1;
        double y_0, y_1;
        double q_0, q_1;
        double s_0, s_1;
        double s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = sum[j];
        }
        if(incv == 1){
          if(incy == 1){

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
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                tmp_BLP.d = v_1;
                tmp_BLP.l |= 1;
                q_1 = s_1 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = s_0 - q_0;
                q_1 = s_1 - q_1;
                v_0 = v_0 + q_0;
                v_1 = v_1 + q_1;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_buffer[((j * 2) + 1)] = s_buffer[((j * 2) + 1)] + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                q_0 = s_0 - q_0;
                v_0 = v_0 + q_0;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              i += 1, v += 1, y += 1;
            }
          }else{

            for(i = 0; i + 2 <= n; i += 2, v += 2, y += (incy * 2)){
              v_0 = v[0];
              v_1 = v[1];
              y_0 = y[0];
              y_1 = y[incy];
              v_0 = v_0 * y_0;
              v_1 = v_1 * y_1;
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                s_1 = s_buffer[((j * 2) + 1)];
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                tmp_BLP.d = v_1;
                tmp_BLP.l |= 1;
                q_1 = s_1 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = s_0 - q_0;
                q_1 = s_1 - q_1;
                v_0 = v_0 + q_0;
                v_1 = v_1 + q_1;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_buffer[((j * 2) + 1)] = s_buffer[((j * 2) + 1)] + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                q_0 = s_0 - q_0;
                v_0 = v_0 + q_0;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              i += 1, v += 1, y += incy;
            }
          }
        }else{
          if(incy == 1){

            for(i = 0; i + 2 <= n; i += 2, v += (incv * 2), y += 2){
              v_0 = v[0];
              v_1 = v[incv];
              y_0 = y[0];
              y_1 = y[1];
              v_0 = v_0 * y_0;
              v_1 = v_1 * y_1;
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                s_1 = s_buffer[((j * 2) + 1)];
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                tmp_BLP.d = v_1;
                tmp_BLP.l |= 1;
                q_1 = s_1 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = s_0 - q_0;
                q_1 = s_1 - q_1;
                v_0 = v_0 + q_0;
                v_1 = v_1 + q_1;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_buffer[((j * 2) + 1)] = s_buffer[((j * 2) + 1)] + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                q_0 = s_0 - q_0;
                v_0 = v_0 + q_0;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              i += 1, v += incv, y += 1;
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
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                tmp_BLP.d = v_1;
                tmp_BLP.l |= 1;
                q_1 = s_1 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                s_buffer[((j * 2) + 1)] = q_1;
                q_0 = s_0 - q_0;
                q_1 = s_1 - q_1;
                v_0 = v_0 + q_0;
                v_1 = v_1 + q_1;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              s_buffer[((j * 2) + 1)] = s_buffer[((j * 2) + 1)] + tmp_BLP.d;
            }
            if(i + 1 <= n){
              v_0 = v[0];
              y_0 = y[0];
              v_0 = v_0 * y_0;
              for(j = 0; j < fold - 1; j++){
                s_0 = s_buffer[(j * 2)];
                tmp_BLP.d = v_0;
                tmp_BLP.l |= 1;
                q_0 = s_0 + tmp_BLP.d;
                s_buffer[(j * 2)] = q_0;
                q_0 = s_0 - q_0;
                v_0 = v_0 + q_0;
              }
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              s_buffer[(j * 2)] = s_buffer[(j * 2)] + tmp_BLP.d;
              i += 1, v += incv, y += incy;
            }
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