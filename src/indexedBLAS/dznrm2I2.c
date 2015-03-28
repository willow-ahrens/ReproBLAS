#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "../config.h"
#include "../Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>

/*[[[cog
import cog
import sys, os
import generate
import dataTypes
import vectorizations
import nrm2I2
]]]*/
//[[[end]]]

#if defined( __AVX__ )
  void dznrm2I2(int n, double complex* v, int incv, double scale, int fold, double complex* sum){
    /*[[[cog
    cog.out(generate.generate(nrm2I2.Nrm2I2(dataTypes.DoubleComplex, vectorizations.AVX), args, params))
    ]]]*/
    __m256d scale_mask = _mm256_set1_pd(scale);
    __m256d mask_BLP; AVX_BLP_MASKD(mask_BLP);
    double complex tmp_cons[2] __attribute__((aligned(32)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        double* sum_base = (double*) sum;
        double* v_base = (double*) v;
        __m256d v_0, v_1, v_2, v_3;
        __m256d q_0, q_1;
        __m256d s_0_0, s_0_1;
        __m256d s_1_0, s_1_1;
        __m256d s_2_0, s_2_1;

        s_0_0 = s_0_1 = _mm256_broadcast_pd((__m128d *)(sum_base));
        s_1_0 = s_1_1 = _mm256_broadcast_pd((__m128d *)(sum_base + 2));
        s_2_0 = s_2_1 = _mm256_broadcast_pd((__m128d *)(sum_base + 4));
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16){
            v_0 = _mm256_mul_pd(_mm256_loadu_pd(v_base), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 4), scale_mask);
            v_2 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 8), scale_mask);
            v_3 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 12), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
            v_2 = _mm256_mul_pd(v_2, v_2);
            v_3 = _mm256_mul_pd(v_3, v_3);
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
          }
          if(i + 4 <= n){
            v_0 = _mm256_mul_pd(_mm256_loadu_pd(v_base), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 4), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
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
            i += 4, v_base += 8;
          }
          if(i + 2 <= n){
            v_0 = _mm256_mul_pd(_mm256_loadu_pd(v_base), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
            q_0 = _mm256_sub_pd(q_0, s_0_0);
            v_0 = _mm256_add_pd(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
            q_0 = _mm256_sub_pd(q_0, s_1_0);
            v_0 = _mm256_add_pd(v_0, q_0);
            s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
            i += 2, v_base += 4;
          }
          if(i < n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(0, 0, v_base[1], v_base[0]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
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

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
            v_0 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), scale_mask);
            v_2 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), scale_mask);
            v_3 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
            v_2 = _mm256_mul_pd(v_2, v_2);
            v_3 = _mm256_mul_pd(v_3, v_3);
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
          }
          if(i + 4 <= n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
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
            i += 4, v_base += (incv * 8);
          }
          if(i + 2 <= n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            q_0 = s_0_0;
            s_0_0 = _mm256_add_pd(s_0_0, _mm256_or_pd(v_0, mask_BLP));
            q_0 = _mm256_sub_pd(q_0, s_0_0);
            v_0 = _mm256_add_pd(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm256_add_pd(s_1_0, _mm256_or_pd(v_0, mask_BLP));
            q_0 = _mm256_sub_pd(q_0, s_1_0);
            v_0 = _mm256_add_pd(v_0, q_0);
            s_2_0 = _mm256_add_pd(s_2_0, _mm256_or_pd(v_0, mask_BLP));
            i += 2, v_base += (incv * 4);
          }
          if(i < n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(0, 0, v_base[1], v_base[0]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
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
        s_0_0 = _mm256_sub_pd(s_0_0, _mm256_set_pd(sum_base[1], sum_base[0], 0, 0));
        q_0 = _mm256_broadcast_pd((__m128d *)(sum_base));
        s_0_0 = _mm256_add_pd(s_0_0, _mm256_sub_pd(s_0_1, q_0));
        _mm256_store_pd((double*)tmp_cons, s_0_0);
        sum[0] = tmp_cons[0] + tmp_cons[1];
        s_1_0 = _mm256_sub_pd(s_1_0, _mm256_set_pd(sum_base[3], sum_base[2], 0, 0));
        q_0 = _mm256_broadcast_pd((__m128d *)(sum_base + 2));
        s_1_0 = _mm256_add_pd(s_1_0, _mm256_sub_pd(s_1_1, q_0));
        _mm256_store_pd((double*)tmp_cons, s_1_0);
        sum[1] = tmp_cons[0] + tmp_cons[1];
        s_2_0 = _mm256_sub_pd(s_2_0, _mm256_set_pd(sum_base[5], sum_base[4], 0, 0));
        q_0 = _mm256_broadcast_pd((__m128d *)(sum_base + 4));
        s_2_0 = _mm256_add_pd(s_2_0, _mm256_sub_pd(s_2_1, q_0));
        _mm256_store_pd((double*)tmp_cons, s_2_0);
        sum[2] = tmp_cons[0] + tmp_cons[1];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        double* sum_base = (double*) sum;
        double* v_base = (double*) v;
        __m256d v_0, v_1, v_2, v_3;
        __m256d q_0, q_1;
        __m256d s_0, s_1;
        __m256d s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = _mm256_broadcast_pd((__m128d *)(sum_base + (j * 2)));
        }
        if(incv == 1){

          for(i = 0; i + 8 <= n; i += 8, v_base += 16){
            v_0 = _mm256_mul_pd(_mm256_loadu_pd(v_base), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 4), scale_mask);
            v_2 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 8), scale_mask);
            v_3 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 12), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
            v_2 = _mm256_mul_pd(v_2, v_2);
            v_3 = _mm256_mul_pd(v_3, v_3);
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
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_2, mask_BLP));
              q_1 = _mm256_add_pd(s_1, _mm256_or_pd(v_3, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm256_sub_pd(s_0, q_0);
              q_1 = _mm256_sub_pd(s_1, q_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
            }
            s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_2, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm256_add_pd(s_buffer[((j * 2) + 1)], _mm256_or_pd(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm256_mul_pd(_mm256_loadu_pd(v_base), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_loadu_pd(v_base + 4), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
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
            i += 4, v_base += 8;
          }
          if(i + 2 <= n){
            v_0 = _mm256_mul_pd(_mm256_loadu_pd(v_base), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm256_sub_pd(s_0, q_0);
              v_0 = _mm256_add_pd(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
            i += 2, v_base += 4;
          }
          if(i < n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(0, 0, v_base[1], v_base[0]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
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

          for(i = 0; i + 8 <= n; i += 8, v_base += (incv * 16)){
            v_0 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), scale_mask);
            v_2 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 10) + 1)], v_base[(incv * 10)], v_base[((incv * 8) + 1)], v_base[(incv * 8)]), scale_mask);
            v_3 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 14) + 1)], v_base[(incv * 14)], v_base[((incv * 12) + 1)], v_base[(incv * 12)]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
            v_2 = _mm256_mul_pd(v_2, v_2);
            v_3 = _mm256_mul_pd(v_3, v_3);
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
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_2, mask_BLP));
              q_1 = _mm256_add_pd(s_1, _mm256_or_pd(v_3, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm256_sub_pd(s_0, q_0);
              q_1 = _mm256_sub_pd(s_1, q_1);
              v_2 = _mm256_add_pd(v_2, q_0);
              v_3 = _mm256_add_pd(v_3, q_1);
            }
            s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_2, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm256_add_pd(s_buffer[((j * 2) + 1)], _mm256_or_pd(v_3, mask_BLP));
          }
          if(i + 4 <= n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), scale_mask);
            v_1 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 6) + 1)], v_base[(incv * 6)], v_base[((incv * 4) + 1)], v_base[(incv * 4)]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            v_1 = _mm256_mul_pd(v_1, v_1);
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
            i += 4, v_base += (incv * 8);
          }
          if(i + 2 <= n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(v_base[((incv * 2) + 1)], v_base[(incv * 2)], v_base[1], v_base[0]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm256_add_pd(s_0, _mm256_or_pd(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm256_sub_pd(s_0, q_0);
              v_0 = _mm256_add_pd(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_or_pd(v_0, mask_BLP));
            i += 2, v_base += (incv * 4);
          }
          if(i < n){
            v_0 = _mm256_mul_pd(_mm256_set_pd(0, 0, v_base[1], v_base[0]), scale_mask);
            v_0 = _mm256_mul_pd(v_0, v_0);
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
        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = _mm256_sub_pd(s_buffer[(j * 2)], _mm256_set_pd(sum_base[((j * 2) + 1)], sum_base[(j * 2)], 0, 0));
          q_0 = _mm256_broadcast_pd((__m128d *)(sum_base + (j * 2)));
          s_buffer[(j * 2)] = _mm256_add_pd(s_buffer[(j * 2)], _mm256_sub_pd(s_buffer[((j * 2) + 1)], q_0));
          _mm256_store_pd((double*)tmp_cons, s_buffer[(j * 2)]);
          sum[j] = tmp_cons[0] + tmp_cons[1];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
    //[[[end]]]
  }
#elif defined( __SSE2__ )
  void dznrm2I2(int n, double complex* v, int incv, double scale, int fold, double complex* sum){
    /*[[[cog
    cog.out(generate.generate(nrm2I2.Nrm2I2(dataTypes.DoubleComplex, vectorizations.SSE), args, params))
    ]]]*/
    __m128d scale_mask = _mm_set1_pd(scale);
    __m128d mask_BLP; SSE_BLP_MASKD(mask_BLP);
    double complex tmp_cons[1] __attribute__((aligned(16)));
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        double* sum_base = (double*) sum;
        double* v_base = (double*) v;
        __m128d v_0, v_1, v_2, v_3;
        __m128d q_0, q_1;
        __m128d s_0_0, s_0_1;
        __m128d s_1_0, s_1_1;
        __m128d s_2_0, s_2_1;

        s_0_0 = s_0_1 = _mm_loadu_pd(sum_base);
        s_1_0 = s_1_1 = _mm_loadu_pd(sum_base + 2);
        s_2_0 = s_2_1 = _mm_loadu_pd(sum_base + 4);
        if(incv == 1){

          for(i = 0; i + 4 <= n; i += 4, v_base += 8){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + 2), scale_mask);
            v_2 = _mm_mul_pd(_mm_loadu_pd(v_base + 4), scale_mask);
            v_3 = _mm_mul_pd(_mm_loadu_pd(v_base + 6), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
            v_2 = _mm_mul_pd(v_2, v_2);
            v_3 = _mm_mul_pd(v_3, v_3);
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
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_2, mask_BLP));
            s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_3, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_0_0);
            q_1 = _mm_sub_pd(q_1, s_0_1);
            v_2 = _mm_add_pd(v_2, q_0);
            v_3 = _mm_add_pd(v_3, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_2, mask_BLP));
            s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_3, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_1_0);
            q_1 = _mm_sub_pd(q_1, s_1_1);
            v_2 = _mm_add_pd(v_2, q_0);
            v_3 = _mm_add_pd(v_3, q_1);
            s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_2, mask_BLP));
            s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_3, mask_BLP));
          }
          if(i + 2 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + 2), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
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
            i += 2, v_base += 4;
          }
          if(i + 1 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            q_0 = s_0_0;
            s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_0_0);
            v_0 = _mm_add_pd(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_1_0);
            v_0 = _mm_add_pd(v_0, q_0);
            s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
            i += 1, v_base += 2;
          }
        }else{

          for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 2)), scale_mask);
            v_2 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 4)), scale_mask);
            v_3 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 6)), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
            v_2 = _mm_mul_pd(v_2, v_2);
            v_3 = _mm_mul_pd(v_3, v_3);
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
            q_0 = s_0_0;
            q_1 = s_0_1;
            s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_2, mask_BLP));
            s_0_1 = _mm_add_pd(s_0_1, _mm_or_pd(v_3, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_0_0);
            q_1 = _mm_sub_pd(q_1, s_0_1);
            v_2 = _mm_add_pd(v_2, q_0);
            v_3 = _mm_add_pd(v_3, q_1);
            q_0 = s_1_0;
            q_1 = s_1_1;
            s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_2, mask_BLP));
            s_1_1 = _mm_add_pd(s_1_1, _mm_or_pd(v_3, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_1_0);
            q_1 = _mm_sub_pd(q_1, s_1_1);
            v_2 = _mm_add_pd(v_2, q_0);
            v_3 = _mm_add_pd(v_3, q_1);
            s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_2, mask_BLP));
            s_2_1 = _mm_add_pd(s_2_1, _mm_or_pd(v_3, mask_BLP));
          }
          if(i + 2 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 2)), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
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
            i += 2, v_base += (incv * 4);
          }
          if(i + 1 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            q_0 = s_0_0;
            s_0_0 = _mm_add_pd(s_0_0, _mm_or_pd(v_0, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_0_0);
            v_0 = _mm_add_pd(v_0, q_0);
            q_0 = s_1_0;
            s_1_0 = _mm_add_pd(s_1_0, _mm_or_pd(v_0, mask_BLP));
            q_0 = _mm_sub_pd(q_0, s_1_0);
            v_0 = _mm_add_pd(v_0, q_0);
            s_2_0 = _mm_add_pd(s_2_0, _mm_or_pd(v_0, mask_BLP));
            i += 1, v_base += (incv * 2);
          }
        }
        q_0 = _mm_loadu_pd(sum_base);
        s_0_0 = _mm_add_pd(s_0_0, _mm_sub_pd(s_0_1, q_0));
        _mm_store_pd((double*)tmp_cons, s_0_0);
        sum[0] = tmp_cons[0];
        q_0 = _mm_loadu_pd(sum_base + 2);
        s_1_0 = _mm_add_pd(s_1_0, _mm_sub_pd(s_1_1, q_0));
        _mm_store_pd((double*)tmp_cons, s_1_0);
        sum[1] = tmp_cons[0];
        q_0 = _mm_loadu_pd(sum_base + 4);
        s_2_0 = _mm_add_pd(s_2_0, _mm_sub_pd(s_2_1, q_0));
        _mm_store_pd((double*)tmp_cons, s_2_0);
        sum[2] = tmp_cons[0];
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        double* sum_base = (double*) sum;
        double* v_base = (double*) v;
        __m128d v_0, v_1, v_2, v_3;
        __m128d q_0, q_1;
        __m128d s_0, s_1;
        __m128d s_buffer[(MAX_FOLD * 2)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 2)] = s_buffer[((j * 2) + 1)] = _mm_loadu_pd(sum_base + (j * 2));
        }
        if(incv == 1){

          for(i = 0; i + 4 <= n; i += 4, v_base += 8){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + 2), scale_mask);
            v_2 = _mm_mul_pd(_mm_loadu_pd(v_base + 4), scale_mask);
            v_3 = _mm_mul_pd(_mm_loadu_pd(v_base + 6), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
            v_2 = _mm_mul_pd(v_2, v_2);
            v_3 = _mm_mul_pd(v_3, v_3);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
              q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_pd(s_0, q_0);
              q_1 = _mm_sub_pd(s_1, q_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_2, mask_BLP));
              q_1 = _mm_add_pd(s_1, _mm_or_pd(v_3, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_pd(s_0, q_0);
              q_1 = _mm_sub_pd(s_1, q_1);
              v_2 = _mm_add_pd(v_2, q_0);
              v_3 = _mm_add_pd(v_3, q_1);
            }
            s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_or_pd(v_2, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm_add_pd(s_buffer[((j * 2) + 1)], _mm_or_pd(v_3, mask_BLP));
          }
          if(i + 2 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + 2), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
              q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_pd(s_0, q_0);
              q_1 = _mm_sub_pd(s_1, q_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
            }
            s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_or_pd(v_0, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm_add_pd(s_buffer[((j * 2) + 1)], _mm_or_pd(v_1, mask_BLP));
            i += 2, v_base += 4;
          }
          if(i + 1 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_pd(s_0, q_0);
              v_0 = _mm_add_pd(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_or_pd(v_0, mask_BLP));
            i += 1, v_base += 2;
          }
        }else{

          for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 2)), scale_mask);
            v_2 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 4)), scale_mask);
            v_3 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 6)), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
            v_2 = _mm_mul_pd(v_2, v_2);
            v_3 = _mm_mul_pd(v_3, v_3);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
              q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_pd(s_0, q_0);
              q_1 = _mm_sub_pd(s_1, q_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_2, mask_BLP));
              q_1 = _mm_add_pd(s_1, _mm_or_pd(v_3, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_pd(s_0, q_0);
              q_1 = _mm_sub_pd(s_1, q_1);
              v_2 = _mm_add_pd(v_2, q_0);
              v_3 = _mm_add_pd(v_3, q_1);
            }
            s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_or_pd(v_2, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm_add_pd(s_buffer[((j * 2) + 1)], _mm_or_pd(v_3, mask_BLP));
          }
          if(i + 2 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_1 = _mm_mul_pd(_mm_loadu_pd(v_base + (incv * 2)), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            v_1 = _mm_mul_pd(v_1, v_1);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              s_1 = s_buffer[((j * 2) + 1)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
              q_1 = _mm_add_pd(s_1, _mm_or_pd(v_1, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              s_buffer[((j * 2) + 1)] = q_1;
              q_0 = _mm_sub_pd(s_0, q_0);
              q_1 = _mm_sub_pd(s_1, q_1);
              v_0 = _mm_add_pd(v_0, q_0);
              v_1 = _mm_add_pd(v_1, q_1);
            }
            s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_or_pd(v_0, mask_BLP));
            s_buffer[((j * 2) + 1)] = _mm_add_pd(s_buffer[((j * 2) + 1)], _mm_or_pd(v_1, mask_BLP));
            i += 2, v_base += (incv * 4);
          }
          if(i + 1 <= n){
            v_0 = _mm_mul_pd(_mm_loadu_pd(v_base), scale_mask);
            v_0 = _mm_mul_pd(v_0, v_0);
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 2)];
              q_0 = _mm_add_pd(s_0, _mm_or_pd(v_0, mask_BLP));
              s_buffer[(j * 2)] = q_0;
              q_0 = _mm_sub_pd(s_0, q_0);
              v_0 = _mm_add_pd(v_0, q_0);
            }
            s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_or_pd(v_0, mask_BLP));
            i += 1, v_base += (incv * 2);
          }
        }
        for(j = 0; j < fold; j += 1){
          q_0 = _mm_loadu_pd(sum_base + (j * 2));
          s_buffer[(j * 2)] = _mm_add_pd(s_buffer[(j * 2)], _mm_sub_pd(s_buffer[((j * 2) + 1)], q_0));
          _mm_store_pd((double*)tmp_cons, s_buffer[(j * 2)]);
          sum[j] = tmp_cons[0];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
    //[[[end]]]
  }
#else
  void dznrm2I2(int n, double complex* v, int incv, double scale, int fold, double complex* sum){
    /*[[[cog
    cog.out(generate.generate(nrm2I2.Nrm2I2(dataTypes.DoubleComplex, vectorizations.SISD), args, params))
    ]]]*/
    double scale_mask = scale;
    long_double tmp_BLP;
    SET_DAZ_FLAG;
    switch(fold){
      case 3:{
        int i;

        double* sum_base = (double*) sum;
        double* v_base = (double*) v;
        double v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        double q_0, q_1, q_2, q_3;
        double s_0_0, s_0_1, s_0_2, s_0_3;
        double s_1_0, s_1_1, s_1_2, s_1_3;
        double s_2_0, s_2_1, s_2_2, s_2_3;

        s_0_0 = s_0_2 = sum_base[0];
        s_0_1 = s_0_3 = sum_base[1];
        s_1_0 = s_1_2 = sum_base[2];
        s_1_1 = s_1_3 = sum_base[3];
        s_2_0 = s_2_2 = sum_base[4];
        s_2_1 = s_2_3 = sum_base[5];
        if(incv == 1){

          for(i = 0; i + 4 <= n; i += 4, v_base += 8){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[2] * scale_mask;
            v_3 = v_base[3] * scale_mask;
            v_4 = v_base[4] * scale_mask;
            v_5 = v_base[5] * scale_mask;
            v_6 = v_base[6] * scale_mask;
            v_7 = v_base[7] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            v_4 = v_4 * v_4;
            v_5 = v_5 * v_5;
            v_6 = v_6 * v_6;
            v_7 = v_7 * v_7;
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_0_0 = s_0_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_0_1 = s_0_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_0_2 = s_0_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_0_3 = s_0_3 + tmp_BLP.d;
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
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_1_0 = s_1_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_1_1 = s_1_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_1_2 = s_1_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_1_3 = s_1_3 + tmp_BLP.d;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_2_0 = s_2_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_2_1 = s_2_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_2_2 = s_2_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_2_3 = s_2_3 + tmp_BLP.d;
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_0_0 = s_0_0 + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_0_1 = s_0_1 + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_0_2 = s_0_2 + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_0_3 = s_0_3 + tmp_BLP.d;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            q_2 = q_2 - s_0_2;
            q_3 = q_3 - s_0_3;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            v_6 = v_6 + q_2;
            v_7 = v_7 + q_3;
            q_0 = s_1_0;
            q_1 = s_1_1;
            q_2 = s_1_2;
            q_3 = s_1_3;
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_1_0 = s_1_0 + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_1_1 = s_1_1 + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_1_2 = s_1_2 + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_1_3 = s_1_3 + tmp_BLP.d;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            v_6 = v_6 + q_2;
            v_7 = v_7 + q_3;
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_2_0 = s_2_0 + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_2_1 = s_2_1 + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_2_2 = s_2_2 + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_2_3 = s_2_3 + tmp_BLP.d;
          }
          if(i + 2 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[2] * scale_mask;
            v_3 = v_base[3] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_0_0 = s_0_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_0_1 = s_0_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_0_2 = s_0_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_0_3 = s_0_3 + tmp_BLP.d;
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
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_1_0 = s_1_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_1_1 = s_1_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_1_2 = s_1_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_1_3 = s_1_3 + tmp_BLP.d;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_2_0 = s_2_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_2_1 = s_2_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_2_2 = s_2_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_2_3 = s_2_3 + tmp_BLP.d;
            i += 2, v_base += 4;
          }
          if(i + 1 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
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
            i += 1, v_base += 2;
          }
        }else{

          for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[(incv * 2)] * scale_mask;
            v_3 = v_base[((incv * 2) + 1)] * scale_mask;
            v_4 = v_base[(incv * 4)] * scale_mask;
            v_5 = v_base[((incv * 4) + 1)] * scale_mask;
            v_6 = v_base[(incv * 6)] * scale_mask;
            v_7 = v_base[((incv * 6) + 1)] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            v_4 = v_4 * v_4;
            v_5 = v_5 * v_5;
            v_6 = v_6 * v_6;
            v_7 = v_7 * v_7;
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_0_0 = s_0_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_0_1 = s_0_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_0_2 = s_0_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_0_3 = s_0_3 + tmp_BLP.d;
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
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_1_0 = s_1_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_1_1 = s_1_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_1_2 = s_1_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_1_3 = s_1_3 + tmp_BLP.d;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_2_0 = s_2_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_2_1 = s_2_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_2_2 = s_2_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_2_3 = s_2_3 + tmp_BLP.d;
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_0_0 = s_0_0 + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_0_1 = s_0_1 + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_0_2 = s_0_2 + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_0_3 = s_0_3 + tmp_BLP.d;
            q_0 = q_0 - s_0_0;
            q_1 = q_1 - s_0_1;
            q_2 = q_2 - s_0_2;
            q_3 = q_3 - s_0_3;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            v_6 = v_6 + q_2;
            v_7 = v_7 + q_3;
            q_0 = s_1_0;
            q_1 = s_1_1;
            q_2 = s_1_2;
            q_3 = s_1_3;
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_1_0 = s_1_0 + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_1_1 = s_1_1 + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_1_2 = s_1_2 + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_1_3 = s_1_3 + tmp_BLP.d;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_4 = v_4 + q_0;
            v_5 = v_5 + q_1;
            v_6 = v_6 + q_2;
            v_7 = v_7 + q_3;
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_2_0 = s_2_0 + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_2_1 = s_2_1 + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_2_2 = s_2_2 + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_2_3 = s_2_3 + tmp_BLP.d;
          }
          if(i + 2 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[(incv * 2)] * scale_mask;
            v_3 = v_base[((incv * 2) + 1)] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            q_0 = s_0_0;
            q_1 = s_0_1;
            q_2 = s_0_2;
            q_3 = s_0_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_0_0 = s_0_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_0_1 = s_0_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_0_2 = s_0_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_0_3 = s_0_3 + tmp_BLP.d;
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
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_1_0 = s_1_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_1_1 = s_1_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_1_2 = s_1_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_1_3 = s_1_3 + tmp_BLP.d;
            q_0 = q_0 - s_1_0;
            q_1 = q_1 - s_1_1;
            q_2 = q_2 - s_1_2;
            q_3 = q_3 - s_1_3;
            v_0 = v_0 + q_0;
            v_1 = v_1 + q_1;
            v_2 = v_2 + q_2;
            v_3 = v_3 + q_3;
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_2_0 = s_2_0 + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_2_1 = s_2_1 + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_2_2 = s_2_2 + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_2_3 = s_2_3 + tmp_BLP.d;
            i += 2, v_base += (incv * 4);
          }
          if(i + 1 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
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
            i += 1, v_base += (incv * 2);
          }
        }
        q_0 = ((double*)sum)[0];
        s_0_0 = s_0_0 + (s_0_2 - q_0);
        q_0 = ((double*)sum)[1];
        s_0_1 = s_0_1 + (s_0_3 - q_0);
        ((double*)sum)[0] = s_0_0;
        ((double*)sum)[1] = s_0_1;
        q_0 = ((double*)sum)[2];
        s_1_0 = s_1_0 + (s_1_2 - q_0);
        q_0 = ((double*)sum)[3];
        s_1_1 = s_1_1 + (s_1_3 - q_0);
        ((double*)sum)[2] = s_1_0;
        ((double*)sum)[3] = s_1_1;
        q_0 = ((double*)sum)[4];
        s_2_0 = s_2_0 + (s_2_2 - q_0);
        q_0 = ((double*)sum)[5];
        s_2_1 = s_2_1 + (s_2_3 - q_0);
        ((double*)sum)[4] = s_2_0;
        ((double*)sum)[5] = s_2_1;
        RESET_DAZ_FLAG
        return;
      }
      default:{
        int i, j;

        double* sum_base = (double*) sum;
        double* v_base = (double*) v;
        double v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
        double q_0, q_1, q_2, q_3;
        double s_0, s_1, s_2, s_3;
        double s_buffer[(MAX_FOLD * 4)];

        for(j = 0; j < fold; j += 1){
          s_buffer[(j * 4)] = s_buffer[((j * 4) + 2)] = sum_base[(j * 2)];
          s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 3)] = sum_base[((j * 2) + 1)];
        }
        if(incv == 1){

          for(i = 0; i + 4 <= n; i += 4, v_base += 8){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[2] * scale_mask;
            v_3 = v_base[3] * scale_mask;
            v_4 = v_base[4] * scale_mask;
            v_5 = v_base[5] * scale_mask;
            v_6 = v_base[6] * scale_mask;
            v_7 = v_base[7] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            v_4 = v_4 * v_4;
            v_5 = v_5 * v_5;
            v_6 = v_6 * v_6;
            v_7 = v_7 * v_7;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              tmp_BLP.d = v_2;
              tmp_BLP.l |= 1;
              q_2 = s_2 + tmp_BLP.d;
              tmp_BLP.d = v_3;
              tmp_BLP.l |= 1;
              q_3 = s_3 + tmp_BLP.d;
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
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              tmp_BLP.d = v_4;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_5;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              tmp_BLP.d = v_6;
              tmp_BLP.l |= 1;
              q_2 = s_2 + tmp_BLP.d;
              tmp_BLP.d = v_7;
              tmp_BLP.l |= 1;
              q_3 = s_3 + tmp_BLP.d;
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              q_2 = s_2 - q_2;
              q_3 = s_3 - q_3;
              v_4 = v_4 + q_0;
              v_5 = v_5 + q_1;
              v_6 = v_6 + q_2;
              v_7 = v_7 + q_3;
            }
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 2)] + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 3)] = s_buffer[((j * 4) + 3)] + tmp_BLP.d;
          }
          if(i + 2 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[2] * scale_mask;
            v_3 = v_base[3] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              tmp_BLP.d = v_2;
              tmp_BLP.l |= 1;
              q_2 = s_2 + tmp_BLP.d;
              tmp_BLP.d = v_3;
              tmp_BLP.l |= 1;
              q_3 = s_3 + tmp_BLP.d;
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
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 2)] + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 3)] = s_buffer[((j * 4) + 3)] + tmp_BLP.d;
            i += 2, v_base += 4;
          }
          if(i + 1 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
            }
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.d;
            i += 1, v_base += 2;
          }
        }else{

          for(i = 0; i + 4 <= n; i += 4, v_base += (incv * 8)){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[(incv * 2)] * scale_mask;
            v_3 = v_base[((incv * 2) + 1)] * scale_mask;
            v_4 = v_base[(incv * 4)] * scale_mask;
            v_5 = v_base[((incv * 4) + 1)] * scale_mask;
            v_6 = v_base[(incv * 6)] * scale_mask;
            v_7 = v_base[((incv * 6) + 1)] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            v_4 = v_4 * v_4;
            v_5 = v_5 * v_5;
            v_6 = v_6 * v_6;
            v_7 = v_7 * v_7;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              tmp_BLP.d = v_2;
              tmp_BLP.l |= 1;
              q_2 = s_2 + tmp_BLP.d;
              tmp_BLP.d = v_3;
              tmp_BLP.l |= 1;
              q_3 = s_3 + tmp_BLP.d;
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
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              tmp_BLP.d = v_4;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_5;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              tmp_BLP.d = v_6;
              tmp_BLP.l |= 1;
              q_2 = s_2 + tmp_BLP.d;
              tmp_BLP.d = v_7;
              tmp_BLP.l |= 1;
              q_3 = s_3 + tmp_BLP.d;
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              s_buffer[((j * 4) + 2)] = q_2;
              s_buffer[((j * 4) + 3)] = q_3;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              q_2 = s_2 - q_2;
              q_3 = s_3 - q_3;
              v_4 = v_4 + q_0;
              v_5 = v_5 + q_1;
              v_6 = v_6 + q_2;
              v_7 = v_7 + q_3;
            }
            tmp_BLP.d = v_4;
            tmp_BLP.l |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.d;
            tmp_BLP.d = v_5;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.d;
            tmp_BLP.d = v_6;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 2)] + tmp_BLP.d;
            tmp_BLP.d = v_7;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 3)] = s_buffer[((j * 4) + 3)] + tmp_BLP.d;
          }
          if(i + 2 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_2 = v_base[(incv * 2)] * scale_mask;
            v_3 = v_base[((incv * 2) + 1)] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            v_2 = v_2 * v_2;
            v_3 = v_3 * v_3;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              s_2 = s_buffer[((j * 4) + 2)];
              s_3 = s_buffer[((j * 4) + 3)];
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              tmp_BLP.d = v_2;
              tmp_BLP.l |= 1;
              q_2 = s_2 + tmp_BLP.d;
              tmp_BLP.d = v_3;
              tmp_BLP.l |= 1;
              q_3 = s_3 + tmp_BLP.d;
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
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.d;
            tmp_BLP.d = v_2;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 2)] = s_buffer[((j * 4) + 2)] + tmp_BLP.d;
            tmp_BLP.d = v_3;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 3)] = s_buffer[((j * 4) + 3)] + tmp_BLP.d;
            i += 2, v_base += (incv * 4);
          }
          if(i + 1 <= n){
            v_0 = v_base[0] * scale_mask;
            v_1 = v_base[1] * scale_mask;
            v_0 = v_0 * v_0;
            v_1 = v_1 * v_1;
            for(j = 0; j < fold - 1; j++){
              s_0 = s_buffer[(j * 4)];
              s_1 = s_buffer[((j * 4) + 1)];
              tmp_BLP.d = v_0;
              tmp_BLP.l |= 1;
              q_0 = s_0 + tmp_BLP.d;
              tmp_BLP.d = v_1;
              tmp_BLP.l |= 1;
              q_1 = s_1 + tmp_BLP.d;
              s_buffer[(j * 4)] = q_0;
              s_buffer[((j * 4) + 1)] = q_1;
              q_0 = s_0 - q_0;
              q_1 = s_1 - q_1;
              v_0 = v_0 + q_0;
              v_1 = v_1 + q_1;
            }
            tmp_BLP.d = v_0;
            tmp_BLP.l |= 1;
            s_buffer[(j * 4)] = s_buffer[(j * 4)] + tmp_BLP.d;
            tmp_BLP.d = v_1;
            tmp_BLP.l |= 1;
            s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + tmp_BLP.d;
            i += 1, v_base += (incv * 2);
          }
        }
        for(j = 0; j < fold; j += 1){
          q_0 = ((double*)sum)[(j * 2)];
          s_buffer[(j * 4)] = s_buffer[(j * 4)] + (s_buffer[((j * 4) + 2)] - q_0);
          q_0 = ((double*)sum)[((j * 2) + 1)];
          s_buffer[((j * 4) + 1)] = s_buffer[((j * 4) + 1)] + (s_buffer[((j * 4) + 3)] - q_0);
          ((double*)sum)[(j * 2)] = s_buffer[(j * 4)];
          ((double*)sum)[((j * 2) + 1)] = s_buffer[((j * 4) + 1)];
        }
        RESET_DAZ_FLAG
        return;
      }
    }
    //[[[end]]]
  }
#endif
