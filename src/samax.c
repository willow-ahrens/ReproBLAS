#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "Common/Common.h"
#include <immintrin.h>
#include <emmintrin.h>


#if defined( __AVX__ )
  float samax(int n, float* v, int incv){
    __m256 mask_ABS; AVX_ABS_MASKS(mask_ABS);
    float tmp_max[8] __attribute__((aligned(32)));
    int i;
    float max;
    float* max_ptr = (float*) &max;

    __m256 v_0, v_1, v_2, v_3;
    __m256 m_0;
    m_0 = _mm256_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v + 8), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_loadu_ps(v + 16), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_loadu_ps(v + 24), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
      }
      if(i + 16 <= n){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_loadu_ps(v + 8), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        i += 16, v += 16;
      }
      if(i + 8 <= n){
        v_0 = _mm256_and_ps(_mm256_loadu_ps(v), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        i += 8, v += 8;
      }
      if(i < n){
        v_0 = _mm256_and_ps(_mm256_set_ps(0, (n - i)>6?v[6]:0, (n - i)>5?v[5]:0, (n - i)>4?v[4]:0, (n - i)>3?v[3]:0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
      }
    }else{

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32)){
        v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
        v_2 = _mm256_and_ps(_mm256_set_ps(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)], v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]), mask_ABS);
        v_3 = _mm256_and_ps(_mm256_set_ps(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)], v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        m_0 = _mm256_max_ps(m_0, v_2);
        m_0 = _mm256_max_ps(m_0, v_3);
      }
      if(i + 16 <= n){
        v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm256_and_ps(_mm256_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)], v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        m_0 = _mm256_max_ps(m_0, v_1);
        i += 16, v += (incv * 16);
      }
      if(i + 8 <= n){
        v_0 = _mm256_and_ps(_mm256_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)], v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        m_0 = _mm256_max_ps(m_0, v_0);
        i += 8, v += (incv * 8);
      }
      if(i < n){
        v_0 = _mm256_and_ps(_mm256_set_ps(0, (n - i)>6?v[(incv * 6)]:0, (n - i)>5?v[(incv * 5)]:0, (n - i)>4?v[(incv * 4)]:0, (n - i)>3?v[(incv * 3)]:0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
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
  float samax(int n, float* v, int incv){
    __m128 mask_ABS; SSE_ABS_MASKS(mask_ABS);
    float tmp_max[4] __attribute__((aligned(16)));
    int i;
    float max;
    float* max_ptr = (float*) &max;

    __m128 v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7;
    __m128 m_0;
    m_0 = _mm_setzero_ps();

    if(incv == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32){
        v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v + 4), mask_ABS);
        v_2 = _mm_and_ps(_mm_loadu_ps(v + 8), mask_ABS);
        v_3 = _mm_and_ps(_mm_loadu_ps(v + 12), mask_ABS);
        v_4 = _mm_and_ps(_mm_loadu_ps(v + 16), mask_ABS);
        v_5 = _mm_and_ps(_mm_loadu_ps(v + 20), mask_ABS);
        v_6 = _mm_and_ps(_mm_loadu_ps(v + 24), mask_ABS);
        v_7 = _mm_and_ps(_mm_loadu_ps(v + 28), mask_ABS);
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
        v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v + 4), mask_ABS);
        v_2 = _mm_and_ps(_mm_loadu_ps(v + 8), mask_ABS);
        v_3 = _mm_and_ps(_mm_loadu_ps(v + 12), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        i += 16, v += 16;
      }
      if(i + 8 <= n){
        v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
        v_1 = _mm_and_ps(_mm_loadu_ps(v + 4), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        i += 8, v += 8;
      }
      if(i + 4 <= n){
        v_0 = _mm_and_ps(_mm_loadu_ps(v), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        i += 4, v += 4;
      }
      if(i < n){
        v_0 = _mm_and_ps(_mm_set_ps(0, (n - i)>2?v[2]:0, (n - i)>1?v[1]:0, v[0]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
      }
    }else{

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32)){
        v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
        v_2 = _mm_and_ps(_mm_set_ps(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
        v_3 = _mm_and_ps(_mm_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]), mask_ABS);
        v_4 = _mm_and_ps(_mm_set_ps(v[(incv * 19)], v[(incv * 18)], v[(incv * 17)], v[(incv * 16)]), mask_ABS);
        v_5 = _mm_and_ps(_mm_set_ps(v[(incv * 23)], v[(incv * 22)], v[(incv * 21)], v[(incv * 20)]), mask_ABS);
        v_6 = _mm_and_ps(_mm_set_ps(v[(incv * 27)], v[(incv * 26)], v[(incv * 25)], v[(incv * 24)]), mask_ABS);
        v_7 = _mm_and_ps(_mm_set_ps(v[(incv * 31)], v[(incv * 30)], v[(incv * 29)], v[(incv * 28)]), mask_ABS);
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
        v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
        v_2 = _mm_and_ps(_mm_set_ps(v[(incv * 11)], v[(incv * 10)], v[(incv * 9)], v[(incv * 8)]), mask_ABS);
        v_3 = _mm_and_ps(_mm_set_ps(v[(incv * 15)], v[(incv * 14)], v[(incv * 13)], v[(incv * 12)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        m_0 = _mm_max_ps(m_0, v_2);
        m_0 = _mm_max_ps(m_0, v_3);
        i += 16, v += (incv * 16);
      }
      if(i + 8 <= n){
        v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        v_1 = _mm_and_ps(_mm_set_ps(v[(incv * 7)], v[(incv * 6)], v[(incv * 5)], v[(incv * 4)]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        m_0 = _mm_max_ps(m_0, v_1);
        i += 8, v += (incv * 8);
      }
      if(i + 4 <= n){
        v_0 = _mm_and_ps(_mm_set_ps(v[(incv * 3)], v[(incv * 2)], v[incv], v[0]), mask_ABS);
        m_0 = _mm_max_ps(m_0, v_0);
        i += 4, v += (incv * 4);
      }
      if(i < n){
        v_0 = _mm_and_ps(_mm_set_ps(0, (n - i)>2?v[(incv * 2)]:0, (n - i)>1?v[incv]:0, v[0]), mask_ABS);
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
  float samax(int n, float* v, int incv){
    int i;
    float max;
    float* max_ptr = (float*) &max;

    float v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9, v_10, v_11, v_12, v_13, v_14, v_15, v_16, v_17, v_18, v_19, v_20, v_21, v_22, v_23, v_24, v_25, v_26, v_27, v_28, v_29, v_30, v_31;
    float m_0;
    m_0 = 0;

    if(incv == 1){

      for(i = 0; i + 32 <= n; i += 32, v += 32){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        v_2 = fabs(v[2]);
        v_3 = fabs(v[3]);
        v_4 = fabs(v[4]);
        v_5 = fabs(v[5]);
        v_6 = fabs(v[6]);
        v_7 = fabs(v[7]);
        v_8 = fabs(v[8]);
        v_9 = fabs(v[9]);
        v_10 = fabs(v[10]);
        v_11 = fabs(v[11]);
        v_12 = fabs(v[12]);
        v_13 = fabs(v[13]);
        v_14 = fabs(v[14]);
        v_15 = fabs(v[15]);
        v_16 = fabs(v[16]);
        v_17 = fabs(v[17]);
        v_18 = fabs(v[18]);
        v_19 = fabs(v[19]);
        v_20 = fabs(v[20]);
        v_21 = fabs(v[21]);
        v_22 = fabs(v[22]);
        v_23 = fabs(v[23]);
        v_24 = fabs(v[24]);
        v_25 = fabs(v[25]);
        v_26 = fabs(v[26]);
        v_27 = fabs(v[27]);
        v_28 = fabs(v[28]);
        v_29 = fabs(v[29]);
        v_30 = fabs(v[30]);
        v_31 = fabs(v[31]);
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
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        v_2 = fabs(v[2]);
        v_3 = fabs(v[3]);
        v_4 = fabs(v[4]);
        v_5 = fabs(v[5]);
        v_6 = fabs(v[6]);
        v_7 = fabs(v[7]);
        v_8 = fabs(v[8]);
        v_9 = fabs(v[9]);
        v_10 = fabs(v[10]);
        v_11 = fabs(v[11]);
        v_12 = fabs(v[12]);
        v_13 = fabs(v[13]);
        v_14 = fabs(v[14]);
        v_15 = fabs(v[15]);
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
        i += 16, v += 16;
      }
      if(i + 8 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        v_2 = fabs(v[2]);
        v_3 = fabs(v[3]);
        v_4 = fabs(v[4]);
        v_5 = fabs(v[5]);
        v_6 = fabs(v[6]);
        v_7 = fabs(v[7]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        i += 8, v += 8;
      }
      if(i + 4 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        v_2 = fabs(v[2]);
        v_3 = fabs(v[3]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        i += 4, v += 4;
      }
      if(i + 2 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[1]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        i += 2, v += 2;
      }
      if(i + 1 <= n){
        v_0 = fabs(v[0]);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += 1;
      }
    }else{

      for(i = 0; i + 32 <= n; i += 32, v += (incv * 32)){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        v_2 = fabs(v[(incv * 2)]);
        v_3 = fabs(v[(incv * 3)]);
        v_4 = fabs(v[(incv * 4)]);
        v_5 = fabs(v[(incv * 5)]);
        v_6 = fabs(v[(incv * 6)]);
        v_7 = fabs(v[(incv * 7)]);
        v_8 = fabs(v[(incv * 8)]);
        v_9 = fabs(v[(incv * 9)]);
        v_10 = fabs(v[(incv * 10)]);
        v_11 = fabs(v[(incv * 11)]);
        v_12 = fabs(v[(incv * 12)]);
        v_13 = fabs(v[(incv * 13)]);
        v_14 = fabs(v[(incv * 14)]);
        v_15 = fabs(v[(incv * 15)]);
        v_16 = fabs(v[(incv * 16)]);
        v_17 = fabs(v[(incv * 17)]);
        v_18 = fabs(v[(incv * 18)]);
        v_19 = fabs(v[(incv * 19)]);
        v_20 = fabs(v[(incv * 20)]);
        v_21 = fabs(v[(incv * 21)]);
        v_22 = fabs(v[(incv * 22)]);
        v_23 = fabs(v[(incv * 23)]);
        v_24 = fabs(v[(incv * 24)]);
        v_25 = fabs(v[(incv * 25)]);
        v_26 = fabs(v[(incv * 26)]);
        v_27 = fabs(v[(incv * 27)]);
        v_28 = fabs(v[(incv * 28)]);
        v_29 = fabs(v[(incv * 29)]);
        v_30 = fabs(v[(incv * 30)]);
        v_31 = fabs(v[(incv * 31)]);
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
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        v_2 = fabs(v[(incv * 2)]);
        v_3 = fabs(v[(incv * 3)]);
        v_4 = fabs(v[(incv * 4)]);
        v_5 = fabs(v[(incv * 5)]);
        v_6 = fabs(v[(incv * 6)]);
        v_7 = fabs(v[(incv * 7)]);
        v_8 = fabs(v[(incv * 8)]);
        v_9 = fabs(v[(incv * 9)]);
        v_10 = fabs(v[(incv * 10)]);
        v_11 = fabs(v[(incv * 11)]);
        v_12 = fabs(v[(incv * 12)]);
        v_13 = fabs(v[(incv * 13)]);
        v_14 = fabs(v[(incv * 14)]);
        v_15 = fabs(v[(incv * 15)]);
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
        i += 16, v += (incv * 16);
      }
      if(i + 8 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        v_2 = fabs(v[(incv * 2)]);
        v_3 = fabs(v[(incv * 3)]);
        v_4 = fabs(v[(incv * 4)]);
        v_5 = fabs(v[(incv * 5)]);
        v_6 = fabs(v[(incv * 6)]);
        v_7 = fabs(v[(incv * 7)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        m_0 = (m_0 > v_4? m_0: v_4);
        m_0 = (m_0 > v_5? m_0: v_5);
        m_0 = (m_0 > v_6? m_0: v_6);
        m_0 = (m_0 > v_7? m_0: v_7);
        i += 8, v += (incv * 8);
      }
      if(i + 4 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        v_2 = fabs(v[(incv * 2)]);
        v_3 = fabs(v[(incv * 3)]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        m_0 = (m_0 > v_2? m_0: v_2);
        m_0 = (m_0 > v_3? m_0: v_3);
        i += 4, v += (incv * 4);
      }
      if(i + 2 <= n){
        v_0 = fabs(v[0]);
        v_1 = fabs(v[incv]);
        m_0 = (m_0 > v_0? m_0: v_0);
        m_0 = (m_0 > v_1? m_0: v_1);
        i += 2, v += (incv * 2);
      }
      if(i + 1 <= n){
        v_0 = fabs(v[0]);
        m_0 = (m_0 > v_0? m_0: v_0);
        i += 1, v += incv;
      }
    }
    max_ptr[0] = m_0;
    return max;
  }
#endif