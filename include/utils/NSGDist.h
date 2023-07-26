// taken from ParlayANN 6/13/2023
//
// Created by 付聪 on 2017/6/21.
//
//Edited

#ifndef EFANNA2E_DISTANCE_H
#define EFANNA2E_DISTANCE_H

#include <x86intrin.h>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <math.h>
#include "parlay/parallel.h"
#include "parlay/primitives.h"

#include "threadlocal.h"

extern bool report_stats;

namespace efanna2e{

  // atomic_sum_counter<size_t> distance_calls;
  
  enum Metric{
    L2 = 0,
    INNER_PRODUCT = 1,
    FAST_L2 = 2,
    PQ = 3
  };
    class Distance {
    public:
        virtual float compare(const float* a, const float* b, unsigned length) const = 0;
        virtual ~Distance() {}
    };

    class DistanceL2 : public Distance{
    public:
        float compare(const float* a, const float* b, unsigned size) const {
          //std::cout << "inside compare" << std::endl;
            float result = 0;

#ifdef __GNUC__
//std::cout << "gnuc defined\n";
#ifdef __AVX__
//std::cout << "avx defined\n";
  #define AVX_L2SQR(addr1, addr2, dest, tmp1, tmp2) \
      tmp1 = _mm256_loadu_ps(addr1);\
      tmp2 = _mm256_loadu_ps(addr2);\
      tmp1 = _mm256_sub_ps(tmp1, tmp2); \
      tmp1 = _mm256_mul_ps(tmp1, tmp1); \
      dest = _mm256_add_ps(dest, tmp1);

      __m256 sum;
      __m256 l0, l1;
      __m256 r0, r1;
      unsigned D = (size + 7) & ~7U;
      unsigned DR = D % 16;
      unsigned DD = D - DR;
      const float *l = a;
      const float *r = b;
      const float *e_l = l + DD;
      const float *e_r = r + DD;
      float unpack[8] __attribute__ ((aligned (32))) = {0, 0, 0, 0, 0, 0, 0, 0};

      sum = _mm256_loadu_ps(unpack);
      if(DR){AVX_L2SQR(e_l, e_r, sum, l0, r0);}

      for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
      	AVX_L2SQR(l, r, sum, l0, r0);
      	AVX_L2SQR(l + 8, r + 8, sum, l1, r1);
      }
      _mm256_storeu_ps(unpack, sum);

      // std::cout << "here" << std::endl;
      // std::cout << unpack[0] << " " << unpack[1] << " " << unpack[2] << " " << unpack[3] << " " << unpack[4] << std::endl;
      result = unpack[0] + unpack[1] + unpack[2] + unpack[3] + unpack[4] + unpack[5] + unpack[6] + unpack[7];

#else
//std::cout << "avx not defined\n";
#ifdef __SSE2__
//std::cout << "sse defined\n";
  #define SSE_L2SQR(addr1, addr2, dest, tmp1, tmp2) \
          tmp1 = _mm_load_ps(addr1);\
          tmp2 = _mm_load_ps(addr2);\
          tmp1 = _mm_sub_ps(tmp1, tmp2); \
          tmp1 = _mm_mul_ps(tmp1, tmp1); \
          dest = _mm_add_ps(dest, tmp1);

  __m128 sum;
  __m128 l0, l1, l2, l3;
  __m128 r0, r1, r2, r3;
  unsigned D = (size + 3) & ~3U;
  unsigned DR = D % 16;
  unsigned DD = D - DR;
  const float *l = a;
  const float *r = b;
  const float *e_l = l + DD;
  const float *e_r = r + DD;
  float unpack[4] __attribute__ ((aligned (16))) = {0, 0, 0, 0};

  sum = _mm_load_ps(unpack);
  switch (DR) {
      case 12:
      SSE_L2SQR(e_l+8, e_r+8, sum, l2, r2);
      case 8:
      SSE_L2SQR(e_l+4, e_r+4, sum, l1, r1);
      case 4:
      SSE_L2SQR(e_l, e_r, sum, l0, r0);
    default:
      break;
  }
  for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
      SSE_L2SQR(l, r, sum, l0, r0);
      SSE_L2SQR(l + 4, r + 4, sum, l1, r1);
      SSE_L2SQR(l + 8, r + 8, sum, l2, r2);
      SSE_L2SQR(l + 12, r + 12, sum, l3, r3);
  }
  _mm_storeu_ps(unpack, sum);

  //std::cout << "unpack 4 " << unpack[0] << " " << unpack[1] << " " << unpack[2] << " " << unpack[3] << std::endl;
  result += unpack[0] + unpack[1] + unpack[2] + unpack[3];

//normal distance
#else
     // std::cout << "normal roll loop" << std::endl;
      float diff0, diff1, diff2, diff3;
      const float* last = a + size;
      const float* unroll_group = last - 3;

      /* Process 4 items with each loop for efficiency. */
      while (a < unroll_group) {
          diff0 = a[0] - b[0];
          diff1 = a[1] - b[1];
          diff2 = a[2] - b[2];
          diff3 = a[3] - b[3];
          result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
          a += 4;
          b += 4;
      }
      /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
      while (a < last) {
          diff0 = *a++ - *b++;
          result += diff0 * diff0;
      }
#endif
#endif
#endif
//std::cout << "about to return\n";
            return result;
        }
    };

  class DistanceInnerProduct : public Distance{
  public:
    float compare(const float* a, const float* b, unsigned size) const {
      float result = 0;
#ifdef __GNUC__
#ifdef __AVX__
      #define AVX_DOT(addr1, addr2, dest, tmp1, tmp2) \
          tmp1 = _mm256_loadu_ps(addr1);\
          tmp2 = _mm256_loadu_ps(addr2);\
          tmp1 = _mm256_mul_ps(tmp1, tmp2); \
          dest = _mm256_add_ps(dest, tmp1);

	  __m256 sum;
   	  __m256 l0, l1;
   	  __m256 r0, r1;
      unsigned D = (size + 7) & ~7U;
      unsigned DR = D % 16;
      unsigned DD = D - DR;
   	  const float *l = a;
   	  const float *r = b;
      const float *e_l = l + DD;
   	  const float *e_r = r + DD;
      float unpack[8] __attribute__ ((aligned (32))) = {0, 0, 0, 0, 0, 0, 0, 0};

      sum = _mm256_loadu_ps(unpack);
      if(DR){AVX_DOT(e_l, e_r, sum, l0, r0);}

      for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
	    AVX_DOT(l, r, sum, l0, r0);
	    AVX_DOT(l + 8, r + 8, sum, l1, r1);
      }
      _mm256_storeu_ps(unpack, sum);
      result = unpack[0] + unpack[1] + unpack[2] + unpack[3] + unpack[4] + unpack[5] + unpack[6] + unpack[7];

#else
#ifdef __SSE2__
      #define SSE_DOT(addr1, addr2, dest, tmp1, tmp2) \
          tmp1 = _mm128_loadu_ps(addr1);\
          tmp2 = _mm128_loadu_ps(addr2);\
          tmp1 = _mm128_mul_ps(tmp1, tmp2); \
          dest = _mm128_add_ps(dest, tmp1);
      __m128 sum;
      __m128 l0, l1, l2, l3;
      __m128 r0, r1, r2, r3;
      unsigned D = (size + 3) & ~3U;
      unsigned DR = D % 16;
      unsigned DD = D - DR;
      const float *l = a;
      const float *r = b;
      const float *e_l = l + DD;
      const float *e_r = r + DD;
      float unpack[4] __attribute__ ((aligned (16))) = {0, 0, 0, 0};

      sum = _mm_load_ps(unpack);
      switch (DR) {
          case 12:
          SSE_DOT(e_l+8, e_r+8, sum, l2, r2);
          case 8:
          SSE_DOT(e_l+4, e_r+4, sum, l1, r1);
          case 4:
          SSE_DOT(e_l, e_r, sum, l0, r0);
        default:
          break;
      }
      for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
          SSE_DOT(l, r, sum, l0, r0);
          SSE_DOT(l + 4, r + 4, sum, l1, r1);
          SSE_DOT(l + 8, r + 8, sum, l2, r2);
          SSE_DOT(l + 12, r + 12, sum, l3, r3);
      }
      _mm_storeu_ps(unpack, sum);
      result += unpack[0] + unpack[1] + unpack[2] + unpack[3];
#else

      float dot0, dot1, dot2, dot3;
      const float* last = a + size;
      const float* unroll_group = last - 3;

      /* Process 4 items with each loop for efficiency. */
      while (a < unroll_group) {
          dot0 = a[0] * b[0];
          dot1 = a[1] * b[1];
          dot2 = a[2] * b[2];
          dot3 = a[3] * b[3];
          result += dot0 + dot1 + dot2 + dot3;
          a += 4;
          b += 4;
      }
      /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
      while (a < last) {
          result += *a++ * *b++;
      }
#endif
#endif
#endif
      return result;
    }

  };
  class DistanceFastL2 : public DistanceInnerProduct{
   public:
    float norm(const float* a, unsigned size) const{
      float result = 0;
#ifdef __GNUC__
#ifdef __AVX__
#define AVX_L2NORM(addr, dest, tmp) \
    tmp = _mm256_loadu_ps(addr); \
    tmp = _mm256_mul_ps(tmp, tmp); \
    dest = _mm256_add_ps(dest, tmp);

    __m256 sum;
   	__m256 l0, l1;
    unsigned D = (size + 7) & ~7U;
    unsigned DR = D % 16;
    unsigned DD = D - DR;
    const float *l = a;
    const float *e_l = l + DD;
    float unpack[8] __attribute__ ((aligned (32))) = {0, 0, 0, 0, 0, 0, 0, 0};

    sum = _mm256_loadu_ps(unpack);
    if(DR){AVX_L2NORM(e_l, sum, l0);}
	for (unsigned i = 0; i < DD; i += 16, l += 16) {
      AVX_L2NORM(l, sum, l0);
      AVX_L2NORM(l + 8, sum, l1);
    }
    _mm256_storeu_ps(unpack, sum);
    result = unpack[0] + unpack[1] + unpack[2] + unpack[3] + unpack[4] + unpack[5] + unpack[6] + unpack[7];
#else
#ifdef __SSE2__
#define SSE_L2NORM(addr, dest, tmp) \
    tmp = _mm128_loadu_ps(addr); \
    tmp = _mm128_mul_ps(tmp, tmp); \
    dest = _mm128_add_ps(dest, tmp);

    __m128 sum;
    __m128 l0, l1, l2, l3;
    unsigned D = (size + 3) & ~3U;
    unsigned DR = D % 16;
    unsigned DD = D - DR;
    const float *l = a;
    const float *e_l = l + DD;
    float unpack[4] __attribute__ ((aligned (16))) = {0, 0, 0, 0};

    sum = _mm_load_ps(unpack);
    switch (DR) {
        case 12:
        SSE_L2NORM(e_l+8, sum, l2);
        case 8:
        SSE_L2NORM(e_l+4, sum, l1);
        case 4:
        SSE_L2NORM(e_l, sum, l0);
      default:
        break;
    }
    for (unsigned i = 0; i < DD; i += 16, l += 16) {
        SSE_L2NORM(l, sum, l0);
        SSE_L2NORM(l + 4, sum, l1);
        SSE_L2NORM(l + 8, sum, l2);
        SSE_L2NORM(l + 12, sum, l3);
    }
    _mm_storeu_ps(unpack, sum);
    result += unpack[0] + unpack[1] + unpack[2] + unpack[3];
#else
    float dot0, dot1, dot2, dot3;
    const float* last = a + size;
    const float* unroll_group = last - 3;

    /* Process 4 items with each loop for efficiency. */
    while (a < unroll_group) {
        dot0 = a[0] * a[0];
        dot1 = a[1] * a[1];
        dot2 = a[2] * a[2];
        dot3 = a[3] * a[3];
        result += dot0 + dot1 + dot2 + dot3;
        a += 4;
    }
    /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
    while (a < last) {
        result += (*a) * (*a);
        a++;
    }
#endif
#endif
#endif
      return result;
    }
    using DistanceInnerProduct::compare;
    float compare(const float* a, const float* b, float norm, unsigned size) const {//not implement
      float result = -2 * DistanceInnerProduct::compare(a, b, size);
      result += norm;
      return result;
    }
  };
}

class Distance{
  public:
    virtual std::string id(){return "generic";}
    virtual float distance(uint8_t *p, uint8_t *q, unsigned d){return 0;}
    virtual float distance(int8_t *p, int8_t *q, unsigned d){return 0;}
    virtual float distance(float *p, float *q, unsigned d){return 0;}
    //virtual ~Distance() = 0; //destructor?
    virtual float distance(uint8_t *p, float *q, unsigned d){return 0;}
    virtual float distance(int8_t *p, float *q, unsigned d){return 0;}
    virtual float distance(float *p, uint8_t *q, unsigned d){return 0;}
    virtual float distance(float *p, int8_t *q, unsigned d){return 0;}

};

struct Mips_Distance : public Distance{

  std::string id(){return "mips";}

  float distance(uint8_t *p, uint8_t *q, unsigned d){
    int result = 0;
    for(unsigned i=0; i<d; i++){ //changing type of i to unsigned to avoid compiler warnings
      result += ((int32_t) q[i]) *
                    ((int32_t) p[i]);
    }
    return -((float) result);
  }

  float distance(int8_t *p, int8_t *q, unsigned d){
    int result = 0;
    for(unsigned i=0; i<d; i++){ //changed type of i to unsigned to avoid compiler warnings
      result += ((int32_t) q[i]) *
                    ((int32_t) p[i]);
    }
    return -((float) result);
  }

  float distance(float *p, float *q, unsigned d){
      float result = 0;
      for(unsigned i=0; i<d; i++){
        result += (q[i]) * (p[i]);
      }
      return -result;
  }

};

struct EuclideanDistance : public Distance{

 threadlocal::buffer<float>* buf = nullptr;

  std::string id(){return "euclidean";}

  float distance(uint8_t *p, uint8_t *q, unsigned d){
   // std::cout << "calling uint8_t dist\n";
    int result = 0;
    for(unsigned i=0; i<d; i++){
      result += ((int32_t)((int16_t) q[i] - (int16_t) p[i])) *
                    ((int32_t)((int16_t) q[i] - (int16_t) p[i]));
    }
    return (float) result;
  }

  float distance(int8_t *p, int8_t *q, unsigned d){
  //  std::cout << "calling int8_t dist" << std::endl;
    int result = 0;
    for(unsigned i=0; i<d; i++){
      result += ((int32_t)((int16_t) q[i] - (int16_t) p[i])) *
                    ((int32_t)((int16_t) q[i] - (int16_t) p[i]));
    }
    return (float) result;
  }

  float distance(float *p, float *q, unsigned d){
      efanna2e::DistanceL2 distfunc;
   //   std::cout << "comparing distfunc" << std::endl;
      return distfunc.compare(p, q, d);
  }

  float distance(uint8_t *p, float *q, unsigned d){
      if (buf == nullptr || buf->length < d) {
        buf = new threadlocal::buffer<float>(d);
      }
      buf->write(p);
      return distance(buf->begin(), q, d);
  }

  float distance(int8_t *p, float *q, unsigned d){
      if (buf == nullptr || buf->length < d) {
        buf = new threadlocal::buffer<float>(d);
      }
      buf->write(p);
      return distance(buf->begin(), q, d);
  }

  float distance(float *p, uint8_t *q, unsigned d){
      return distance(q, p, d);
  }

  float distance(float *p, int8_t *q, unsigned d){
      return distance(q, p, d);
  }

};

//Euclidian distance to use if d < 36 (I believe that 30 is the bound, but 
//just to be safe)
struct EuclideanDistanceSmall : public Distance {
    std::string id() {return "euclidean_small";}

    float distance(uint8_t *p, uint8_t *q, unsigned d){
          // std::cout << "uint8: d: " << static_cast<int>(d) << std::endl;

    int result = 0;
    for(unsigned i=0; i<d; i++){
      result += ((int32_t)((int16_t) q[i] - (int16_t) p[i])) *
                    ((int32_t)((int16_t) q[i] - (int16_t) p[i]));
    }
    return (float) result;
  }

  float distance(int8_t *p, int8_t *q, unsigned d){
    int result = 0;
    for(unsigned i=0; i<d; i++){
      result += ((int32_t)((int16_t) q[i] - (int16_t) p[i])) *
                    ((int32_t)((int16_t) q[i] - (int16_t) p[i]));
    }
    return (float) result;
  }
  float distance(float* p, float* q, unsigned d) {
    float result = 0;
    for (unsigned i = 0; i < d; i++) {
      result += (p[i]-q[i])*(p[i]-q[i]);
    }
    return result;
  }

  float distance(uint8_t *p, float *q, unsigned d){
    std::cout << "Mixed distance call aborting1" << std::endl;
    abort(); return 0;}
  float distance(int8_t *p, float *q, unsigned d){
    std::cout << "Mixed distance call aborting2" << std::endl;

    abort();
  return 0;}
  float distance(float *p, uint8_t *q, unsigned d){
    std::cout << "Mixed distance call aborting3" << std::endl;
    abort();
    return 0;}
  float distance(float *p, int8_t *q, unsigned d){
    std::cout << "Mixed distance call aborting4" << std::endl;
    abort();

    return 0;}



};


#endif //EFANNA2E_DISTANCE_H
