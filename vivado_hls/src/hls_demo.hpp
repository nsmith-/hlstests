#include <cmath>
#include "ap_int.h"
#include "ap_fixed.h"
#include "hls_half.h"
#include "hls_math.h"

// Bit depth of input coordinates
#define BITD 10

typedef struct {
	ap_int<BITD> x;
	ap_int<BITD> y;
} in_cartesian_t;

typedef struct {
	ap_int<BITD+1> r; // ceil(BITD+0.5), LSB = LSB of input
	ap_fixed<BITD, 3> phi; // radians
} out_polar_t;

// implementations of the coordinate transform
typedef struct {
  out_polar_t cordic;
  out_polar_t hlsmath;
  out_polar_t fapprox;
} output_struct;
void to_polar_impl_cordic(in_cartesian_t input, out_polar_t& output);
void to_polar_impl_hlsmath(in_cartesian_t input, out_polar_t& output);
void to_polar_impl_fapprox(in_cartesian_t input, out_polar_t& output);

void hls_demo(in_cartesian_t input, output_struct& output);

// Preprocessor crap to make function usable to BITD = 17
// But a customized routine could save space since not all are needed
#define CORDIC_ITER BITD+1-2
typedef ap_fixed<BITD+1+4, BITD+1> cordic_r_t;
typedef ap_fixed<BITD+1+4, 3> cordic_phi_t;
const cordic_phi_t cordic_angles[16] = {
  // python
  // >> import math
  // >> [math.atan(pow(2, -x-1)) for x in range(16)]
  0.46364760900080615, 0.24497866312686414, 0.12435499454676144, 0.06241880999595735, 0.031239833430268277, 0.015623728620476831, 0.007812341060101111, 0.0039062301319669718, 0.0019531225164788188, 0.0009765621895593195, 0.0004882812111948983, 0.00024414062014936177, 0.00012207031189367021, 6.103515617420877e-05, 3.0517578115526096e-05, 1.5258789061315762e-05
};
const ap_ufixed<BITD+1, 0> cordic_scales[16] = {
  // python
  // >> from operator import mul
  // >> a = [pow(1+pow(2, -2*(x+1)), -.5) for x in range(16)]
  // >> [reduce(mul, a[:n+1]) for n in range(16)]
  0.8944271909999159, 0.8677218312746247, 0.8610211763152799, 0.8593444051498349, 0.8589251104651273, 0.8588202804030459, 0.8587940724877404, 0.8587875204839219, 0.8587858824814052, 0.8587854729806784, 0.8587853706054906, 0.8587853450116933, 0.858785338613244, 0.8587853370136317, 0.8587853366137286, 0.8587853365137528
};
