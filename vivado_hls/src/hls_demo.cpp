#include "hls_demo.hpp"
#include <iostream>

void hls_demo(in_cartesian_t input, output_struct& output)
{
#pragma HLS PIPELINE II=6
#pragma HLS INTERFACE ap_none port=output

  to_polar_impl_cordic(input, output.cordic);
  to_polar_impl_hlsmath(input, output.hlsmath);
  to_polar_impl_fapprox(input, output.fapprox);

	return;
}

void to_polar_impl_cordic(in_cartesian_t input, out_polar_t& output) {
#pragma HLS INLINE off
  assert( BITD<=17 );

  cordic_r_t x, y;
  cordic_phi_t phi;

  // Coarse rotate to [-pi/4,pi/4)
  if ( input.x-input.y >= 0 ) {
    if ( input.x+input.y >= 0 ) {
      // East (Correct) quadrant
      x = input.x;
      y = input.y;
      phi = 0.;
    }
    else {
      // South, rotate by +pi/2
      x = -input.y;
      y = input.x;
      phi = M_PI/2.;
    }
  }
  else {
    if ( input.x+input.y >= 0 ) {
      // North, rotate by -pi/2
      x = input.y;
      y = -input.x;
      phi = -M_PI/2.;
    }
    else {
      // West, rotate by pi
      x = -input.x;
      y = -input.y;
      phi = M_PI;
    }
  }

  // std::cout << "Starting cordic" << std::endl;
  // std::cout << "  Iter 0, x: " << x << ", y: " << y << ", phi: " << phi << std::endl;

  cordic_r_t xtmp, ytmp;
  for(int i=1; i<=CORDIC_ITER; ++i) {
    xtmp = (y>=0) ? x + (y>>i) : x - (y>>i);
    ytmp = (y>=0) ? y - (x>>i) : y + (x>>i);
    phi = (y>=0) ? phi-cordic_angles[i-1] : phi+cordic_angles[i-1];
    x = xtmp;
    y = ytmp;
    // std::cout << "  Iter " << i << ", x: " << x << ", y: " << y << ", phi: " << phi << std::endl;
  }

  output.r = x * cordic_scales[CORDIC_ITER];
  output.phi = (phi > cordic_phi_t(M_PI)) ? cordic_phi_t(2*M_PI)-phi: -phi;
  // std::cout << "  Output r: " << output.r << ", phi: " << output.phi << std::endl;
}

void to_polar_impl_hlsmath(in_cartesian_t input, out_polar_t& output) {
#pragma HLS INLINE off
  assert( BITD<(10+1) ); // mantissa of half-precision = 10 bits

  float x_fp = static_cast<float>(input.x);
  float y_fp = static_cast<float>(input.y);

  float r_fp = hls::sqrtf( x_fp*x_fp + y_fp*y_fp );
  float phi_fp = hls::atan2f(y_fp, x_fp);

  output.r = static_cast<ap_int<BITD+1> >(r_fp);
  output.phi = static_cast<ap_fixed<BITD, 2> >(phi_fp);
}

void to_polar_impl_fapprox(in_cartesian_t input, out_polar_t& output) {
#pragma HLS INLINE off
  // TODO
  int32_t x = input.x;
  int32_t y = input.y;
  output.r = hls::sqrt(x*x+y*y);
}
