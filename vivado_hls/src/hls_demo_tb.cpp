#include <cstdio>
#include <cstdlib>
#include "hls_demo.hpp"

int main(void)
{
  std::srand(42);

  in_cartesian_t val;
  output_struct output;
  
  for(size_t i=0; i<100; ++i) {
    uint32_t rnd = std::rand();
    val.x = rnd & ((1<<BITD)-1);
    val.y = (rnd>>BITD) & ((1<<BITD)-1);
    std::printf("In: x=% 4d y=% 4d\n", static_cast<int>(val.x), static_cast<int>(val.y));

    double x = val.x;
    double y = val.y;
    double r = std::sqrt(x*x+y*y);
    double phi = std::atan2(y, x);
    std::printf("  Out FP     : r=% 6.1f phi=% 1.4f\n", r, phi);

    hls_demo(val, output);
    double r_cordic = static_cast<double>(output.cordic.r);
    double phi_cordic = static_cast<double>(output.cordic.phi);
    std::printf("  Out cordic : r=% 4.0f delta/LSB=% 3.1f,  phi=% 1.4f delta/LSB=% 3.1f\n", r_cordic, (r_cordic-r), phi_cordic, (phi_cordic-phi)/pow(2, -BITD+3));
    double r_hlsmath = static_cast<double>(output.hlsmath.r);
    double phi_hlsmath = static_cast<double>(output.hlsmath.phi);
    std::printf("  Out hlsmath: r=% 4.0f delta/LSB=% 3.1f,  phi=% 1.4f delta/LSB=% 3.1f\n", r_hlsmath, (r_hlsmath-r), phi_hlsmath, (phi_hlsmath-phi)/pow(2, -BITD+3));
    double r_lorentz = static_cast<double>(output.lorentz.r);
    double phi_lorentz = static_cast<double>(output.lorentz.phi);
    std::printf("  Out lorentz: r=% 4.0f delta/LSB=% 3.1f,  phi=% 1.4f delta/LSB=% 3.1f\n", r_lorentz, (r_lorentz-r), phi_lorentz, (phi_lorentz-phi)/pow(2, -BITD+3));
  }

  std::printf("Test passed!\n");
  return 0;
}
