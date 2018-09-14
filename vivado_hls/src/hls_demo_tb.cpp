#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "hls_demo.h"

int main(void)
{
  std::srand(42);

  cand_in input[2];
  lorentz output;
  lorentz::xyzt_t m;
  std::cout << std::fixed;

  const size_t nprint = 10;
  std::ofstream csv("deltas.csv");
  csv << "d_px,d_py,d_pz,d_E,d_M" << std::endl;
  
  for(size_t i=0; i<1000; ++i) {
    input[0].pt  = lorentz::xyzt_t(std::rand()*190./RAND_MAX+10.);
    input[0].eta = lorentz::eta_t(std::rand()*6./RAND_MAX-3.);
    input[0].phi = lorentz::phi_t(M_PI*std::rand()*1.98/RAND_MAX-M_PI);

    input[1].pt  = lorentz::xyzt_t(std::rand()*190./RAND_MAX+10.);
    input[1].eta = lorentz::eta_t(std::rand()*4./RAND_MAX-2.);
    input[1].phi = lorentz::phi_t(M_PI*std::rand()*1.98/RAND_MAX-M_PI);

    std::cout << std::setprecision(4);
    if ( i < nprint ) std::cout << "input 0 pt=" << input[0].pt << ", eta=" << input[0].eta << ", phi=" << input[0].phi << std::endl;
    if ( i < nprint ) std::cout << "input 1 pt=" << input[1].pt << ", eta=" << input[1].eta << ", phi=" << input[1].phi << std::endl;

    std::cout << std::setprecision(8);
    hls_demo(input, output, m);

    if ( i < nprint ) std::cout << "output " << output << ", mass = " << m << std::endl;

    double px0 = (double) input[0].pt * cos((double) input[0].phi);
    double py0 = (double) input[0].pt * sin((double) input[0].phi);
    double pz0 = (double) input[0].pt * sinh((double) input[0].eta);
    double E0  = (double) input[0].pt * cosh((double) input[0].eta);

    double px1 = (double) input[1].pt * cos((double) input[1].phi);
    double py1 = (double) input[1].pt * sin((double) input[1].phi);
    double pz1 = (double) input[1].pt * sinh((double) input[1].eta);
    double E1  = (double) input[1].pt * cosh((double) input[1].eta);

    double fp_m = std::sqrt(std::max(pow(E0+E1,2)-pow(px0+px1,2)-pow(py0+py1,2)-pow(pz0+pz1,2),0.));
    
    double lsb = pow(2, lorentz::xyzt_t::iwidth-lorentz::xyzt_t::width);
    double d_px = ((double) output.x_ - (px0 + px1)) / lsb;
    double d_py = ((double) output.y_ - (py0 + py1)) / lsb;
    double d_pz = ((double) output.z_ - (pz0 + pz1)) / lsb;
    double d_E = ((double) output.t_ - (E0 + E1)) / lsb;
    double d_M = ((double) m - fp_m) / lsb;

    std::cout << std::setprecision(1);
    if ( i < nprint ) std::cout << "Delta (/LSB) x: " << d_px << ", y: " << d_py << ", z: " << d_pz << ", t: " << d_E << ", d_m: " << d_M << std::endl;
    if ( i < nprint ) std::cout << std::endl;
    csv << d_px << "," << d_py << "," << d_pz << "," << d_E << "," << d_M << std::endl;
  }

  std::printf("Test passed!\n");
  return 0;
}
