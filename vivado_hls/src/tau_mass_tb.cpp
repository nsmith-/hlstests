#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "tau_mass.h"

int main(void)
{
  cand_in input[3];
  lorentz output;
  lorentz::mag2_t m2;

  std::ifstream csv("3prong.csv");
  if ( not csv.good() ) {
    std::cout << "Could not open input csv file" << std::endl;
    return 1;
  }
  csv.ignore(10000, '\n');

  std::ofstream csv_out("deltas.csv");
  csv_out << "mass2_reco,mass2_fp,mass2_lorentz,mass2_boosted" << std::endl;
  
  while(not csv.eof()) {
    double pt, eta, phi;
    char comma;

    csv >> pt >> comma >> eta >> comma >> phi >> comma;
    input[0].pt  = lorentz::xyzt_t(pt);
    input[0].eta = lorentz::eta_t(eta);
    input[0].phi = lorentz::phi_t(phi);

    csv >> pt >> comma >> eta >> comma >> phi >> comma;
    input[1].pt  = lorentz::xyzt_t(pt);
    input[1].eta = lorentz::eta_t(eta);
    input[1].phi = lorentz::phi_t(phi);

    csv >> pt >> comma >> eta >> comma >> phi >> comma;
    input[2].pt  = lorentz::xyzt_t(pt);
    input[2].eta = lorentz::eta_t(eta);
    input[2].phi = lorentz::phi_t(phi);

    tau_mass(input, output, m2);

    double reco_mass, mass;
    csv >> pt >> comma >> eta >> comma >> phi >> comma >> reco_mass >> comma >> mass;

    csv_out << reco_mass*reco_mass << "," << mass*mass << "," << output.mass2() << "," << m2 << std::endl;
    if ( std::abs(((double) m2)-mass*mass) > .01 ) {
      lorentz::debug_ = true;
      tau_mass(input, output, m2);
      lorentz::debug_ = false;
      std::cout << "Reco pT: " << pt << ", HLS: " << output.pt() << std::endl;
      std::cout << "Reco m2: " << mass*mass << ", HLS: " << m2 << ", HLS2: " << output.mass2() << std::endl;
      std::cout << std::endl;
    }
  }

  std::cout << "Test passed!" << std::endl;
  return 0;
}
