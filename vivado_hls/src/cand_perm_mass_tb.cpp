#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "cand_perm_mass.h"

int main(void)
{
  std::srand(42);

  cand_in input[n_cands];
  lorentz::xyzt_t mass_out[n_perm];

  for(size_t i=0; i<n_cands; ++i) {
    input[i].pt  = lorentz::xyzt_t(std::rand()*190./RAND_MAX+10.);
    input[i].eta = lorentz::eta_t(std::rand()*6./RAND_MAX-3.);
    input[i].phi = lorentz::phi_t(M_PI*std::rand()*1.98/RAND_MAX-M_PI);
  }

  cand_perm_mass(input, mass_out);


  std::printf("Test passed!\n");
  return 0;
}
