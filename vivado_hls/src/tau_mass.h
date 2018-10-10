#include "ap_lorentz.h"

typedef ap_lorentz<ap_fixed<15, 10> > lorentz;

typedef struct {
  lorentz::xyzt_t pt;
  lorentz::eta_t eta;
  lorentz::phi_t phi;
} cand_in;

void tau_mass(cand_in input[3], lorentz& output, lorentz::mag2_t& m2);
