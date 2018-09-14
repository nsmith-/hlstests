#include "ap_lorentz.h"

typedef ap_lorentz<ap_fixed<14, 12> > lorentz;

const size_t n_cands = 8;
// n choose 3
const size_t n_perm = n_cands*(n_cands-1)*(n_cands-2)/6;

typedef struct {
  lorentz::xyzt_t pt;
  lorentz::eta_t eta;
  lorentz::phi_t phi;
} cand_in;

void cand_perm_mass(cand_in input[n_cands], lorentz::xyzt_t mass_out[n_perm]);
