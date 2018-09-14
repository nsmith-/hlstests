#include "ap_lorentz.h"

typedef ap_lorentz<ap_fixed<14, 12> > lorentz;

typedef struct {
  lorentz::xyzt_t pt;
  lorentz::eta_t eta;
  lorentz::phi_t phi;
} cand_in;

void hls_demo(cand_in input[2], lorentz& output, lorentz::xyzt_t& m);
