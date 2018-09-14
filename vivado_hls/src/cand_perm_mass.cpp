#include "cand_perm_mass.h"

void cand_perm_mass(cand_in input[n_cands], lorentz::xyzt_t mass_out[n_perm])
{
#pragma HLS PIPELINE II=6
#pragma HLS INTERFACE ap_none port=mass_out
#pragma HLS ARRAY_PARTITION variable=input complete dim=0
#pragma HLS ARRAY_PARTITION variable=mass_out complete dim=0

  size_t im = 0;
  for(size_t i=0; i<n_cands-2; ++i) {
    lorentz cand_i = lorentz(input[i].pt, input[i].eta, input[i].phi, 0.);
    for(size_t j=i+1; j<n_cands-1; ++j) {
      lorentz cand_j = lorentz(input[j].pt, input[j].eta, input[j].phi, 0.);
      lorentz cand_ij = cand_i + cand_j;
      for(size_t k=j+1; k<n_cands; ++k) {
        lorentz cand_k = lorentz(input[k].pt, input[k].eta, input[k].phi, 0.);
        lorentz cand_ijk = cand_ij + cand_k;
        mass_out[im++] = cand_ijk.mass();
      }
    }
  }
}

