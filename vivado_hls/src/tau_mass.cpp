#include "tau_mass.h"

void tau_mass(cand_in input[3], lorentz& output, lorentz::mag2_t& m2)
{
#pragma HLS PIPELINE II=6
#pragma HLS INTERFACE ap_none port=output
#pragma HLS ARRAY_PARTITION variable=input complete dim=0

  lorentz cand0 = lorentz(input[0].pt, input[0].eta, input[0].phi, 0.);
  lorentz cand1 = lorentz(input[1].pt, input[1].eta, input[1].phi, 0.);
  lorentz cand2 = lorentz(input[2].pt, input[2].eta, input[2].phi, 0.);
  output = cand0 + cand1 + cand2;

  lorentz::mag2_t m2_01 = lorentz::mass2_boosted2cand(input[0].pt, input[0].eta, input[0].phi, input[1].pt, input[1].eta, input[1].phi);
  lorentz::mag2_t m2_02 = lorentz::mass2_boosted2cand(input[0].pt, input[0].eta, input[0].phi, input[2].pt, input[2].eta, input[2].phi);
  lorentz::mag2_t m2_12 = lorentz::mass2_boosted2cand(input[1].pt, input[1].eta, input[1].phi, input[2].pt, input[2].eta, input[2].phi);
  m2 = m2_01 + m2_02 + m2_12;
}

