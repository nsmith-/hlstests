#include "hls_demo.h"

void hls_demo(cand_in input[2], lorentz& output, lorentz::mag2_t& m2)
{
#pragma HLS PIPELINE II=6
#pragma HLS INTERFACE ap_none port=output
#pragma HLS ARRAY_PARTITION variable=input complete dim=0

  lorentz cand0 = lorentz(input[0].pt, input[0].eta, input[0].phi, 0.);
  lorentz cand1 = lorentz(input[1].pt, input[1].eta, input[1].phi, 0.);
  output = cand0 + cand1;
  // m2 = output.mass2();
  m2 = lorentz::mass2_boosted2cand(input[0].pt, input[0].eta, input[0].phi, input[1].pt, input[1].eta, input[1].phi);
}

