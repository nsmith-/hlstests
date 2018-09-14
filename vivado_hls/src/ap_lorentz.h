/*
 * ap_lorentz: A class for Lorentz vectors in HLS
 *
 * Internal storage is cartesian representation, with methods
 * to construct from both cartesian and polar/hyperbolic inputs
 *   i.e. (x,y,z,t) or (pt,eta,phi,E)
 *
 * The bitwidth of the storage is templated, so make sure the
 * chosen ap_fixed type is sufficient, e.g. ap_fixed<12, 10> will
 * provide values up to 1024 with an LSB of 0.25.
 *
 * This file is self-contained, except for the tables in
 * lorentz_tables.h, which will need to be in the include path.
 * 
 * - Nick Smith <nick.smith@cern.ch>
 */
#ifndef AP_LORENTZ_H
#define AP_LORENTZ_H

#include <ostream>
#include "ap_fixed.h"
#include "hls_dsp.h"
#include "hls/utils/x_hls_utils.h"

#include "lorentz_tables.h"


// T must be ap_fixed type
template<typename T>
class ap_lorentz {
public:
  typedef T xyzt_t;
  // 3-magnitude at most sqrt(3) larger than max xyzt_t
  typedef ap_fixed<xyzt_t::width+1, xyzt_t::width+1> mag_t;
  // eta range: [-8, 8)
  typedef ap_fixed<xyzt_t::width, 4> eta_t;
  // phi range: [-4, 4) (obviously |phi|<pi in reality)
  typedef ap_fixed<xyzt_t::width, 3> phi_t;
  // Squared quantities can be quite large
  typedef ap_ufixed<xyzt_t::width*2-1, xyzt_t::iwidth*2-1> mag2_t;
  // TODO: allow also scaled eta/phi (units of pi*1 rad)
  // This would make delta-phi a trivial subtraction

  xyzt_t x_, y_, z_, t_;

  // Default
  ap_lorentz() : x_(0.), y_(0.), z_(0.), t_(0.) {};

  // Construct from cartesian
  // (x,y,z,t) or (px, py, pz, E)
  ap_lorentz(xyzt_t x, xyzt_t y, xyzt_t z, xyzt_t t) {
    x_ = x;
    y_ = y;
    z_ = z;
    t_ = t;
  };

  // Construct from (pt, eta, phi, E)
  //  * Eta and phi must be same bit depth as pt, but use
  //    only 4 or 3 integer bits, respectively
  //    (i.e. -8 < eta <= 8, -4 < phi <= 4)
  //
  //  * If E=0, it is set to |p| (i.e. assume M=0)
  //
  //  * Warning: z_ and t_ can easily overlow if pt and |eta| are large
  //    Up to ceil(log2(cosh(8))) = 4 extra bits would be needed
  ap_lorentz(xyzt_t pt, eta_t eta, phi_t phi, xyzt_t E = 0.) {
    polar_to_cart(pt, phi, x_, y_);
    xyzt_t mag;
    hyperb_to_cart(pt, eta, mag, z_);
    t_ = (E<mag) ? mag : E;
  };

  // Return transverse component, i.e. p_T or rho
  // If you want phi also, use pt_phi instead
  mag_t pt() {
    mag_t pt_out;
    phi_t _drop;
    cart_to_polar(x_, y_, pt_out, _drop);
    return pt_out;
  };

  // Calculate p_T and phi angle at same time
  void pt_phi(mag_t& pt_out, phi_t& phi_out) {
    cart_to_polar(x_, y_, pt_out, phi_out);
  };

  // Return pseudorapidity
  // If you want pT as well, use pt_eta_phi
  eta_t eta() {
    // TODO
    return eta_t(0.);
  };

  // Calculate polar/hyperbolic coordinates
  void pt_eta_phi(mag_t& pt_out, eta_t& eta_out, phi_t& phi_out) {
    pt_phi(pt_out, phi_out);
    // TODO: eta = asinh(pz/pt)
  };

  // Rapidity
  eta_t rapidity() {
    // TODO: rapidity = atanh(pz/E)
  };

  // Invariant mass squared
  // This is cheaper than mass, so use it for a cut
  mag2_t mass2() {
    return mag2_t(t_*t_-x_*x_-y_*y_-z_*z_);
  };

  // Invariant mass
  xyzt_t mass() {
    #pragma HLS INLINE off
    return sqrt(mass2());
  };

  // Add two 4-vectors
  ap_lorentz<T> operator+(ap_lorentz<T> rhs) {
    return ap_lorentz(xyzt_t(x_+rhs.x_), xyzt_t(y_+rhs.y_), xyzt_t(z_+rhs.z_), xyzt_t(t_+rhs.t_));
  };


// semi-private:

  // Guard bits = log_2(cordic iterations ~ input width)
  typedef BitWidth<xyzt_t::width> gwidth;
  typedef ap_fixed<mag_t::width+gwidth::Value, mag_t::width> cordic_r_t;
  typedef ap_fixed<phi_t::width+gwidth::Value, phi_t::iwidth> cordic_phi_t;
  typedef ap_fixed<eta_t::width+gwidth::Value, eta_t::iwidth> cordic_eta_t;

  static void polar_to_cart(xyzt_t r_in, phi_t phi_in, xyzt_t& x_out, xyzt_t& y_out) {
    #pragma HLS PIPELINE II=1
    cordic_r_t x, y;
    cordic_phi_t phi;

    // std::cout << "Starting cordic" << std::endl;
    // std::cout << "r_in: " << r_in << ", phi_in: " << phi_in << std::endl;

    // pre-scale
    cordic_r_t r_tmp = r_in * ap_ufixed<cordic_r_t::width, 0>(ap_lorentz_cordic_scales[mag_t::width-2-1]);

    // Coarse rotate
    bool phi_pos = phi_in > 0;
    phi_t absphi_in = (phi_pos) ? phi_in : phi_t(-phi_in);
    if ( absphi_in <= phi_t(M_PI/4.) ) {
      // East quadrant
      x = r_tmp;
      y = 0.;
      phi = phi_in;
    }
    else if ( phi_in <= phi_t(3*M_PI/4.) ) {
      if ( phi_pos ) {
        // North, rotate by -pi/2
        x = 0.;
        y = r_tmp;
        phi = phi_in - phi_t(M_PI/2.);
      }
      else {
        // South, rotate by +pi/2
        x = 0.;
        y = -r_tmp;
        phi = phi_in + phi_t(M_PI/2.);
      }
    }
    else {
      // West, rotate by pi
      x = -r_tmp;
      y = 0.;
      phi = (phi_pos) ? cordic_phi_t(phi_in) - cordic_phi_t(M_PI) : cordic_phi_t(phi_in) + cordic_phi_t(M_PI);
    }

    // std::cout << "  Iter 0, x: " << x << ", y: " << y << ", phi: " << phi << std::endl;

    cordic_r_t xtmp, ytmp;
    for(size_t i=1; i<=mag_t::width-2; ++i) {
      xtmp = (phi>=0) ? x - (y>>i) : x + (y>>i);
      ytmp = (phi>=0) ? y + (x>>i) : y - (x>>i);
      phi = (phi>=0) ? phi-cordic_phi_t(ap_lorentz_cordic_angles[i-1]) : phi+cordic_phi_t(ap_lorentz_cordic_angles[i-1]);
      x = xtmp;
      y = ytmp;
      // std::cout << "  Iter " << i << ", x: " << x << ", y: " << y << ", phi: " << phi << std::endl;
    }

    x_out = x;
    y_out = y;
    // std::cout << "  Output x: " << x_out << ", y: " << y_out << std::endl;
  };

  static void cart_to_polar(xyzt_t x_in, xyzt_t y_in, mag_t& r_out, phi_t& phi_out) {
    #pragma HLS PIPELINE II=1
    cordic_r_t x, y;
    cordic_phi_t phi;

    // Coarse rotate to [-pi/4,pi/4)
    if ( x_in-y_in >= 0 ) {
      if ( x_in+y_in >= 0 ) {
        // East (Correct) quadrant
        x = x_in;
        y = y_in;
        phi = 0.;
      }
      else {
        // South, rotate by +pi/2
        x = -y_in;
        y = x_in;
        phi = M_PI/2.;
      }
    }
    else {
      if ( x_in+y_in >= 0 ) {
        // North, rotate by -pi/2
        x = y_in;
        y = -x_in;
        phi = -M_PI/2.;
      }
      else {
        // West, rotate by pi
        x = -x_in;
        y = -y_in;
        phi = M_PI;
      }
    }

    // std::cout << "Starting cordic" << std::endl;
    // std::cout << "  Iter 0, x: " << x << ", y: " << y << ", phi: " << phi << std::endl;

    cordic_r_t xtmp, ytmp;
    for(size_t i=1; i<=mag_t::width-2; ++i) {
      xtmp = (y>=0) ? x + (y>>i) : x - (y>>i);
      ytmp = (y>=0) ? y - (x>>i) : y + (x>>i);
      phi = (y>=0) ? phi-cordic_phi_t(ap_lorentz_cordic_angles[i-1]) : phi+cordic_phi_t(ap_lorentz_cordic_angles[i-1]);
      x = xtmp;
      y = ytmp;
      // std::cout << "  Iter " << i << ", x: " << x << ", y: " << y << ", phi: " << phi << std::endl;
    }

    r_out = x * ap_ufixed<mag_t::width, 0>(ap_lorentz_cordic_scales[mag_t::width-2-1]);
    phi_out = (phi > cordic_phi_t(M_PI)) ? cordic_phi_t(2*M_PI)-phi: -phi;
    // std::cout << "  Output r: " << r_out << ", phi: " << phi_out << std::endl;
  };

  static void hyperb_to_cart(xyzt_t r_in, eta_t eta_in, xyzt_t& x_out, xyzt_t& y_out) {
    #pragma HLS PIPELINE II=1
    // Prevent from using BRAM for such tiny arrays
    #pragma HLS ARRAY_PARTITION variable=ap_lorentz_cordic_xexact_hyp complete dim=0
    #pragma HLS ARRAY_PARTITION variable=ap_lorentz_cordic_yexact_hyp complete dim=0
    #pragma HLS ARRAY_PARTITION variable=ap_lorentz_cordic_xcoarse_hyp complete dim=0
    #pragma HLS ARRAY_PARTITION variable=ap_lorentz_cordic_ycoarse_hyp complete dim=0
    // std::cout << "Starting cordic" << std::endl;
    // std::cout << "r_in: " << r_in << ", eta_in: " << eta_in << std::endl;

    bool eta_pos = eta_in > 0;
    eta_t abseta = (eta_pos) ? eta_in : eta_t(-eta_in);
    cordic_eta_t eta = abseta & ~eta_t(7.5);
    ap_uint<4> ieta = abseta<<1;
    // std::cout << "abseta: " << abseta << ", eta: " << eta << ", ieta: " << ieta << std::endl;

    // sum(ap_lorentz_cordic_angles_hyp[1:]) = 0.2488, so the
    // ap_lorentz_cordic_xexact_hyp[0]=0.2554 step away from 0 or 0.5 will
    // never be corrected for.  This conditional keeps the numerical error low at i*0.5
    bool eta_zero = eta == 0;
    cordic_r_t x = r_in * ap_ufixed<cordic_r_t::width, 10>( (eta_zero) ? ap_lorentz_cordic_xexact_hyp[ieta] : ap_lorentz_cordic_xcoarse_hyp[mag_t::width-3-1][ieta]);
    cordic_r_t y = r_in * ap_ufixed<cordic_r_t::width, 10>( (eta_zero) ? ap_lorentz_cordic_yexact_hyp[ieta] : ap_lorentz_cordic_ycoarse_hyp[mag_t::width-3-1][ieta]);
    if ( not eta_zero ) {
      // std::cout << "  Scaled, x: " << x << ", y: " << y << ", eta: " << eta << std::endl;

      cordic_r_t xtmp, ytmp;
      for(size_t i=0; i<mag_t::width-3; ++i) {
        xtmp = (eta>=0) ? x + (y>>(i+2)) : x - (y>>(i+2));
        ytmp = (eta>=0) ? y + (x>>(i+2)) : y - (x>>(i+2));
        eta = (eta>=0) ? eta-cordic_eta_t(ap_lorentz_cordic_angles_hyp[i]) : eta+cordic_eta_t(ap_lorentz_cordic_angles_hyp[i]);
        x = xtmp;
        y = ytmp;
        // std::cout << "  Iter " << i << ", x: " << x << ", y: " << y << ", eta: " << eta << std::endl;
      }
    }

    x_out = x;
    y_out = (eta_pos) ? xyzt_t(y) : xyzt_t(-y);
    // std::cout << "  Output x: " << x_out << ", y: " << y_out << std::endl;
  };

  static xyzt_t sqrt(mag2_t xsq_in) {
    mag2_t x = 0;
    mag2_t xsq = xsq_in;
    mag2_t test = 0;
    test[mag2_t::width-1] = 1;
    // std::cout << "xsq: " << xsq << ", x: " << x << ", test: " << test << std::endl;
    for(size_t i=0; i<xyzt_t::iwidth; ++i) {
      bool gt = xsq >= (test+x);
      xsq = (gt) ? mag2_t(xsq - (test+x)) : xsq;
      x = (gt) ? mag2_t((x>>1) + test) : (x>>1);
      test = test >> 2;
      // std::cout << "xsq: " << xsq << ", x: " << x << ", test: " << test << std::endl;
    }
    test <<= 1;
    for(size_t i=0; i<xyzt_t::width-xyzt_t::iwidth; ++i) {
      bool gt = xsq >= (test+(x>>i));
      xsq = (gt) ? mag2_t(xsq - (test+(x>>i))) : xsq;
      x = (gt) ? mag2_t(x + test) : x;
      test = test >> 1;
      // std::cout << "*xsq: " << xsq << ", x: " << x << ", test: " << test << std::endl;
    }
    return xyzt_t(x);
  };

  static xyzt_t hls_sqrt(mag2_t xsq_in) {
    typename hls::sqrt_input<mag2_t::width, hls::CORDIC_FORMAT_USIG_INT>::in r2;
    r2.in = xsq_in;
    typename hls::sqrt_output<xyzt_t::width, hls::CORDIC_FORMAT_USIG_INT>::out r;
    hls::sqrt<hls::CORDIC_FORMAT_USIG_INT, mag2_t::width, xyzt_t::width, hls::CORDIC_ROUND_TRUNCATE>(r2, r);
    return r.out;
  };
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const ap_lorentz<T>& v) {
  os << "ap_lorentz(x=" << v.x_ << ", y=" << v.y_ << ", z=" << v.z_ << ", t=" << v.t_ << ")";
};

#endif // AP_LORENTZ_H
