#ifndef AP_LORENTZ_H
#define AP_LORENTZ_H

#include <ostream>
#include "ap_fixed.h"
// BitWidth and Power template classes
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

  // Return transverse component
  // i.e. p_T or rho
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

  // Add two 4-vectors
  ap_lorentz<T> operator+(ap_lorentz<T> rhs) {
    return ap_lorentz(xyzt_t(x_+rhs.x_), xyzt_t(y_+rhs.y_), xyzt_t(z_+rhs.z_), xyzt_t(t_+rhs.t_));
  };


private:
  // Guard bits = log_2(cordic iterations ~ input width)
  typedef BitWidth<xyzt_t::width> gwidth;
  // Number of coarse steps for given eta range
  typedef Power<2, eta_t::iwidth-1> maxEta;
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

    cordic_r_t x = r_in * ap_ufixed<cordic_r_t::width, 10>( (eta==0) ? ap_lorentz_cordic_xexact_hyp[ieta] : ap_lorentz_cordic_xcoarse_hyp[mag_t::width-3-1][ieta]);
    cordic_r_t y = r_in * ap_ufixed<cordic_r_t::width, 10>( (eta==0) ? ap_lorentz_cordic_yexact_hyp[ieta] : ap_lorentz_cordic_ycoarse_hyp[mag_t::width-3-1][ieta]);
    if ( eta != 0 ) {
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

  // TODO: eta = asinh(pz/pt)
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const ap_lorentz<T>& v) {
  os << "ap_lorentz(x=" << v.x_ << ", y=" << v.y_ << ", z=" << v.z_ << ", t=" << v.t_ << ")";
};

#endif // AP_LORENTZ_H
