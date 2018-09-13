#ifndef AP_LORENTZ_H
#define AP_LORENTZ_H

#include "ap_fixed.h"
// BitWidth and Power template classes
#include "hls/utils/x_hls_utils.h"

const double ap_lorentz_cordic_angles[16] = {
  // python
  // >> import math
  // >> [math.atan(pow(2, -x-1)) for x in range(16)]
  0.46364760900080615, 0.24497866312686414, 0.12435499454676144, 0.06241880999595735, 0.031239833430268277, 0.015623728620476831, 0.007812341060101111, 0.0039062301319669718, 0.0019531225164788188, 0.0009765621895593195, 0.0004882812111948983, 0.00024414062014936177, 0.00012207031189367021, 6.103515617420877e-05, 3.0517578115526096e-05, 1.5258789061315762e-05
};
const double ap_lorentz_cordic_scales[16] = {
  // python
  // >> from operator import mul
  // >> a = [pow(1+pow(2, -2*(x+1)), -.5) for x in range(16)]
  // >> [reduce(mul, a[:n+1]) for n in range(16)]
  0.8944271909999159, 0.8677218312746247, 0.8610211763152799, 0.8593444051498349, 0.8589251104651273, 0.8588202804030459, 0.8587940724877404, 0.8587875204839219, 0.8587858824814052, 0.8587854729806784, 0.8587853706054906, 0.8587853450116933, 0.858785338613244, 0.8587853370136317, 0.8587853366137286, 0.8587853365137528
};
const double ap_lorentz_cordic_angles_hyp[16] = {
  // python
  // >> import math
  // >> [math.atanh(pow(2, -x-1)) for x in range(16)]
  0.5493061443340549, 0.2554128118829953, 0.12565721414045306, 0.06258157147700301, 0.03126017849066699, 0.01562627175205221, 0.007812658951540421, 0.003906269868396826, 0.0019531274835325498, 0.000976562810441036, 0.0004882812888051128, 0.0002441406298506386, 0.00012207031310632982, 6.103515632579122e-05, 3.05175781344739e-05, 1.5258789063684237e-05
};
const double ap_lorentz_cordic_scales_hyp[16] = {
  // python
  // >> from operator import mul
  // >> a = [pow(1-pow(2, -2*(x+1)), -.5) for x in range(16)]
  // >> [reduce(mul, a[:n+1]) for n in range(16)]
  1.1547005383792515, 1.1925695879998877, 1.2019971622805568, 1.204351713336805, 1.2049402067573813, 1.205087321122957, 1.205124099153, 1.205133293625434, 1.2051355922413503, 1.2051361668951923, 1.2051363105586443, 1.2051363464745068, 1.2051363554534724, 1.2051363576982137, 1.205136358259399, 1.2051363583996955
};

// T must be ap_fixed type
template<typename T>
class ap_lorentz {
public:
typedef T xyzt_t;
  // 3-magnitude at most sqrt(3) larger than max xyzt_t
  typedef ap_fixed<xyzt_t::width+1, xyzt_t::width+1> mag_t;
  typedef ap_fixed<xyzt_t::width, 4> eta_t;
  typedef ap_fixed<xyzt_t::width, 3> phi_t;

  xyzt_t x_, y_, z_, t_;

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
    t_ = (E==0.) ? mag : E;
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

private:
  // Guard bits = log(cordic iterations ~ input width)
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
    bool phi_pos = phi_in > 0.;
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
    #pragma HLS PIPELINE II=2
    std::cout << "Starting cordic" << std::endl;
    std::cout << "r_in: " << r_in << ", eta_in: " << eta_in << std::endl;

    cordic_r_t x = r_in;
    cordic_r_t y = 0.;
    cordic_eta_t eta = eta_in;
    
    std::cout << "  Iter 0, x: " << x << ", y: " << y << ", eta: " << eta << std::endl;

    cordic_r_t xtmp, ytmp;
    // Step 1: steps of atanh(0.5) to guarantee |eta|<0.5
    for(size_t i=0; i<size_t(1.82*maxEta::Value)-1; ++i) {
      xtmp = (eta>=0) ? x + (y>>1) : x - (y>>1);
      ytmp = (eta>=0) ? y + (x>>1) : y - (x>>1);
      eta = (eta>=0) ? eta-cordic_eta_t(ap_lorentz_cordic_angles_hyp[0]) : eta+cordic_eta_t(ap_lorentz_cordic_angles_hyp[0]);
      x = xtmp;
      y = ytmp;
      std::cout << "  Step iter " << i << ", x: " << x << ", y: " << y << ", eta: " << eta << std::endl;
    }
    // Step 2: binary search
    for(size_t i=1; i<=mag_t::width-2; ++i) {
      xtmp = (eta>=0) ? x + (y>>i) : x - (y>>i);
      ytmp = (eta>=0) ? y + (x>>i) : y - (x>>i);
      eta = (eta>=0) ? eta-cordic_eta_t(ap_lorentz_cordic_angles_hyp[i-1]) : eta+cordic_eta_t(ap_lorentz_cordic_angles_hyp[i-1]);
      x = xtmp;
      y = ytmp;
      std::cout << "  Search iter " << i << ", x: " << x << ", y: " << y << ", eta: " << eta << std::endl;
    }

    x_out = x * ap_ufixed<mag_t::width, 1>(ap_lorentz_cordic_scales_hyp[mag_t::width-2]);
    y_out = y * ap_ufixed<mag_t::width, 1>(ap_lorentz_cordic_scales_hyp[mag_t::width-2]);
    std::cout << "  Output x: " << x_out << ", y: " << y_out << std::endl;
  };

  // TODO: eta = asinh(pz/pt)
};

#endif // AP_LORENTZ_H
