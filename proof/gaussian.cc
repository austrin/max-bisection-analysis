#include "gaussian.h"
#include "bvnl.h"
#include <algorithm>
#include <boost/math/special_functions/erf.hpp>

const I epsilon_interval(-std::numeric_limits<double>::epsilon(),
						 std::numeric_limits<double>::epsilon());

// Add a certain number of epsilons as error to an interval
I nudge(const I& x, int epsilons) {
  I eps = (double)epsilons * epsilon_interval;
  return x * (1.0 + eps) + eps;
}

I erf(const I& x) { 
	// Built in erf should have an error of <= epsilon.
	// We overestimate with 10 epsilons.
	return nudge(I(erf(x.lower()), erf(x.upper())), 10);
}

I erf_inv(const I& x) { 
	// erf_inv has an error of <= 2 epsilons according to Boost docs.
	// We overestimate with 20 epsilons.
	return nudge(I(boost::math::erf_inv(x.lower()),
				   boost::math::erf_inv(x.upper())), 20);
}

const I sqrt2 = sqrt(I(2.0));

I Phi(const I& x) {
  return (1.0 + erf(x/sqrt2))/2.0;
}

I Phi_inv(const I& x) {
  return sqrt2*erf_inv(2.0*x-1.0);
}


I bvnl(const I& dh, const I& dk, const I& r) {
  if (empty(dh) || empty(dk) || empty(r)) return I::empty();
  // Assumes Fact: bvnl is monotone in all three parameters
  I ans = I(bvnl_down(dh.lower(), dk.lower(), r.lower()), 
			bvnl_up(dh.upper(), dk.upper(), r.upper()));
  return ans;
}

I Gamma(const I& q1, const I& q2, const I& rho) {
  return bvnl(Phi_inv(q1), Phi_inv(q2), rho);
}

double Gamma_up(double q1, double q2, double rho) {
	// Assumes Fact: Gamma is monotone in all three parameters
	return bvnl_up(Phi_inv(q1).upper(), Phi_inv(q2).upper(), rho);
}


// Naive implementation of Lambda_{\trho}(r1, r2)
// Has unnecessary loss of precision due to repeated occurrences of r1 and r2
I Lambda_naive(const I& r1, const I &r2, const I& trho) {
  return 2.0*Gamma((1.0-r1)/2.0, (1.0-r2)/2.0, trho) + (r1+r2)/2.0;
}


double Lambda_up(double r1, double r2, double trho) {
	// TODO assumes something about error
	return 2.0*Gamma_up((1.0-r1)/2.0 + 1e-15, (1.0-r2)/2.0 + 1e-15, trho) + (r1+r2)/2.0 + 1e-15;
}


// Upper bound on Lambda which gives a good approximation for trho close to 0.
// Precondition: trho.upper() >= 0
double Lambda_up_near_zero(const I& r1, const I& r2, const I& trho) {
	// Assumes Lemma 2.7, which implies Lambda_trho(r1, r2) <= (1+r1*r2)/2 + 4*|trho|
	return ((1.0+r1*r2) / 2.0 + 4.0 * trho).upper();
}

// The "g" function from Lemma 5.5 of the paper
I Lambda_g(const I& r, const I& trho) {
  return 1.0 - 2.0*Phi(Phi_inv((1.0-r)/2.0) / trho);
}

// More accurate implementation of Lambda_{\trho}(r1, r2).
// Uses Lemma 5.5 of paper which characterizes the extreme points of
// Lambda_{\trho}(I_1, I_2).  In fact for performance reasons we only
// use it for the upper bound, which is what we need a good estimate
// on in order to get a good lower bound on alpha.  For the lower
// bound on Lambda we just use the naive bound.
I Lambda_precise(const I& r1, const I &r2, const I& trho) {
  I ans = Lambda_naive(r1, r2, trho);

  double r1_lo = r1.lower(), r1_hi = r1.upper();
  double r2_lo = r2.lower(), r2_hi = r2.upper();

  // The four combinations of extreme points for r1, r2.
  // Assumes Fact: Lambda is monotone in trho
  double ub = std::max(std::max(Lambda_up(r1_lo, r2_lo, trho.upper()),
								Lambda_up(r1_lo, r2_hi, trho.upper())),
					   std::max(Lambda_up(r1_hi, r2_lo, trho.upper()),
								Lambda_up(r1_hi, r2_hi, trho.upper())));

  if (posgt(trho, 0.0)) {
	  // When trho is (possibly) positive, Lambda is convex and there
	  // are five more possibilities for the upper bound.
	  I z;
	  
	  // r1 at extreme point, r2 = g(r1)
	  z = hull(z, Lambda_naive(r1_lo, 
							   intersect(r2, Lambda_g(r1_lo, trho.upper())), 
							   trho.upper()));
	  z = hull(z, Lambda_naive(r1_hi, 
							   intersect(r2, Lambda_g(r1_hi, trho.upper())), 
							   trho.upper()));
	  
	  // r2 at extreme point, r1 = g(r2)
	  z = hull(z, Lambda_naive(intersect(r1, Lambda_g(r2_lo, trho.upper())),
							   r2_lo,
							   trho.upper()));
	  z = hull(z, Lambda_naive(intersect(r1, Lambda_g(r2_hi, trho.upper())),
							   r2_hi,
							   trho.upper()));
	  
	  // (0, 0)
	  if (poseq(r1, 0.0) && poseq(r2, 0.0))
		  z = hull(z, Lambda_up(0.0, 0.0, trho.upper()));
	  
	  // If trho is close to zero, the computation of the "g" function
	  // of Lemma 5.5 is quite unstable and sometimes gives poor
	  // bounds.  To safeguard against these cases we also use the
	  // upper bound provided by Lemma 2.7 which gives good bounds for
	  // trho close to 0.
	  ub = std::max(ub, std::min(z.upper(), Lambda_up_near_zero(r1, r2, trho)));
  }
  
  return intersect(ans, I(0.0, std::min(ub, 1.0)));
}


I Lambda(const I& r1, const I& r2, const I& trho) {
	return Lambda_precise(r1, r2, trho);
}
