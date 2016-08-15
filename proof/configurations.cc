#include "configurations.h"

Conf::Conf() {}

Conf::Conf(const I& mu1, const I& mu2, const I& rho): 
  mu1(mu1), mu2(mu2), rho(rho) { }


I trho_naive(const I& mu1, const I& mu2, const I& rho) {
	I num = (rho-mu1*mu2);
	I den = sqrt((1.0-square(mu1))*(1.0-square(mu2)));
	if (poseq(den, 0.0)) {
		// If the denominator is zero then trho is defined to be 0
		I ans = I(0.0);
		if (posgt(den, 0.0)) {
			// If the numerator can be positive and denominator
			// arbitrarily small, need to add trho = 1 as possibility.
			// Similarly if numerator can be negative.
			if (posgt(num, 0.0))
				ans = hull(ans, 1.0);
			if (poslt(num, 0.0))
				ans = hull(ans, -1.0);
		}
		return ans;
	} 
	// Denominator non-zero, just divide and truncate to [-1,1]
	return intersect(num/den, I(-1.0, 1.0));
}


I Conf::trho() const {
	I mu1_lo(mu1.lower()), mu1_hi(mu1.upper());
	I mu2_lo(mu2.lower()), mu2_hi(mu2.upper());
	// Assumes Fact: any extreme points for trho has one of mu1 and
	// mu2 at the endpoint of its interval.
	I ans = hull(hull(trho_naive(mu1_lo, mu2, rho),
					  trho_naive(mu1_hi, mu2, rho)),
				 hull(trho_naive(mu1, mu2_lo, rho),
					  trho_naive(mu1, mu2_hi, rho)));
	return ans;
}

void Conf::include(const Conf& X) {
	mu1 = hull(mu1, X.mu1);
	mu2 = hull(mu2, X.mu2);
	rho = hull(rho, X.rho);
}


void Conf::split(std::queue<Conf> &dest) const {
  std::pair<I, I> sub_mu1 = bisect(mu1);
  std::pair<I, I> sub_mu2 = bisect(mu2);
  std::pair<I, I> sub_rho = bisect(rho);
  Conf(sub_mu1.first, sub_mu2.first, sub_rho.first).add_to(dest);
  Conf(sub_mu1.first, sub_mu2.first, sub_rho.second).add_to(dest);
  Conf(sub_mu1.first, sub_mu2.second, sub_rho.first).add_to(dest);
  Conf(sub_mu1.first, sub_mu2.second, sub_rho.second).add_to(dest);
  Conf(sub_mu1.second, sub_mu2.first, sub_rho.first).add_to(dest);
  Conf(sub_mu1.second, sub_mu2.first, sub_rho.second).add_to(dest);
  Conf(sub_mu1.second, sub_mu2.second, sub_rho.first).add_to(dest);
  Conf(sub_mu1.second, sub_mu2.second, sub_rho.second).add_to(dest);
}

void Conf::add_to(std::queue<Conf> &dest) const {
	if (canonic() && valid()) dest.push(*this);
}


const double triangle_ineq[4][3] = {
  {-1, -1, -1},
  {-1,  1,  1},
  { 1, -1,  1},
  { 1,  1, -1}};

// Returns:
//  true: some part of the cube is inside the polytope
//  false: the cube is completely outside the polytope
bool Conf::valid() const {
	for (int i = 0; i < 4; ++i) {
		I sum = (triangle_ineq[i][0]*mu1 + triangle_ineq[i][1]*mu2 + triangle_ineq[i][2]*rho);
		if (sum.lower() > 1.0) {
			return false;
		}
	}
	return true;
}

bool Conf::canonic() const {
	if (mu1.upper() < mu2.lower()) return false;
	if (mu1.upper() < 0.0) return false;
	if ((-mu2).lower() > mu1.upper()) return false;
	return true;
}
