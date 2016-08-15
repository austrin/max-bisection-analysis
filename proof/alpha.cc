#include "alpha.h"
#include "gaussian.h"


Pairing::Pairing(double c, double A, double B): c(c), A(A), B(B) {}
  
I Pairing::f(const I& x) const {
	assert(x >= 0.0);
	return A*max(0.0, x-B);
}

I Pairing::r(const I& mu) const {
	if (poslt(mu, 0.0)) {
		if (posgt(mu, 0.0)) {
			return hull(-r(I(0.0, -mu.lower())), r(I(0.0, mu.upper())));
		} else {
			return -r(-mu);
		}
	}
	// At this point, mu is guaranteed to be positive
	return I((c*mu).lower(), (c*mu + (1-c)*f(mu)).upper());
}



I alpha(const Conf& X, const I& r1, const I& r2) {
	I trh = X.trho();
	return 2.0*(1.0-Lambda(r1, r2, trh)) / (1.0 - X.rho);
}


I alpha_linear(const Conf& X, double c) {
	return alpha(X, c*X.mu1, c*X.mu2);
}


I alpha_pairing(const Conf &X, const Pairing &P) {
	I r1 = P.r(X.mu1);
	I r2 = P.r(X.mu2);
	I boost = (1-P.c)*P.f(min(abs(X.mu1), abs(X.mu2)));
	if (X.mu1.lower() >= 0.0 && X.mu2.upper() <= 0.0 && boost.lower() > 1e-15) {
		// Opposite signs and there is a boost, use stronger bounds
		// Note: there is no need to consider the case when mu1 <= 0 and
		// mu2 >= 0 as we assume the configuration to be in canonical
		// form, implying mu1 >= mu2.
		return hull(alpha(X, r1, I(r2.lower(), r2.upper() - boost.lower() + 1e-15)),
					alpha(X, I(r1.lower() + boost.lower() - 1e-15, r1.upper()), r2));
	} 
	// (Possibly) equal signs, or no boost, use weak bounds
	return alpha(X, r1, r2);
}
