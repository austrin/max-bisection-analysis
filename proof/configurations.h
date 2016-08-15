/* A class for dealing with configurations
 */
#ifndef __CONFIGURATIONS_H
#define __CONFIGURATIONS_H

#include "intervals.h"
#include <queue>

class Conf {
public:
	I mu1, mu2, rho;
	Conf();
	Conf(const I& mu1, const I& mu2, const I& rho);
  
	// Computes \tilde{rho} = (rho - mu1*mu2)/sqrt((1-mu1^2)*(1-mu2^2))
	I trho() const;
  
	// Sets this cube to the smallest cube containing itself and X
	void include(const Conf& X);
  
	// Split this configuration into 8 equal sub-cubes and
	// add those subcubes which are valid to dest
	void split(std::queue<Conf> &dest) const;

	// Check if a configuration satisfies triangle inequalities.
	// Return true if some part of the cube is inside the polytope.
	// The function only checks each triangle inequality individually
	// and therefore has false positives, but that's OK (whereas false
	// negatives would not be acceptable).  The false positives turn
	// out to be irrelevant in practice.
	bool valid() const;

	// Check if a configuration is in canonic form.  If not, it can be
	// discarded.  (mu1, mu2, rho) is in canonic form if
	// (1) mu >= 0 and
	// (2) mu1 >= mu2 and
	// (3) mu1 >= -mu2.
	// A configuration which is not in canonic form can
	// be transformed into one using the symmetries:
	// (mu1, mu2, rho) is equivalent to (mu2, mu1, rho)
	// and to (-mu1, -mu2, rho).
	// The function returns true iff some configuration in the cube is
	// canonic.
	bool canonic() const;

private:
	void add_to(std::queue<Conf> &dest) const;
};

#endif
