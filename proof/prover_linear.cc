#include <queue>
#include <iostream>
#include "alpha.h"
#include "configurations.h"
#include "search.h"

// Approximation ratio we want to prove
const double target_ratio = 0.87362;
// rounding function used as (c, A, B) where f(x) is A*max(0, x-B)
const Pairing rounding(0.86451, 0.0, 0.0);

// Discard configurations where abs(mu_i) > max_mu
const double max_mu = 1-1e-5;
// Discard configurations where rho > max_rho
const double max_rho = 1-1e-5;

// Maximum depth to search to, width of intervals at this depth is
// 2/(2^max_depth).
const int max_depth = 30;

int main(int argc, char **argv) {
  std::cout.precision(10);

  Conf init(I(-max_mu, max_mu), I(-max_mu, max_mu), I(-max_rho, max_rho));
  int res = search(init, rounding, target_ratio, max_depth);
  if (res < 1) {
	  std::cerr << "PROOF FAILED\n";
  }
  
  return 0;
  
}
