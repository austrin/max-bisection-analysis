#include <iostream>
#include "search.h"

const size_t MAX_QSIZE = 100000000;

int search(const Conf& initbox, const Pairing &rounding, 
		   double target_ratio, int max_depth) {
	std::queue<Conf> Q[2];
	Q[0].push(initbox);
	std::cerr << "Searching in box (" << initbox.mu1 << ", " << initbox.mu2 << ", " << initbox.rho << ")\n";

	double level_lb = 0;
	double worst_lb = 1;
	double worst_ub = 1;
	int tot_inspected = 0;
	int inconclusive = 0;

	for (int level = 0; !Q[level % 2].empty(); ++level) {
		int cur_level_inspected = 0;
		int cur_level_size = Q[level % 2].size();
		level_lb = 1;
		Conf curbox = Conf(I(), I(), I());

		while (!Q[level % 2].empty()) {
			Conf X = Q[level % 2].front();
			Q[level % 2].pop();

			I a = alpha_pairing(X, rounding);
			if (a.upper() < target_ratio) {
				std::cerr << "ERROR!\n";
				std::cerr << "Case:\n  " << X.mu1 << "\n  " << X.mu2 << "\n  " << X.rho << "\n";
				std::cerr << "Alpha: " << a.lower() << " up to " << a.upper() << "\n";
				return -1;
			}

			worst_ub = std::min(worst_ub, a.upper());
			level_lb = std::min(level_lb, a.lower());
			++tot_inspected;
			if (++cur_level_inspected % 50000 == 0) {
				std::cerr << "Level " << level << ": " << cur_level_size << " cases, " << (100.0*cur_level_inspected/cur_level_size) << "% done. Next level: " << Q[(level+1)%2].size() << ".          \r";
			}
			curbox.include(X);
		
			
			if (a.lower() >= target_ratio) {
				worst_lb = std::min(worst_lb, a.lower());
				// ==> Lower bound strong enough, done!
				continue;
			}
			if (level == max_depth || Q[(level+1)%2].size() > MAX_QSIZE) {
				std::cout << "Inconclusive case at:\n  (" << X.mu1 << ", " << X.mu2 << ", " << X.rho << ")\n";
				I tr = X.trho();
				std::cout << "  trho: " << tr << "\n";
				std::cout << "  alph: " << a << "\n";
				++inconclusive;
				continue;
			}
			X.split(Q[(level+1)%2]);
		}
		std::cerr << "Level " << level << ": " << cur_level_size << " cases. ";
		std::cerr << "LB: " << worst_lb << ". UB: " << worst_ub << ". ";
		std::cerr << "Bbox: (" << curbox.mu1 << ", " << curbox.mu2 << ", " << curbox.rho << ")\n";
	}
	
	std::cerr << "Search complete. ";
	std::cerr << "A total of " << tot_inspected << " cases were inspected\n";

	if (inconclusive) {
		std::cerr << "Of these, " << inconclusive << " were inconclusive\n";
		return 0;
	}
	std::cerr << "LB: " << worst_lb << ". UB: " << worst_ub << ".\n";
	std::cerr << "\n";
	return 1;
}

