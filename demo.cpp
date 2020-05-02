#include "cf.h"
#include <iostream>
#include <time.h>

// Computes pi, simple version using C++ operator overloading
void pi_simple(int Terms) {
	CF term(2, 1);
	CF pi(term);

	double start = clock();
	std::cout << "Computing series (" << Terms << " terms, no precision control)...";
	for (int n = 1; n < Terms; ++n) {
		term *= n;
		term /= 2 * n + 1;
		pi += term;
	}
	std::cout << std::endl;
	double end = clock() - start;

	std::cout << "CF of PI: " << pi.to_str() << std::endl;
	std::cout << "Took " << end / CLOCKS_PER_SEC << "s" << std::endl;
}

// Computes pi, precision control with c-style functions
void pi_precise(int Terms) {
	const unsigned long long Precision = Terms;
	CF term(2, 1);
	CF pi(term);

	double start = clock();
	std::cout << "Computing series (" << Terms << " terms, precision capped at " << Precision << " bits)...";
	for (int n = 1; n < Terms; ++n) {
		term = cf_mul_frac(term, n, (2 * n + 1), Precision);
		pi = cf_add(pi, term, Precision);
	}
	std::cout << std::endl;
	double end = clock() - start;

	std::cout << "CF of PI: " << pi.to_str() << std::endl;
	std::cout << "Took " << end / CLOCKS_PER_SEC << "s" << std::endl;
}

// Calculate m * arctan(x) with precision capped at "Precision" bits
CF arctan(long long int x, long long int m, int Precision) {
	// If splitting is computationally less expensive, do it
	if ((1. / log(2 * x) + 1. / log(x * (4 * x * x + 3))) < 1. / log(x))
		return cf_sub(arctan(2 * x, 2 * m, Precision), arctan(x * (4 * x * x + 3), m, Precision), Precision);

	std::cout << "ArcTan[1 / " << x << "]...\n";
	int Terms = Precision * (log(2) / log(x * x));

	std::vector<CF> series;
	series.resize(Terms);
	mpz_class X = x;
	series[0] = CF(1, x);

	double start = clock();
	std::cout << "Computing series [" << Terms << " terms, precision capped at " << Precision << " bits)...";
	for (int n = 1; n < Terms; ++n) {
		X *= x * x;
		series[n].add_term(0);
		series[n].add_term(X * (2 * n + 1));
	}
	std::cout << std::endl;

	std::cout << "Summation";

	// First pass, alternating terms are negative so subtract
	for (int n = 0, m = 0; n < series.size() - 1; n += 2, ++m)
		series[m] = cf_sub(series[n], series[n + 1], Precision);
	if (series.size() % 2 == 1) {
		series[series.size() / 2] = series[series.size() - 1];
		series.resize(series.size() / 2 + 1);
	}
	else series.resize(series.size() / 2);

	std::cout << '.' << std::flush;

	while (series.size() != 1) {
		for (int n = 0, m = 0; n < series.size() - 1; n += 2, ++m) {
			series[m] = cf_add(series[n], series[n + 1], Precision);
		}

		if (series.size() % 2 == 1) {
			series[series.size() / 2] = series[series.size() - 1];
			series.resize(series.size() / 2 + 1);
		}
		else series.resize(series.size() / 2);

		std::cout << '.' << std::flush;
	}
	double end = clock() - start;
	std::cout << "\nTook " << end / CLOCKS_PER_SEC << "s" << std::endl;

	return cf_mul_frac(series[0], m, 1, Precision);
}

// Calculate pi = 4 * arctan(1)
void pi_arctan_1(int Precision) {
	double start = clock();
	CF pi = arctan(1, 4, Precision);
	double end = clock();

	std::cout << "CF of PI: " << pi.to_str() << std::endl;
	std::cout << "Took " << (end - start) / CLOCKS_PER_SEC << "s" << std::endl;
}

// Calculate pi using pi = 176 * arctan(1 / 57) + 28 * arctan(1 / 239) - 48 * arctan(1 / 682) + 96 * arctan(1 / 12943)
void pi_arctan_fast(int Precision) {
	double start = clock();

	CF pi = cf_add(arctan(57, 176, Precision), arctan(12943, 96, Precision), Precision);
	pi = cf_add(pi, cf_sub(arctan(239, 28, Precision), arctan(682, 48, Precision), Precision), Precision);

	double end = clock();

	std::cout << "CF of PI: " << pi.to_str() << std::endl;
	std::cout << "Took " << (end - start) / CLOCKS_PER_SEC << "s" << std::endl;
}

int main() {
	pi_simple(1000);
	//pi_precise(1000);
	//pi_arctan_1(1000);
	//pi_arctan_fast(1000);

	return 0;
}
