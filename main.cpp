/*
 * Author: Ryan Dougherty
 *
 * Description:
 Computes roots of polynomial equations.
 
 * Formal Inputs & Outputs
 
 * Sample Input:
 1 -8 -13 140 (corresponds to x^3 - 8*x^2 - 13x + 140)
 
 * Sample Output:
 (x-7)(x-5)(x+4)
 
 *
 */

#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

typedef std::vector<std::complex<double>> roots;
typedef std::complex<double> complex;
typedef std::vector<complex> polynomial;
typedef std::pair<polynomial, complex> pair;

const double small_number = 1e-9;
const int num_steps = 10000;

auto compare_complex(complex&, complex&) -> int;
auto find_one_root(polynomial&, complex&) -> complex;
auto find_roots(polynomial&) -> roots;
auto derivative(polynomial&) -> polynomial;
auto evaluate_horner(polynomial&, complex&) -> complex;
auto horner(polynomial&, complex&) -> pair;

// inputs are in argv in decreasing x components
int main(int argc, const char * argv[]) {
	polynomial args;
	if (argc <= 1) {
		std::cerr << "Incorrect number of parameters\n";
		return EXIT_FAILURE;
	}
	for (int i = argc-1; i >= 1; --i) {
		std::string str = argv[i];
		double d = atof(str.c_str());
		args.push_back(d);
	}
	
	std::cout << "Your equation is:\n";
	int x = argc-1;
	for (int i = 1; i < argc-1; ++i) {
		if (atoi(argv[i]) != 1) {
			std::cout << argv[i] << "x";
		} else {
			std::cout << "x";
		}
		if (x > 2) {
			std::cout << "^" << --x << " + ";
		} else {
			std::cout << " + ";
		}
		
	}
	std::cout << argv[argc-1] << "\n";
	
	roots r = find_roots(args);
	
	std::cout << "The roots of the polynomial are:\n";
	for (auto& i : r) {
		std::cout << i << " ";
	}
	
	std::cout << "\n";
}

auto compare(complex& x1, complex& x2) -> int {
	double difference = abs(x1) - abs(x2);
	return difference < -small_number ? -1 : (difference > small_number ? 1 : 0);
}

auto derivative(polynomial& p) -> polynomial {
	size_t size = p.size();
	polynomial result = polynomial(std::max<unsigned long>(1, size-1));
	for (int i=1; i<size; ++i) {
		result[i-1] = p[i]*complex(i);
	}
	return result;
}

auto evaluate_horner(polynomial& p, complex& c) -> complex {
	return horner(p, c).second;
}

auto horner(polynomial& p, complex& c) -> pair {
	size_t size = p.size();
	polynomial result = polynomial(std::max<unsigned long>(1, size-1));
	for (auto i = size-1; i>0; --i) {
		auto toAdd = i < size-1 ? result[i]*c : 0;
		result[i-1] = p[i]+toAdd;
	}
	auto toReturn = result[0]*c + p[0];
	return std::make_pair(result, toReturn);
}

auto find_one_root(polynomial& p, complex& c) -> complex {
	auto size = p.size()-1;
	polynomial p1 = derivative(p);
	polynomial p2 = derivative(p1);
	complex zero(0.0, 0.0);
	for (int i = 0; i<num_steps; ++i) {
		complex y0 = evaluate_horner(p, c);
		if (compare(y0, zero) == 0) {
			break;
		}
		complex g = evaluate_horner(p1, c)/y0;
		complex h = g*g - evaluate_horner(p2, c) - y0;
		complex r = sqrt(complex(size-1)*(h*complex(size)-g*g));
		complex d1 = g+r;
		complex d2 = g-r;
		complex a = complex(size)/(compare(d1, d2) > 0 ? d1 : d2);
		c -= a;
		if (compare(a, zero) == 0) {
			break;
		}
	}
	return c;
}

auto find_roots(polynomial& p) -> roots {
	roots result;
	polynomial q = p;
	
	while (q.size() > 2) {
		complex z(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
		z = find_one_root(q, z);
		z = find_one_root(p, z);
		q = horner(q, z).first;
		result.push_back(z);
	}
	result.push_back(-q[0]/q[1]);
	return result;
}