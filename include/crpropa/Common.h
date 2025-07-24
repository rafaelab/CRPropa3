#ifndef CRPROPA_COMMON_H
#define CRPROPA_COMMON_H


#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <span>
#include <string>
#include <vector>

#include "kiss/logger.h"



/**
 * @file
 * @brief Common helper functions
 */

namespace crpropa {
	
/**
 * \addtogroup Core
 * @{
 */

// Returns the full path to a CRPropa data file
std::string getDataPath(std::string filename);

// Returns the install prefix
std::string getInstallPrefix();

// Returns a certain digit from a given integer
inline int digit(const int& value, const int& d) {
	return (value % (d * 10)) / d;
}

// Return value xclip which is the closest to x, so that lower <= xclip <= upper
template <typename T>
T clip(const T& x, const T& lower, const T& upper) {
	return std::max(lower, std::min(x, upper));
}

// Perform linear interpolation on a set of n tabulated data points X[0 .. n-1] -> Y[0 .. n-1]
// Returns Y[0] if x < X[0] and Y[n-1] if x > X[n-1]
double interpolate(double x, std::span<const double> X, std::span<const double> Y);

// Perform bilinear interpolation on a set of (n,m) tabulated data points X[0 .. n-1], Y[0 .. m-1] -> Z[0.. n-1*m-1]
// Returns 0 if x < X[0] or x > X[n-1] or y < Y[0] or y > Y[m-1]
double interpolate2d(double x, double y,  std::span<const double> X,  std::span<const double> Y,  std::span<const double> Z);

// Perform linear interpolation on equidistant tabulated data
// Returns Y[0] if x < lo and Y[n-1] if x > hi
double interpolateEquidistant(double x, double lo, double hi, const  std::span<const double> Y);

// Find index of value in a sorted vector X that is closest to x
size_t closestIndex(double x, const std::span<const double> X);



// pow implementation as template for integer exponents pow_integer<2>(x)
// evaluates to x*x
template<unsigned int exponent>
constexpr double pow_integer(double base) {
	if constexpr (exponent == 0) {
		return 1.;
	} else if constexpr (exponent == 1) {
		return base;
	} else {
		double result = 1.;
		double current = base;
		unsigned int exp = exponent;
		while (exp > 0) {
			if (exp & 1) {
				result *= current;
			}
			current *= current;
			exp >>= 1;
		}
		return result;
	}
}

/** @}*/

// - input:  function over which to integrate, integration limits A and B
// - output: 8-points Gauß-Legendre integral
constexpr std::array<double, 8> X = {.0950125098, .2816035507, .4580167776, .6178762444, .7554044083, .8656312023, .9445750230, .9894009349};
constexpr std::array<double, 8> W = {.1894506104, .1826034150, .1691565193, .1495959888, .1246289712, .0951585116, .0622535239, .0271524594};
template<typename Integrand>
[[nodiscard]] constexpr double gaussInt(Integrand&& integrand, double A, double B) {
	const double XM = 0.5 * (B + A);
	const double XR = 0.5 * (B - A);
	double SS = 0.;

	for (size_t i = 0; i < X.size(); ++i) {
		double DX = XR * X[i];
		SS += W[i] * (integrand(XM + DX) + integrand(XM - DX));
	}
	return XR * SS;
}

} // namespace crpropa

#endif // CRPROPA_COMMON_H
