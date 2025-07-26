#include "crpropa/massDistribution/Cordes.h"


namespace crpropa {


Cordes::Cordes() {
	densityHII = TargetMedium("HII", "H", 1.0);
}

constexpr bool Cordes::getContainsHI() const noexcept {
	return false;
}

constexpr bool Cordes::getContainsHII() const noexcept {
	return true;
}

constexpr bool Cordes::getContainsH2() const noexcept {
	return false;
}

double Cordes::getDensity(const Vector3d& position, const double& z) const {
	double n = getDensityHII(position);
	return n * pow_integer<3>(1 + z);
}

double Cordes::getNucleonDensity(const Vector3d& position, const double& z) const {
	return getDensityHII(position) * densityHII.getWeight() * pow_integer<3>(1 + z);
}

std::string Cordes::getDescription() const {
	std::ostringstream oss;
	oss << "Cordes density distribution in the Milky Way:\n";
	oss << ". contains HI: " << (getContainsHI() ? "yes" : "no") << "\n";
	oss << ". contains HII: " << (getContainsHII() ? "yes" : "no") << "\n";
	oss << ". contains H2: " << (getContainsH2() ? "yes" : "no") << "\n";
	return oss.str();
}

double Cordes::getDensityHII(const Vector3d& position) {
	double n = 0;  // in m^-3

	double z = position.z;
	double R = hypot(position.x, position.y); 

	 // galactocentric component
	n += 0.025 / cm3 * exp(-fabs(z) / (1 * kpc)) * exp(- pow_integer<2>(R / (20 * kpc))); 

	 // anular component
	n += 0.2 / cm3 * exp(-fabs(z) / (0.15 * kpc)) * exp(- pow_integer<2>((R - 4 * kpc) / (2 * kpc))); 

	return n;
}


}  // namespace crpropa
