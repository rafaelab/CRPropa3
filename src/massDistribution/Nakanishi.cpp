#include "crpropa/massDistribution/Nakanishi.h"


namespace crpropa {

Nakanishi::Nakanishi(bool HI, bool H2) {
	if (not HI and not H2) 
		throw std::runtime_error("Nakanishi: At least one type of hydrogen must be activated.");

	if (HI)
		setContainsHI(HI);
	if (H2)
		setContainsH2(H2);
}

void Nakanishi::setContainsHI(bool HI) {
	containsHI = HI;
	densityHI = TargetMedium("HI", "H", 1.0);
}

void Nakanishi::setContainsH2(bool H2) {
	containsH2 = H2;
	densityH2 = TargetMedium("H2", "H", 2.0);
}

bool Nakanishi::getContainsHI() const noexcept {
	return containsHI;
}

bool Nakanishi::getContainsH2() const noexcept {
	return containsH2;
}

bool Nakanishi::getContainsHII() const noexcept {
	// Nakanishi does not model HII
	return false;  
}

double Nakanishi::getDensity(const Vector3d& position, const double& z) const {
	double totalDensity = 0;
	if (containsHI)
		totalDensity += getDensityHI(position);
	if (containsH2)
		totalDensity += getDensityH2(position);

	return totalDensity * pow_integer<3>(1 + z);
}

double Nakanishi::getNucleonDensity(const Vector3d& position, const double& z) const {
	double totalDensity = 0;

	if (containsHI)
		totalDensity += getDensityHI(position) * densityHI.getWeight();
	if (containsH2)
		totalDensity += getDensityH2(position) * densityH2.getWeight();

	return totalDensity * pow_integer<3>(1 + z);
}

std::string Nakanishi::getDescription() const {
	std::ostringstream oss;
	oss << "Nakanishi density distribution in the Milky Way:\n";
	oss << ". contains HI: " << (getContainsHI() ? "yes" : "no") << "\n";
	oss << ". contains HII: " << (getContainsHII() ? "yes" : "no") << "\n";
	oss << ". contains H2: " << (getContainsH2() ? "yes" : "no") << "\n";
	return oss.str();
}

double Nakanishi::getDensityHI(const Vector3d& position) {
	double R = hypot(position.x, position.y);  // radius in galactic plane
	double vertical = 1.06 * pc * (116.3 + 19.3 * R / kpc + 4.1 * pow_integer<2>(R / kpc) - 0.05 * pow_integer<3>(R / kpc));
	double plane = 0.94 / cm3 * (0.6 * exp(-R / (2.4 * kpc)) + 0.24 * exp(- pow_integer<2>((R - 9.5 * kpc) / (4.8 * kpc))));

	return plane * pow(0.5, pow_integer<2>(position.z / vertical)); // m^-3
}

double Nakanishi::getDensityH2(const Vector3d& position) {
	double R = hypot(position.x, position.y);  // radius in galactic plane
	double vertical = 1.06 * pc * (10.8 * exp(0.28 * R / kpc) + 42.78);
	double plane = 0.94 / cm3 * (11.2 * exp(-R * R / (0.874 * kpc * kpc)) + 0.83 * exp(- pow_integer<2>((R - 4 * kpc) / (3.2 * kpc))));

	return plane * pow(0.5, pow_integer<2>(position.z / vertical)); // m^-3
}


}  // namespace crpropa
