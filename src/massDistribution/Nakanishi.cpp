#include "crpropa/massDistribution/Nakanishi.h"


namespace crpropa {


NakanishiHI::NakanishiHI() noexcept {
}

double NakanishiHI::getDensity(const Vector3d& position, const double& z) const {
	double R = hypot(position.x, position.y);  // radius in galactic plane
	double vertical = 1.06 * pc * (116.3 + 19.3 * R / kpc + 4.1 * pow_integer<2>(R / kpc) - 0.05 * pow_integer<3>(R / kpc));
	double plane = 0.94 / cm3 * (0.6 * exp(-R / (2.4 * kpc)) + 0.24 * exp(- pow_integer<2>((R - 9.5 * kpc) / (4.8 * kpc))));
	double n = plane * pow(0.5, pow_integer<2>(position.z / vertical)); // m^-3

	return n * pow_integer<3>(1 + z);
}

std::string NakanishiHI::getDescription() const {
	return "Nakanishi HI density distribution in the Milky Way.\n";
}

TargetMedium NakanishiHI::getTargetMedium() const {
	return TargetMedium("HI", "H", 1.0);
}



///////////////////////////////////////////////////////////////////////////////

NakanishiH2::NakanishiH2() noexcept {
}

double NakanishiH2::getDensity(const Vector3d& position, const double& z) const {
	double R = hypot(position.x, position.y);  // radius in galactic plane
	double vertical = 1.06 * pc * (10.8 * exp(0.28 * R / kpc) + 42.78);
	double plane = 0.94 / cm3 * (11.2 * exp(-R * R / (0.874 * kpc * kpc)) + 0.83 * exp(- pow_integer<2>((R - 4 * kpc) / (3.2 * kpc))));
	double n = plane * pow(0.5, pow_integer<2>(position.z / vertical)); // m^-3

	return n * pow_integer<3>(1 + z);
}

std::string NakanishiH2::getDescription() const {
	return "Nakanishi H2 density distribution in the Milky Way.";
}

TargetMedium NakanishiH2::getTargetMedium() const {
	return TargetMedium("H2", "H", 2.0);  // factor 2 for molecular hydrogen
}


}  // namespace crpropa
