#include "crpropa/massDistribution/Cordes.h"


namespace crpropa {


CordesHII::CordesHII() noexcept {
	densityHII = TargetMedium("HII", "H", 1.0);
}

double CordesHII::getDensity(const Vector3d& position, const double& z) const {
	double n = 0;  // in m^-3

	double h = position.z;
	double R = hypot(position.x, position.y); 

	 // galactocentric component
	n += 0.025 / cm3 * exp(-fabs(h) / (1 * kpc)) * exp(- pow_integer<2>(R / (20 * kpc)));

	 // anular component
	n += 0.2 / cm3 * exp(-fabs(h) / (0.15 * kpc)) * exp(- pow_integer<2>((R - 4 * kpc) / (2 * kpc))); 

	return n * pow_integer<3>(1 + z);
}

std::string CordesHII::getDescription() const {
	return "CordesHII density distribution in the Milky Way.";
}

TargetMedium CordesHII::getTargetMedium() const {
	return TargetMedium("HII", "H", 1.0);
}

}  // namespace crpropa
