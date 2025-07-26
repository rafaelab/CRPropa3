#include "crpropa/massDistribution/ConstantDensity.h"



namespace crpropa {

ConstantDensity::ConstantDensity(double density, TargetMediumPtrS target) 
	: densityValue(density), targetMedium(target) {
}

void ConstantDensity::setDensityValue(double density) noexcept{
	densityValue = density;
}

double ConstantDensity::getDensityValue() const noexcept {
	return densityValue;
}

[[nodiscard]] double ConstantDensity::getDensity(const Vector3d& position, const double& z) const {
	return densityValue * targetMedium->getWeight();
}

[[nodiscard]] std::string ConstantDensity::getDescription() const {
	std::stringstream ss;
	ss << "Constant density of " << densityValue << " particles/m^3\n";
	ss << targetMedium->getDescription();
	return ss.str();
}


}  // namespace crpropa
