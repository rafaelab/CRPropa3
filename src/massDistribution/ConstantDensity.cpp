#include "crpropa/massDistribution/ConstantDensity.h"



namespace crpropa {

ConstantDensity::ConstantDensity(double density, ref_ptr<MediumComposition> target) {
	Density::setTargetMedium(target);
	setDensityValue(density);
}

void ConstantDensity::setDensityValue(double density) {
	densityValue = density;
}

double ConstantDensity::getDensityValue() const {
	return densityValue;
}

double ConstantDensity::getDensity(const Vector3d& position, const double& z) const {
	return densityValue * target->getCompositionWeight();
}

std::string ConstantDensity::getDescription() const {
	std::stringstream ss;
	ss << "Constant density of " << densityValue << " particles/m^3";
	if (target) {
		ss << " for target: " << target->getDescription();
	} else {
		ss << " for no target";
	}
	return ss.str();
}


}  // namespace crpropa
