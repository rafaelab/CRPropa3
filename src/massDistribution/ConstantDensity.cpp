#include "crpropa/massDistribution/ConstantDensity.h"



namespace crpropa {

ConstantDensity::ConstantDensity(double HI, double HII, double H2, double He) {
	// set all types active which are not equal 0 and change number density
	if (HI != 0)
		setHI(true, HI);
	if (HII != 0)
		setHII(true, HII);
	if (H2 != 0)
		setH2(true, H2);
	if (He != 0)
		setHe(true, He);
}

double ConstantDensity::getDensity(const Vector3d& position) const {
	double n = 0;

	if (isHI)
		n += densityHI;
	if (isHII)
		n += densityHII;
	if (isH2)
		n += densityH2;
	if (isHe)
		n += densityHe;

	// check if all densities are deactivated and raise warning if so
	if ((isHI or isHII or isH2 or isHe) == false) {
		KISS_LOG_WARNING << "Called getNucleonDensity on fully deactivated ConstantDensity gas density model. In this case the density is always set to 0." << std::endl;
	}

	return n;
}

double ConstantDensity::getNucleonDensity(const Vector3d& position) const {
	double n = 0;

	if (isHI)
		n += densityHI;
	if (isHII)
		n += densityHII;
	if (isH2)
		n += 2 * densityH2;
	if (isHe)
		n += 4 * densityHe; // Helium has 4 nucleons

	// check if all densities are deactivated and raise warning if so
	if ((isHI or isHII or isH2 or isHe) == false) {
		KISS_LOG_WARNING << "Called getNucleonDensity on fully deactivated ConstantDensity gas density model. In this case the density is always set to 0." << std::endl;
	}

	return n;
}

double ConstantDensity::getHIDensity(const Vector3d& position) const {
	return densityHI;
}

double ConstantDensity::getHIIDensity(const Vector3d& position) const{
	return densityHII;
}

double ConstantDensity::getH2Density(const Vector3d& position) const{
	return densityH2;
}

bool ConstantDensity::isActiveHI() {
	return isHI;
}

bool ConstantDensity::isActiveHII() {
	return isHII;
}

bool ConstantDensity::isActiveH2() {
	return isH2;
}

void ConstantDensity::setHI(bool activate, double density) {
	isHI = activate;
	densityHI = density;
}

void ConstantDensity::setHI(bool activate) {
	setHI(activate, densityHI);
}

void ConstantDensity::setHI(double density) {
	setHI(isHI, density);
}

void ConstantDensity::setHII(bool activate, double density) {
	isHII = activate;
	densityHII = density;
}

void ConstantDensity::setHII(bool activate) {
	setHII(activate, densityHII);
}

void ConstantDensity::setHII(double density) {
	setHII(isHII, density);
}

void ConstantDensity::setH2(bool activate, double density) {
	isH2 = activate;
	densityH2 = density;
}

void ConstantDensity::setH2(bool activate) {
	setH2(activate, densityH2);
}

void ConstantDensity::setH2(double density) {
	setH2(isH2, density);
}

void ConstantDensity::setHe(bool activate, double density) {
	isHe = activate;
	densityHe = density;
}

void ConstantDensity::setHe(bool activate) {
	setHe(activate, densityHe);
}

void ConstantDensity::setHe(double density) {
	setHe(isHe, density);
}

std::string ConstantDensity::getDescription() {
	std::stringstream s;
	s << "ConstantDensity:" << std::endl;

	s << "HI component is ";
	if (! isHI)
		s << "not ";
	s << "active and has a density of " << densityHI << " m^-3" << std::endl; 

	s << "HII component is ";
	if (! isHII)
		s << "not ";
	s << "active and has a density of " << densityHII << " m^-3" << std::endl;

	s << "H2 component is ";
	if (! isH2)
		s << "not ";
	s << "active and has a density of " << densityH2 << " m^-3" << std::endl;

	s << "He component is ";
	if (! isHe)
		s << "not ";
	s << "active and has a density of " << densityHe << " m^-3" << std::endl;

	return s.str();
}


}  // namespace crpropa
