#include "crpropa/massDistribution/Ferriere.h"


namespace crpropa {

	
Ferriere::Ferriere(bool HI, bool HII, bool H2) {
	if (not HI and not HII and not H2) 
		throw std::runtime_error("Ferrière: At least one type of hydrogen must be activated.");

	if (HI)
		setContainsHI(HI);
	if (HII)
		setContainsHII(HII);
	if (H2)
		setContainsH2(H2);
}

void Ferriere::setContainsHI(bool HI) {
	containsHI = HI;
	densityHI = TargetMedium("HI", "H", 1.0);
}

void Ferriere::setContainsHII(bool HII) {
	containsHII = HII;
	densityHII = TargetMedium("HII", "H", 1.0);
}

void Ferriere::setContainsH2(bool H2) {
	containsH2 = H2;
	containsH2 = H2;
	densityH2 = TargetMedium("H2", "H", 2.0);
}

bool Ferriere::getContainsHI() const noexcept {
	return containsHI;
}

bool Ferriere::getContainsHII() const noexcept {
	return containsHII;
}

bool Ferriere::getContainsH2() const noexcept {
	return containsH2;
}

double Ferriere::getDensity(const Vector3d& position, const double& z) const {
	double totalDensity = 0;
	if (containsHI)
		totalDensity += getDensityHI(position);
	if (containsHII)
		totalDensity += getDensityHII(position);
	if (containsH2)
		totalDensity += getDensityH2(position);

	return totalDensity * pow_integer<3>(1 + z);
}

double Ferriere::getNucleonDensity(const Vector3d& position, const double& z) const {
	double totalDensity = 0;
	if (containsHI)
		totalDensity += getDensityHI(position) * densityHI.getWeight();
	if (containsHII)
		totalDensity += getDensityHII(position) * densityHII.getWeight();
	if (containsH2)
		totalDensity += getDensityH2(position) * densityH2.getWeight();

	return totalDensity * pow_integer<3>(1 + z);
}

[[nodiscard]] std::string Ferriere::getDescription() const {
	std::stringstream ss;
	ss << "Ferriere density distribution in the Milky Way:\n";
	ss << ". contains HI: " << (containsHI ? "yes" : "no") << "\n";
	ss << ". contains HII: " << (containsHII ? "yes" : "no") << "\n";
	ss << ". contains H2: " << (containsH2 ? "yes" : "no") << "\n";
	return ss.str();
}

double Ferriere::getDensityHI(const Vector3d& position) {
	double n = 0;
	double R = hypot(position.x, position.y);

	if (R < 3 * kpc) {   // inner
		double x = position.x / pc;
		double y = position.y / pc;
		double z = position.z / pc;

		// warm interstellar matter
		double H = (x * x + pow_integer<2>(y + 10.)) / (145 * 145);
		double nWIM = exp(-H)* exp(- pow_integer<2>(z + 20.) / (26 * 26));
		nWIM += 0.009 * exp(-pow_integer<2>((R / pc - 3700) / (0.5 * 3700))) * 1 / pow_integer<2>(cosh(z / 140.));
		nWIM += 0.005 * cos(std::numbers::pi * R / pc * 0.5 / 17000) * 1 / pow_integer<2>(cosh(z / 950.));
		nWIM *= 8.0 / cm3;  // rescaling with density at center

		//very hot interstellar matter
		double alphaVH = 21. * deg;  // angle for very hot IM in radiant
		double cosA = cos(alphaVH);
		double sinA = sin(alphaVH);
		double etta = y * cosA + z * sinA;  // coordinate transformation for VHIM along major axis
		double chi = -y * sinA + z * cosA;

		double nVHIM = 0.29 / cm3 * exp(-((x * x + etta * etta) / (162. * 162.) + chi * chi / (90 * 90)));

		n = nWIM + nVHIM;
	} else {  // outer region
		double z = position.z;

		double nWarm = 0.0237 / cm3 * exp(-(R * R - Rsun * Rsun) / pow_integer<2>(37 * kpc)) * exp(-fabs(z) / (1 * kpc));
		nWarm += 0.0013 / cm3 * exp(-(pow_integer<2>(R - 4 * kpc) - pow_integer<2>(Rsun - 4 * kpc)) /pow_integer<2>(2 * kpc)) * exp(-fabs(z) / (150 * pc));

		double nHot = 0.12 * exp(-(R - Rsun) / (4.9 * kpc));
		nHot += 0.88 * exp(-(pow_integer<2>(R - 4.5 * kpc) - pow_integer<2>(Rsun - 4.5 * kpc)) / pow_integer<2>(2.9 * kpc));
		nHot *= pow(R / Rsun, -1.65);
		nHot *= exp(- fabs(z) / ((1500 * pc) * pow(R / Rsun, 1.65)));
		nHot *= 4.8e-4 / cm3;

		n = nWarm + nHot;
	}

	return n;
}

double Ferriere::getDensityHII(const Vector3d& position) {
	double n = 0;
	double R = hypot(position.x, position.y);

	if (R < 3 * kpc) {   // inner
		double x = position.x / pc;
		double y = position.y / pc;
		double z = position.z / pc;

		// warm interstellar matter
		double H = (x * x + pow_integer<2>(y + 10.)) / (145 * 145);
		double nWIM = exp(-H)* exp(-pow_integer<2>(z + 20.) / (26 * 26));
		nWIM += 0.009 * exp(-pow_integer<2>((R / pc - 3700) / (0.5 * 3700))) * 1 / pow_integer<2>(cosh(z / 140.));
		nWIM += 0.005 * cos(std::numbers::pi * R / pc * 0.5 / 17000) * 1 / pow_integer<2>(cosh(z / 950.));
		nWIM *= 8.0 / cm3;  // rescaling with density at center

		//very hot interstellar matter
		double alphaVH = 21. * deg;  // angle for very hot IM in radiant
		double cosA = cos(alphaVH);
		double sinA = sin(alphaVH);
		double etta = y * cosA + z * sinA;  // coordinate transformation for VHIM along major axis
		double chi = -y * sinA + z * cosA;

		double nVHIM = 0.29 / cm3 * exp(-((x * x + etta * etta) / (162. * 162.) + chi * chi / (90 * 90)));

		n = nWIM + nVHIM;
	} else {  // outer region
		double z = position.z;

		double nWarm = 0.0237 / cm3 * exp(-(R * R - Rsun * Rsun) / pow_integer<2>(37 * kpc)) * exp(-fabs(z) / (1 * kpc));
		nWarm += 0.0013 / cm3 * exp(-(pow_integer<2>(R - 4 * kpc) - pow_integer<2>(Rsun - 4 * kpc)) /pow_integer<2>(2 * kpc)) * exp(-fabs(z) / (150 * pc));

		double nHot = 0.12 * exp(-(R - Rsun) / (4.9 * kpc));
		nHot += 0.88 * exp(-(pow_integer<2>(R - 4.5 * kpc) - pow_integer<2>(Rsun - 4.5*kpc)) / pow_integer<2>(2.9 * kpc));
		nHot *= pow(R / Rsun, -1.65);
		nHot *= exp(-fabs(z) / ((1500 * pc) * pow(R / Rsun, 1.65)));
		nHot *= 4.8e-4 / cm3;

		n = nWarm + nHot;
	}
	return n;
}

double Ferriere::getDensityH2(const Vector3d& position) {
	double n = 0;
	double R = hypot(position.x, position.y);

	if(R < 3 * kpc) {
		// density at center
		Vector3d pos = transformCoordinatesCMZ(position);
		double x = pos.x / pc;
		double y = pos.y / pc;
		double z = pos.z / pc;

		double A = hypot(x, 2.5 * y, 2); // ellipticity
		double nCMZ = exp(-pow((A - 125.) / 137., 4)) * exp(- pow(z / 18., 2));
		nCMZ *= 150 / cm3;  // rescaling

		// density in disk
		pos = transformCoordinatesDisk(position);
		x = pos.x / pc;
		y = pos.y / pc;
		z = pos.z / pc;

		A = hypot(x, 3.1 * y);
		double nDisk = exp(- pow_integer<4>((A - 1200) / 438)) * exp(-pow_integer<2>(z / 42));
		nDisk *= 4.8 / ccm;  // rescaling

		n = nCMZ + nDisk;
	} else {  // outer region
		double z = position.z / pc;
		n = pow(R / Rsun, -0.58);
		n *= exp(- (pow_integer<2>(R - 4.5*kpc)-pow_integer<2>(Rsun - 4.5 * kpc)) / pow_integer<2>(2.9 * kpc));
		n *= exp(- pow_integer<2>(z / (81 * pow(R / Rsun, 0.58))));
		n *= 0.58 / cm3;  // rescaling
	}

	return n;
}

Vector3d Ferriere::transformCoordinatesCMZ(const Vector3d& position) {
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistent with JF12 implementation
	double x = -position.x;
	double y = -position.y;

	double xC = -50 * pc; // offset
	double yC = 50 * pc;
	double sinTc = sin(70. * deg);
	double cosTc = cos(70. * deg);

	Vector3d pos;
	pos.x = (x - xC) * cosTc + (y - yC) * sinTc;
	pos.y = -(x - xC) * sinTc + (y - yC) * cosTc;
	pos.z = position.z;

	return pos;
}

Vector3d Ferriere::transformCoordinatesDisk(const Vector3d& position) {
	// set galactocentric coordinate system with the Sun at (-8.5,0.,0.) instead of (8.5, 0, 0) to be consistent with JF12 implementation
	double x = -position.x;
	double y = - position.y;
	double z = position.z;

	double alphaD = 13.5 * deg;  // rotation around x-axis
	double sinAd = sin(alphaD);
	double cosAd = cos(alphaD);
	double betaD = 20. * deg;  // rotation around y'-axis
	double sinBd = sin(betaD);
	double cosBd = cos(betaD);
	double TettaD = 48.5 * deg;  // rotation around x"-axis
	double sinTd = sin(TettaD);
	double cosTd = cos(TettaD);

	Vector3d pos;

	pos.x = x * cosBd * cosTd - y * (sinAd * sinBd * cosTd - cosAd * sinTd) - z * (cosAd * sinBd * cosTd + sinAd * sinTd);

	pos.y =  -x * cosBd * sinTd;
	pos.y += y * (sinAd * sinBd * sinTd + cosAd * cosTd);
	pos.y += z * (cosAd * sinBd * sinTd - sinAd * cosTd);

	pos.z = x * sinBd;
	pos.z += y * sinAd * cosBd;
	pos.z += z * cosAd * cosBd;

	return pos;
}


}  // namespace crpropa 