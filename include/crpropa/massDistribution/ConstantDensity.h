#ifndef CRPROPA_CONSTANTDENSITY_H
#define CRPROPA_CONSTANTDENSITY_H


#include <cmath>
#include <sstream>
#include <string>

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/massDistribution/MediumComposition.h"
#include "kiss/logger.h"


namespace crpropa {

/**
 * @class ConstantDensity
 * @brief Set a constant density
*/
class ConstantDensity: public Density {
	protected:
		ref_ptr<MediumComposition> target;
		double densityValue;

	public:
		ConstantDensity(double value, ref_ptr<MediumComposition> target);
		void setDensityValue(double density);
		void setTargetMedium(ref_ptr<MediumComposition> target);
		double getDensityValue() const;
		ref_ptr<MediumComposition> getTargetMedium() const;
		double getDensity(const Vector3d& position, const double& z = 0) const;
		std::string getDescription() const;
};

}  // namespace crpropa

#endif  // CRPROPA_CONSTANTDENSITY_H



