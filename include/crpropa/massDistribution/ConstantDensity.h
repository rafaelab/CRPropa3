#ifndef CRPROPA_CONSTANTDENSITY_H
#define CRPROPA_CONSTANTDENSITY_H


#include <cmath>
#include <sstream>
#include <string>

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/massDistribution/TargetMedium.h"
#include "kiss/logger.h"


namespace crpropa {

/**
 * @addtogroup MassDistribution
 * @{
 */

/**
 * @class ConstantDensity
 * @brief Set a constant density
*/
class ConstantDensity: public Density {
	protected:
		double densityValue;
		TargetMedium targetMedium;

	public:
		ConstantDensity(double value, TargetMedium target);
		void setDensityValue(double density) noexcept;
		double getDensityValue() const noexcept;
		void setTargetMedium(const TargetMedium& target);
		[[nodiscard]] TargetMedium getTargetMedium() const;
		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		[[nodiscard]] std::string getDescription() const override;
};

}  // namespace crpropa


/** @} 
 */


#endif  // CRPROPA_CONSTANTDENSITY_H



