#pragma once

#include <cmath>
#include <string>
#include <sstream>

#include "crpropa/Common.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/massDistribution/TargetMedium.h"


namespace crpropa {

/**
 * @addtogroup MassDistribution
 * @{
 */


/**
 * @class CordesHII
 * @brief Cylindrical symetrical model of the density of ionised hydrogen (HII) of the Milky Way
 * Cordes et al., 1991, Nature 353,737
 */
class CordesHII: public Density {
	private:
		TargetMedium densityHII;

	public:
		CordesHII() noexcept;
		double getDensity(const Vector3d& position, const double& z = 0) const override;
		std::string getDescription() const override;
		TargetMedium getTargetMedium() const;
};


/** @} 
 */


}  // namespace crpropa

