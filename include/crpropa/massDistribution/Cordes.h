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
 * @class Cordes
 * @brief Cylindrical symetrical model of the density of ionised hydrogen (HII) of the Milky Way
 * Cordes et al., 1991, Nature 353,737
 */
class Cordes: public Density {
	private:
		TargetMedium densityHII;

	public:
		Cordes();
		constexpr bool getContainsHI() const noexcept;
		constexpr bool getContainsHII() const noexcept;
		constexpr bool getContainsH2() const noexcept;

		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		[[nodiscard]] double getNucleonDensity(const Vector3d& position, const double& z = 0) const;
		[[nodiscard]] std::string getDescription() const override;

		[[nodiscard]] static double getDensityHII(const Vector3d& position);
};


/** @} 
 */


}  // namespace crpropa

