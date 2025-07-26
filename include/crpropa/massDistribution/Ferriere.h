#pragma once

#include <cmath>
#include <numbers>
#include <string>
#include <sstream>

#include "crpropa/Common.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/massDistribution/TargetMedium.h"
#include "kiss/logger.h"


namespace crpropa {


/**
 * @addtogroup MassDistribution
 * @{
 */

/**
 * @class Ferriere
 * @brief model of the distribution of hydrogen in the Milky Way
 * The model for the outer galaxy (R > 3 kpc) is based on:
 * 	 "Global model of the interstellar medium in our Galaxy with new constraints on the hot gas component"
 *   K. Ferriere. 
 *   Astrophys. J. 497 (1998) 759.
 * The model for the inner galaxy (R < 3 kpc) is based on:
 *   "Spatial distribution of interstellar gas in the innermost 3 kpc of our galaxy"
 *   K. Ferriere, W. Gillard, P. Jean.
 *   Astron. Astrophys. 467 (2007) 611. 
 */
class Ferriere : public Density {
	private:
		static constexpr double Rsun = 8.5 * kpc;  // distance Solar System -- Galactic Centre
		bool containsHI;
		bool containsHII;
		bool containsH2;
		TargetMedium densityHI;
		TargetMedium densityHII;
		TargetMedium densityH2;

	public:
		Ferriere(bool HI = true, bool HII = true, bool H2 = true);
		void setContainsHI(bool HI);
		void setContainsHII(bool HII);
		void setContainsH2(bool H2);
		bool getContainsHI() const noexcept;
		bool getContainsHII() const noexcept;
		bool getContainsH2() const noexcept;
		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		[[nodiscard]] double getNucleonDensity(const Vector3d& position, const double& z = 0) const;
		[[nodiscard]] std::string getDescription() const override;

		[[nodiscard]] static double getDensityHI(const Vector3d& position);
		[[nodiscard]] static double getDensityHII(const Vector3d& position);
		[[nodiscard]] static double getDensityH2(const Vector3d& position);
		[[nodiscard]] static Vector3d transformCoordinatesCMZ(const Vector3d& position);
		[[nodiscard]] static Vector3d transformCoordinatesDisk(const Vector3d& position);
};

/** @} 
 */

}  // namespace crpropa


