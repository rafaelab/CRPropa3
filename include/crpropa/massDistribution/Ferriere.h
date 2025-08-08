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
 * @class FerriereBase
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
class FerriereBase {
	protected:
		static constexpr double Rsun = 8.5 * kpc;  // distance Solar System -- Galactic Centre

	public:
		/**
		 * @brief Transform coordinates from the Galactic coordinate system to the CMZ coordinate system.
		 * @param position The position in Galactic coordinates.
		 * @return The transformed position in CMZ coordinates.
		 */
		static Vector3d transformCoordinatesCMZ(const Vector3d& position);

		/**
		 * @brief Transform coordinates from the Galactic coordinate system to the Disk coordinate system.
		 * @param position The position in Galactic coordinates.
		 * @return The transformed position in Disk coordinates.
		 */
		static Vector3d transformCoordinatesDisk(const Vector3d& position);

		/**
		 * @brief Return the target medium for the density distribution.
		 * @return The target medium.
		 */
		virtual TargetMedium getTargetMedium() const = 0;

};


/**
 * @class FerriereHI
 * @brief HI density distribution in the Milky Way.
 */
class FerriereHI : public Density, public FerriereBase {
	public:
		FerriereHI() noexcept;
		double getDensity(const Vector3d& position, const double& z = 0) const override;
		std::string getDescription() const override;
		TargetMedium getTargetMedium() const override;
};
	

/**
 * @class FerriereHII
 * @brief HII density distribution in the Milky Way.
 */
class FerriereHII : public Density, public FerriereBase {
	public:
		FerriereHII() noexcept;
		double getDensity(const Vector3d& position, const double& z = 0) const override;
		std::string getDescription() const override;
		TargetMedium getTargetMedium() const override;
};
	

/**
 * @class FerriereH2
 * @brief H2 density distribution in the Milky Way.
 */
class FerriereH2 : public Density, public FerriereBase {
	public:
		FerriereH2() noexcept;
		double getDensity(const Vector3d& position, const double& z = 0) const override;
		std::string getDescription() const override;
		TargetMedium getTargetMedium() const override;
};
	



/** @} 
 */

}  // namespace crpropa


