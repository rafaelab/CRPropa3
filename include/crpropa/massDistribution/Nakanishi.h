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
 * @class NakanishiHI
 * @brief Cylindrical symmetrical model of the density distribution of the Milky Way for atomic (HI) hydrogen.
 * 	Model for HI arXiv:astro-ph/0304338
 *	Fit to the models given in arXiv:1607.07886
 * 	The scale height over the galactic plane of atomic hydrogen is fitted by a polynomial of degree 3.
 *  The scale height over the galactic plane of molecular hydrogen is fitted by an exponential function.
 */
class NakanishiHI : public Density {
	public:
		NakanishiHI() noexcept;
		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		[[nodiscard]] std::string getDescription() const override;
		[[nodiscard]] TargetMedium getTargetMedium() const;
};


/**
 * @class NakanishiH2
 * @brief Cylindrical symmetrical model of the density distribution of the Milky Way for molecular (H2) hydrogen.
 * 	Model for H2 arXiv:astro-ph/0610769
 *	Fit to the models given in arXiv:1607.07886
 * 	The scale height over the galactic plane of atomic hydrogen is fitted by a polynomial of degree 3.
 *  The scale height over the galactic plane of molecular hydrogen is fitted by an exponential function.
 */
class NakanishiH2 : public Density {
	public:
		NakanishiH2() noexcept;
		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		[[nodiscard]] std::string getDescription() const override;
		[[nodiscard]] TargetMedium getTargetMedium() const;
};



/** @} 
 */

}  // namespace crpropa




