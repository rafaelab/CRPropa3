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
 * @class Nakanishi
 * @brief Cylindrical symetrical model of the density distribution of the Milky Way for atomic (HI) and molecular (H2) hydrogen
 * 	Model for HI arXiv:astro-ph/0304338
 *	Model for H2 arxiv:astro-ph/0610769
 *	fit of the models given in arXiv:1607.07886
 * 	The scale height over the galactic plane of atomic hydrogen is fitted by a polynomial of degree 3.
 *  The scale height over the galactic plane of molecular hydrogen is fitted by an exponential function.
 */
class Nakanishi: public Density {
	private:
		bool containsHI;
		bool containsH2;
		TargetMedium densityHI;
		TargetMedium densityH2;
		
	public:
		Nakanishi(bool HI = true, bool H2 = true);
		void setContainsHI(bool HI);
		void setContainsH2(bool H2);
		[[nodiscard]] bool getContainsHI() const noexcept;
		[[nodiscard]] bool getContainsH2() const noexcept;
		[[nodiscard]] bool getContainsHII() const noexcept;

		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		[[nodiscard]] double getNucleonDensity(const Vector3d& position, const double& z = 0) const;
		[[nodiscard]] std::string getDescription() const override;

		[[nodiscard]] static double getDensityHI(const Vector3d& position);
		[[nodiscard]] static double getDensityH2(const Vector3d& position);
};


/** @} 
 */

}  // namespace crpropa




