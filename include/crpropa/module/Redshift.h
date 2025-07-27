#pragma once

#include <limits>

#include "crpropa/Cosmology.h"
#include "crpropa/Module.h"
#include "crpropa/Units.h"

namespace crpropa {

/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 * @class Redshift
 * @brief Updates redshift and applies adiabatic energy loss according to the travelled distance.
 */
class Redshift: public Module {
	public:
		Redshift() = default;
		Redshift(const Redshift&) = default;
		Redshift& operator=(const Redshift&) = default;
		Redshift(Redshift&&) = default;
		Redshift& operator=(Redshift&&) = default;
		~Redshift() = default;
		void process(Candidate* candidate) const override;
		std::string getDescription() const override;
};

/**
 * @class FutureRedshift
 * @brief Updates redshift and applies adiabatic energy loss according to the traveled distance. 
 * Extends to negative redshift values to allow for symmetric time windows around z=0.
 */
class FutureRedshift: public Module {
	public:
		FutureRedshift() = default;
		FutureRedshift(const FutureRedshift&) = default;
		FutureRedshift& operator=(const FutureRedshift&) = default;
		FutureRedshift(FutureRedshift&&) = default;
		FutureRedshift& operator=(FutureRedshift&&) = default;
		~FutureRedshift() = default;
		void process(Candidate* candidate) const override;
		std::string getDescription() const override;
};

/** @}*/

} // namespace crpropa

