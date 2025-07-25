#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include <sstream>
#include <vector>

#include "crpropa/Grid.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/massDistribution/TargetMedium.h"
#include "HepPID/ParticleIDMethods.hh"
#include "kiss/logger.h"


namespace crpropa {

/**
 * @class Density
 * @brief Abstract base class for target densities.
 *
 */
class Density: public Referenced {
	protected:
		TargetMedium target;  ///< Target medium composition for the density

	public:
		Density() = default;
		virtual ~Density() = default;
		[[nodiscard]] virtual double getDensity(const Vector3d& position, const double& z = 0) const = 0;
		[[nodiscard]] virtual std::string getDescription() const {
			std::stringstream ss;
			ss << "Density (abstract base class)\n";
			ss << "Target Medium: " << target.getName() << " (Weight: " << target.getWeight() << ")";
			return ss.str();
		}

		void setTargetMedium(TargetMedium t) {
			target = std::move(t);
		}

		TargetMedium getTargetMedium() const {
			return target;
		}

		Density(const Density&) = default;
		Density& operator=(const Density&) = default;
		Density(Density&&) = default;
		Density& operator=(Density&&) = default;
};


/**
 * @class DensityEvolution
 * @brief Decorator for density evolution model.
 * By default, a `Density` object does not evolve in redshift (in reality, it is (1+z)^3).
 * Here a new index is introduced, which enables evolutions like (1+z)^(index).
 * Note that there is no default evolution index, so the user must set it explicitly.
 */
class DensityEvolution : public Density {
	protected:
		ref_ptr<Density> density;  ///< the density object to evolve
		double index;              ///< evolution index
		
	public:
		DensityEvolution(ref_ptr<Density> density, double index);
		void setDensity(ref_ptr<Density> density);
		void setEvolutionIndex(double index);
		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		std::string getDescription() const override;
};

/**
 * @class DensityGrid
 * @brief Wrapper to use a Grid1f for a density
 *
 * The DensityGrid uses a given grid for the chosen density type. More than one type can be chosen to follow the same distribution.
 * If no type is chosen a warning will be raised and all densities are 0.
 */
class DensityGrid: public Density {
	protected: 
		ref_ptr<Grid1f> grid; 

	public:
		DensityGrid(TargetMedium target, ref_ptr<Grid1f> grid);
		void setGrid(ref_ptr<Grid1f> grid);
		[[nodiscard]] double getDensity(const Vector3d& position, const double& z = 0) const override;
		std::string getDescription() const override;
};



}  // namespace crpropa

#endif  // CRPROPA_DENSITY_H
