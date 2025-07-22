#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include <sstream>

#include "crpropa/Grid.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/massDistribution/MediumComposition.h"
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
		ref_ptr<MediumComposition> target;  ///< target medium composition for the density

	public:
		virtual ~Density() = default;

		virtual double getDensity(const Vector3d& position, const double& z = 0) const = 0;

		virtual std::string getDescription() const {
			return "Density (abstract base class)";
		}

		void setTargetMedium(ref_ptr<MediumComposition> t) {
			target = t;
		}

		ref_ptr<MediumComposition> getTargetMedium() const {
			return target;
		}

		bool isAdmixed() const {
			return target->isAdmixed();
		}
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
		double getDensity(const Vector3d& position, const double& z = 0) const;
		std::string getDescription() const;
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
		DensityGrid(ref_ptr<MediumComposition> target, ref_ptr<Grid1f> grid);
		void setGrid(ref_ptr<Grid1f> grid);
		double getDensity(const Vector3d& position, const double& z = 0) const;
		std::string getDescription() const;
};



}  // namespace crpropa

#endif  // CRPROPA_DENSITY_H
