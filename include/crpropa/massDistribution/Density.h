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
	public:
		virtual ~Density() = default;
		virtual double getDensity(const Vector3d& position, const double& z = 0) const = 0;
		virtual std::string getDescription() const {
			return "Density (abstract base class)";
		};
};


/**
 * @class DensityEvolution
 * @brief Density evolution model.
 * By default, a `Density` object evolves as (1+z)^3.
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
		ref_ptr<Density> getDensity() const;
		double getEvolutionIndex() const;
		double getDensity(const Vector3d& position, const double& z = 0) const override;
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
		ref_ptr<MediumComposition> target;

	public:
		DensityGrid(ref_ptr<MediumComposition> target, ref_ptr<Grid1f> grid);

		/** Get density at a given position.
		 * @param position position in Galactic coordinates with Earth at (-8.5 kpc, 0, 0)
		 * @returns Density in particles/m^3, sum up densities from added densities
		 */
		double getDensity(const Vector3d& position, const double& z = 0) const;
	
		/** Change the grid for the density
		 * @param grid (Grid1f) new grid for the density.
		 */
		void setGrid(ref_ptr<Grid1f> grid);

		/** Set the target field for the density
		 * @param target (MediumComposition) new target field for the density.
		 */
		void setTargetMedium(ref_ptr<MediumComposition> target);

		// double getNucleonDensity(const Vector3d& position) const;

		std::string getDescription() const;
};


/**
 * @class DensityList
 * @brief Superposition of density models.
 * The addDensity function adds a new density to the list.
 * The getDensity function handles the activated types in loaded densities, whereas get(type)Density disregards the activation state.
 * The targets are defined in each one of them.
 */
class DensityList: public Density {
	private:
		std::vector<ref_ptr<Density>> densityList ;

	public:
		/** Add new density to list.
		 * @param density density to add
		 */
		void addDensity(ref_ptr<Density> density);

		/** Get density at a given position.
		 * @param position position in a given coordinate system
		 * @param z redshift
		 * @returns Density in particles/m^3, sum up densities from added densities
		 */
		double getDensity(const Vector3d& position, const double& z = 0) const;

		std::string getDescription() const;

};

}  // namespace crpropa

#endif  // CRPROPA_DENSITY_H
