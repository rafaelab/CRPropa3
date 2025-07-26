#include "crpropa/massDistribution/Density.h"


namespace crpropa {


///////////////////////////////////////////////////////////////////////////////

DensityEvolution::DensityEvolution(ref_ptr<Density> density, double index) 
	: density(density), index(index) {
}

void DensityEvolution::setDensity(ref_ptr<Density> d) {
	density = d;
}

void DensityEvolution::setEvolutionIndex(double m) {
	index = m;
}

[[nodiscard]] double DensityEvolution::getDensity(const Vector3d& position, const double& z) const {
	return density->getDensity(position, z) * std::pow(1 + z, index);
}

[[nodiscard]] std::string DensityEvolution::getDescription() const {
	return std::format("Density evolution with index {} for density: \n{}", index, density->getDescription());
}


///////////////////////////////////////////////////////////////////////////////

DensityGrid::DensityGrid( ref_ptr<Grid1f> densityGrid)
	: grid(densityGrid) {
}

void DensityGrid::setGrid(ref_ptr<Grid1f> g) {
	grid = g;
}

[[nodiscard]] double DensityGrid::getDensity(const Vector3d& position, const double& z) const {
	return grid->interpolate(position);
}

[[nodiscard]] std::string DensityGrid::getDescription() const {
	// return std::format("Density field density distribution in a uniform grid: \n{}", grid->getDescription());
	return "Density field density distribution in a uniform grid";
}



}  // namespace crpropa

