#include "crpropa/massDistribution/Density.h"


namespace crpropa {


DensityEvolution::DensityEvolution(ref_ptr<Density> density, double index) {
	setDensity(density);
	setEvolutionIndex(index);
}

void DensityEvolution::setDensity(ref_ptr<Density> d) {
	density = d;
}

void DensityEvolution::setEvolutionIndex(double m) {
	index = m;
}

double DensityEvolution::getDensity(const Vector3d& position, const double& z) const {
	return density->getDensity(position, z) * std::pow(1 + z, index);
}

std::string DensityEvolution::getDescription() const {
	std::stringstream ss;
	ss << "Density evolution with index " << index << " for density: " << density->getDescription();
	return ss.str();
}


///////////////////////////////////////////////////////////////////////////////

DensityGrid::DensityGrid(ref_ptr<MediumComposition> field, ref_ptr<Grid1f> densityGrid) {
	Density::setTargetMedium(field);
	setGrid(densityGrid);
}

void DensityGrid::setGrid(ref_ptr<Grid1f> g) {
	grid = g;
}

double DensityGrid::getDensity(const Vector3d& position, const double& z) const {
	return grid->interpolate(position);
}


std::string DensityGrid::getDescription() const {
	std::stringstream ss;
	ss << "Target field density distribution in a uniform grid" << std::endl; 
	return ss.str();
}


}  // namespace crpropa

