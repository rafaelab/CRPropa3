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

ref_ptr<Density> DensityEvolution::getDensity() const {
	return density;
}

double DensityEvolution::getEvolutionIndex() const {
	return index;
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
	setTargetMedium(field);
	setGrid(densityGrid);
}

double DensityGrid::getDensity(const Vector3d& position, const double& z) const {
	return grid->interpolate(position);
}

void DensityGrid::setGrid(ref_ptr<Grid1f> g) {
	grid = g;
}

void DensityGrid::setTargetMedium(ref_ptr<MediumComposition> t) {
	target = t;
}


std::string DensityGrid::getDescription() const {
	std::stringstream ss;
	ss << "Target field density distribution in a uniform grid" << std::endl; 
	return ss.str();
}


///////////////////////////////////////////////////////////////////////////////

void DensityList::addDensity(ref_ptr<Density> density) {
	densityList.push_back(density);
}

double DensityList::getDensity(const Vector3d& position, const double& z) const {
	double n = 0.;
	for (int i = 0; i < densityList.size(); i++)
		n += densityList[i]->getDensity(position, z);
	return n;
}

std::string DensityList::getDescription() const {
	std::stringstream ss; 
	ss << "DensityList with " << densityList.size() << " modules: \n";
	for (int i = 0; i < densityList.size(); i++) {
		ss << "density " << i + 1 << ": " << densityList[i] -> getDescription();
	}
	
	return ss.str();
}



}  // namespace crpropa

