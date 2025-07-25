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

[[nodiscard]] double DensityEvolution::getDensity(const Vector3d& position, const double& z) const {
	return density->getDensity(position, z) * std::pow(1 + z, index);
}

[[nodiscard]] std::string DensityEvolution::getDescription() const {
	return std::format("Density evolution with index {} for density: \n{}", index, density->getDescription());
}


///////////////////////////////////////////////////////////////////////////////

DensityGrid::DensityGrid(TargetMedium target, ref_ptr<Grid1f> densityGrid) {
	Density::setTargetMedium(target);
	setGrid(densityGrid);
}

void DensityGrid::setGrid(ref_ptr<Grid1f> g) {
	grid = g;
}

[[nodiscard]] double DensityGrid::getDensity(const Vector3d& position, const double& z) const {
	return grid->interpolate(position);
}

[[nodiscard]] std::string DensityGrid::getDescription() const {
	return std::format("Density field density distribution in a uniform grid\n Target medium: {} (Weight: {})", getTargetMedium().getName(), getTargetMedium().getWeight());
}


///////////////////////////////////////////////////////////////////////////////

DensityList::DensityList() {
}

void DensityList::addDensity(ref_ptr<Density> density) {
    densities.push_back(density);
}

[[nodiscard]] double DensityList::getDensity(const Vector3d& position, const double& z) const {
    double totalDensity = 0;
    for (const auto& density : densities) {
        totalDensity += density->getDensity(position, z);
    }
    return totalDensity;
}

[[nodiscard]] std::string DensityList::getDescription() const {
    std::stringstream ss;
    ss << "DensityList containing " << densities.size() << " densities:\n";
    for (const auto& density : densities) {
        ss << density->getDescription() << "\n";
    }
    return ss.str();
}

}  // namespace crpropa

