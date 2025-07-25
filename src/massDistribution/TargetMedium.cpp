#include "crpropa/massDistribution/TargetMedium.h"


namespace crpropa {

TargetMedium::TargetMedium(double weight) noexcept 
	: weight(weight) {
}

TargetMedium::TargetMedium(std::string_view name, double weight) noexcept 
	: name(name), weight(weight) {
}

void TargetMedium::setName(std::string n) noexcept {
	name = std::move(n);
}

void TargetMedium::setWeight(double w) noexcept {
	weight = w;
}

[[nodiscard]] std::string TargetMedium::getName() const noexcept {
	return name;
}

[[nodiscard]] double TargetMedium::getWeight() const noexcept {
	return weight;
}

[[nodiscard]] std::string TargetMedium::getDescription() const noexcept {
	std::ostringstream oss;
	oss << "Target " << name << " (w = " << weight << ")";
	return oss.str();
}

} // namespace crpropa

