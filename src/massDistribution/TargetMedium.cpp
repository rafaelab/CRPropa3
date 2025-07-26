#include "crpropa/massDistribution/TargetMedium.h"


namespace crpropa {

TargetMedium::TargetMedium(double weight) noexcept 
	: weight(weight) {
}

TargetMedium::TargetMedium(std::string name, double weight) noexcept 
	: name(name), tag(name), weight(weight) {
}

TargetMedium::TargetMedium(std::string name, std::string tag, double weight) noexcept 
	: name(name), tag(tag), weight(weight) {
}

TargetMedium::~TargetMedium() {
}

void TargetMedium::setName(std::string n) noexcept {
	name = std::move(n);
}

void TargetMedium::setTag(std::string t) noexcept {
	tag = std::move(t);
}

void TargetMedium::setWeight(double w) noexcept {
	weight = w;
}

[[nodiscard]] std::string TargetMedium::getName() const noexcept {
	return name;
}

[[nodiscard]] std::string TargetMedium::getTag() const noexcept {
	return tag;
}

[[nodiscard]] double TargetMedium::getWeight() const noexcept {
	return weight;
}

[[nodiscard]] std::string TargetMedium::getDescription() const noexcept {
	std::ostringstream oss;
	oss << "Target " << name << "; tag: " << tag << " (w = " << weight << ")";
	return oss.str();
}

} // namespace crpropa

