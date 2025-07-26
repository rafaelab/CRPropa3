#ifndef CRPROPA_TARGETMEDIUM_H
#define CRPROPA_TARGETMEDIUM_H

#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility> 
#include <vector>

#include "crpropa/ParticleID.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"
#include "HepPID/ParticleIDMethods.hh"
#include "kiss/logger.h"


namespace crpropa {


/**
 * @brief Class defining a given target medium
 *
 * Provides virtual methods for name and weight, allowing runtime polymorphism.
 * It has a given name (e.g., HI, HII, H2, etc), and a tag (e.g., H), used for reading data files.
 * The weight is used to scale the density of the target medium.
 */
class TargetMedium {
	protected:
		std::string name = "unknown target"; 
		std::string tag = "unknown";
		double weight = 1.0; 

	public:
		explicit TargetMedium(double weight = 1) noexcept;
		TargetMedium(std::string_view name, double weight = 1) noexcept;
		TargetMedium(std::string_view name, std::string_view tag, double weight = 1) noexcept; 
		~TargetMedium() = default;
		void setName(std::string n) noexcept;
		void setTag(std::string t) noexcept;
		void setWeight(double w) noexcept;
		[[nodiscard]] std::string getName() const noexcept;
		[[nodiscard]] std::string getTag() const noexcept;
		[[nodiscard]] double getWeight() const noexcept;
		[[nodiscard]] std::string getDescription() const noexcept;
};


using TargetMediumPtrS = std::shared_ptr<TargetMedium>;
using TargetMediumPtrU = std::unique_ptr<TargetMedium>;


} // namespace crpropa

#endif // CRPROPA_TARGETMEDIUM_H
