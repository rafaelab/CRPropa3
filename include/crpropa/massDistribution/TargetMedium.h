#ifndef CRPROPA_TARGETMEDIUM_H
#define CRPROPA_TARGETMEDIUM_H


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
 */
class TargetMedium {
	protected:
		std::string name = "unknown target"; ///< Default member initializer
		double weight = 1.0; ///< Default member initializer

	public:
		explicit TargetMedium(double weight = 1) noexcept;
		TargetMedium(std::string_view name, double weight = 1) noexcept; 
		void setName(std::string n) noexcept;
		void setWeight(double w) noexcept;
		[[nodiscard]] std::string getName() const noexcept;
		[[nodiscard]] double getWeight() const noexcept;
		[[nodiscard]] std::string getDescription() const noexcept;
};


static inline const TargetMedium TargetHI = TargetMedium("HI", 1.0);
static inline const TargetMedium TargetHII = TargetMedium("HII", 1.0);
static inline const TargetMedium TargetH2 = TargetMedium("H2", 2.0);
static inline const TargetMedium TargetHe = TargetMedium("He", 4.0);


} // namespace crpropa

#endif // CRPROPA_TARGETMEDIUM_H
