#ifndef CRPROPA_MEDIUMCOMPOSION_H
#define CRPROPA_MEDIUMCOMPOSION_H

#include <cstring>
#include <sstream>

#include "crpropa/ParticleID.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"
#include "HepPID/ParticleIDMethods.hh"
#include "kiss/logger.h"

namespace crpropa {


class MediumComposition: public Referenced {
	public:
		virtual ~MediumComposition() = default;
		virtual bool isMolecular() const = 0;
		virtual bool isIonized() const = 0;
		virtual bool isNeutral() const = 0;
		bool isCharged() const {
			return ! isNeutral();
		}
		virtual unsigned int getNumberOfNucleons() const = 0;
		virtual std::string getDescription() const {
			return "MediumComposition (abstract base class)";
		};
}; 


/**
 * @class MediumCompositionElementary
 * @brief Represents an ionized nucleus with a specific number of electrons knocked off.
 * The ionisation state indicates the level of ionisation: 0 if neutral, > 0 if ionised. 
 */
class MediumCompositionElementary : public MediumComposition {
	protected:
		int particleId; 
	
	public:
		MediumCompositionElementary(int id) {
			setParticleId(id);
		}

		void setParticleId(int id) {
			if (isNucleus(id)) {
				KISS_LOG_WARNING << "MediumCompositionElementary is more appropriate for more fundamental particles. The Id provided is for a nucleus. Try using MediumCompositionAtomic. I will proceed anyway." << std::endl;
			}
			particleId = id;
		}

		int getParticleId() const {
			return particleId;
		}

		bool isIonized() const {
			return false;
		}

		bool isNeutral() const { 
			return HepPID::charge(particleId) == 0;
		}

		bool isMolecular() const {
			return false;
		}

		unsigned int getNumberOfNucleons() const {
			if (particleId == 2212 or particleId == 2112) { // p or n
				return 1;
			} 
			if (isNucleus(particleId)) {
				return massNumber(particleId);
			}
		
			return 0; // no nucleons
		}

		std::string getDescription() const {
			std::stringstream ss;
			ss << "MediumCompositionElementary with particle ID: " << particleId << std::endl;
			return ss.str();
		}
};


/**
 * @class MediumCompositionAtomic
 * @brief Represents an ionized nucleus with a specific number of electrons in the atomic shell.
 */
class MediumCompositionAtomic : public MediumComposition {
	protected:
		int nucleusId;
		int nElectrons;  
	
	public:
		MediumCompositionAtomic(int id, int nElectrons) {
			setNucleusId(id);
			setNumberOfElectrons(nElectrons);
		}

		void setNucleusId(int id) {
			if (! isNucleus(id)) {
				KISS_LOG_WARNING << "MediumCompositionAtomic works for atomic nuclei only. You provided a particle that is not a nucleus." << std::endl;
				throw std::invalid_argument("MediumCompositionAtomic works for atomic nuclei only. You provided a particle that is not a nucleus.");
			} 
			nucleusId = id;
		}

		void setNumberOfElectrons(int n) {
			nElectrons = n;
		}

		int getNucleusId() const {
			return nucleusId;
		}

		int getNumberOfElectrons() const {
			return nElectrons;
		}

		bool isIonized() const {
			return nElectrons != chargeNumber(nucleusId);
		}

		bool isNeutral() const {
			return chargeNumber(nucleusId) == nElectrons;
		}

		bool isMolecular() const {
			return false;
		}

		unsigned int getNumberOfNucleons() const {
			if (nucleusId == 2212 or nucleusId == 2112) { // p or n
				return 1;
			} 
			if (isNucleus(nucleusId)) {
				return massNumber(nucleusId);
			}
		
			return 0; // no nucleons
		}

		std::string getDescription() const {
			std::stringstream ss;
			ss << "MediumCompositionAtomic with nucleus ID: " << nucleusId << " and number of electrons: " << nElectrons << std::endl;
			return ss.str();
		}
};

/**
 * @class MediumCompositionMolecular
 * @brief Represents a molecule consisting of multiple possibly ionized nuclei.
 * The ionisation state indicates the number of extra/missing electrons compared to the neutral state.
 */
class MediumCompositionMolecular : public MediumComposition {
	protected:
		std::vector<int> nucleiIds;
		int nElectrons;

	public:
		MediumCompositionMolecular(const std::vector<int>& ids, int nElectrons) {
			setNucleiIds(ids);
			setNumberOfElectrons(nElectrons);
		}

		void setNucleiIds(const std::vector<int>& ids) {
			nucleiIds = ids;
		}

		void setNumberOfElectrons(int n) {
			nElectrons = n;
		}

		int getNumberOfElectrons() const {
			return nElectrons;
		}

		const std::vector<int>& getNucleiIds() const {
			return nucleiIds;
		}

		bool isIonized() const {
			return nElectrons != 0;
		}

		bool isNeutral() const {
			return nElectrons == 0;
		}

		bool isMolecular() const {
			return true;
		}

		std::string getDescription() const {
			std::stringstream ss;
			ss << "MediumCompositionMolecular with nuclei IDs: ";
			for (const auto& id : nucleiIds) {
				ss << id << " ";
			}
			ss << "and number of electrons: " << nElectrons << std::endl;
			return ss.str();
		}
};


}  // namespace crpropa

#endif  // CRPROPA_MediumComposition_H
