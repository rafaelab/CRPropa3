#ifndef CRPROPA_MEDIUMCOMPOSION_H
#define CRPROPA_MEDIUMCOMPOSION_H

#include <cstring>
#include <numeric>
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
		virtual bool isAdmixed() const = 0;
		virtual bool isMolecular() const = 0;
		virtual bool isIonized() const = 0;
		virtual bool isNeutral() const = 0;
		bool isCharged() const {
			return ! isNeutral();
		}
		virtual double getCompositionWeight() const = 0;
		virtual std::string getDescription() const {
			return "MediumComposition (abstract base class)";
		};
		virtual unsigned int numberOfComponents() const {
			return 1; 
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

		bool isAdmixed() const override {
			return false;
		}

		bool isIonized() const override {
			return false;
		}

		bool isNeutral() const override {
			return HepPID::charge(particleId) == 0;
		}

		bool isMolecular() const override {
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

		double getCompositionWeight() const override {
			return 0;
		}

		std::string getDescription() const override {
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

		bool isAdmixed() const override {
			return false;
		}

		bool isIonized() const override {
			return nElectrons != chargeNumber(nucleusId);
		}

		bool isNeutral() const override {
			return chargeNumber(nucleusId) == nElectrons;
		}

		bool isMolecular() const override {
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

		double getCompositionWeight() const override {
			return 1.0 * getNumberOfNucleons();
		}

		std::string getDescription() const override {
			std::stringstream ss;
			ss << "MediumCompositionAtomic with nucleus ID: " << nucleusId << " and number of electrons: " << nElectrons << std::endl;
			return ss.str();
		}
};

/**
 * @class MediumCompositionMolecular
 * @brief Represents a molecule consisting of multiple possibly ionised nuclei.
 * The ionisation state indicates the number of extra/missing electrons compared to the neutral state.
 */
class MediumCompositionMolecular : public MediumComposition {
	protected:
		std::vector<int> nucleiIds;
		int nElectrons;

	public:
		MediumCompositionMolecular(int nElectrons = 0) {
			nucleiIds.clear();
			setNumberOfElectrons(nElectrons);
		}

		MediumCompositionMolecular(const std::vector<int>& ids, int nElectrons) {
			for (const auto& id : ids) {
				add(id);
			}
			setNumberOfElectrons(nElectrons);
		}

		void add(int id) {
			if (! isNucleus(id)) {
				KISS_LOG_WARNING << "MediumCompositionMolecular works for atomic nuclei only. You provided a particle that is not a nucleus." << std::endl;
				throw std::invalid_argument("MediumCompositionMolecular works for atomic nuclei only. You provided a particle that is not a nucleus.");
			}
			nucleiIds.push_back(id);
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

		bool isAdmixed() const override{
			return false;
		}

		bool isIonized() const override {
			return nElectrons != 0;
		}

		bool isNeutral() const override {
			return nElectrons == 0;
		}

		bool isMolecular() const override {
			return true;
		}

		unsigned int getNumberOfNucleons() const {
			int nNucleons = 0;
			for (const auto& nucleusId : nucleiIds)
				nNucleons += massNumber(nucleusId);
			return nNucleons;
		}

		double getCompositionWeight() const override {
			return 1.0 * getNumberOfNucleons();
		}

		std::string getDescription() const override {
			std::stringstream ss;
			ss << "MediumCompositionMolecular with nuclei IDs: ";
			for (const auto& id : nucleiIds) {
				ss << id << " ";
			}
			ss << "and number of electrons: " << nElectrons << std::endl;
			return ss.str();
		}
};


/**
 * @class MediumCompositionList
 * @brief Represents a list of medium compositions, allowing for complex admixtures.
 * This class can be used to represent a mixture of different medium compositions.
 * The weights indicate the relative abundance of each composition in the mixture (not weighted by the number of nucleons).
 */
class MediumCompositionList : public MediumComposition {
	protected:
		std::vector<ref_ptr<MediumComposition>> compositions;
		std::vector<double> weights;

	public:
		MediumCompositionList() {
		}

		MediumCompositionList(const std::vector<ref_ptr<MediumComposition>>& comps, const std::vector<double>& weights) {
			if (comps.size() != weights.size()) {
				KISS_LOG_WARNING << "MediumCompositionList: Compositions and weights must have the same size." << std::endl;
				throw std::invalid_argument("MediumCompositionList: Compositions and weights must have the same size.");
			}
			
			for (size_t i = 0; i < comps.size(); ++i)
				add(comps[i], weights[i]);
		}

		void add(ref_ptr<MediumComposition> composition, double weight = 1.0) {
			compositions.push_back(composition);
			weights.push_back(weight);
		}

		const std::vector<ref_ptr<MediumComposition>>& getCompositions() const {
			return compositions;
		}

		const std::vector<double>& getWeights() const {
			return weights;
		}

		bool isAdmixed() const override {
			return compositions.size() > 1;
		}

		bool isMolecular() const override {
			return false; 
		}

		bool isIonized() const override {
			return false;
		}

		bool isNeutral() const override {
			return false;
		}

		unsigned int numberOfComponents() const override {
			return compositions.size();
		}

		double getCompositionWeight() const override {
			double weight = 0.0;
			double sumWeights = accumulate(weights.begin(), weights.end(), 0.);
			for (size_t i = 0; i < numberOfComponents(); ++i) {
				weight += (weights[i] * compositions[i]->getCompositionWeight() / sumWeights);
			}

			return weight;
		}

		std::string getDescription() const override {
			std::stringstream ss;
			ss << "MediumCompositionList with " << compositions.size() << " components." << std::endl;
			for (const auto& comp : compositions)
				ss << comp->getDescription();

			return ss.str();
		}
};


}  // namespace crpropa

#endif  // CRPROPA_MEDIUMCOMPOSITION_H
