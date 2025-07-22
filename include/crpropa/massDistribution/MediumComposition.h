#ifndef CRPROPA_MEDIUMCOMPOSION_H
#define CRPROPA_MEDIUMCOMPOSION_H


#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "crpropa/ParticleID.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"
#include "HepPID/ParticleIDMethods.hh"
#include "kiss/logger.h"



namespace crpropa {

/**
 * @class MediumComposition
 * @brief Abstract base class for medium compositions.
 * 
 * This class serves as an interface for different types of medium compositions, such as elementary particles, atomic nuclei, and molecular compositions.
 */
class MediumComposition: public Referenced {
	public:
		virtual ~MediumComposition() = default;
		virtual bool isAdmixed() const = 0;
		virtual bool isMolecular() const = 0;
		virtual bool isIonized() const = 0;
		virtual bool isNeutral() const = 0;
		bool isCharged() const {
			return ! isNeutral();
		};
		virtual double getCompositionWeight() const = 0;
		virtual std::string getDescription() const = 0;
		virtual unsigned int numberOfComponents() const {
			return 1;
		};
};


/**
 * @class MediumCompositionElementary
 * @brief A class representing an elementary medium composition.
 * 
 * This class allows for the representation of fundamental particles, such as electrons, with their particle ID.
 */
class MediumCompositionElementary : public MediumComposition {
	protected:
		int particleId;

	public:
		MediumCompositionElementary();
		MediumCompositionElementary(int id);
		void setParticleId(int id);
		int getParticleId() const;
		bool isAdmixed() const override;
		bool isIonized() const override;
		bool isNeutral() const override;
		bool isMolecular() const override;
		unsigned int getNumberOfNucleons() const;
		double getCompositionWeight() const override;
		std::string getDescription() const override;
};


/**
 * @class MediumCompositionAtomic
 * @brief A class representing an atomic medium composition.
 * 
 * This class allows for the representation of atomic nuclei, including their number of electrons.
 */
class MediumCompositionAtomic : public MediumComposition {
	protected:
		int nucleusId;
		int nElectrons;

	public:
		MediumCompositionAtomic();
		MediumCompositionAtomic(int id, int nElectrons = 0);
		void setNucleusId(int id);
		void setNumberOfElectrons(int n);
		int getNucleusId() const;
		int getNumberOfElectrons() const;
		bool isAdmixed() const override;
		bool isIonized() const override;
		bool isNeutral() const override;
		bool isMolecular() const override;
		unsigned int getNumberOfNucleons() const;
		double getCompositionWeight() const override;
		std::string getDescription() const override;
};



/**
 * @class MediumCompositionMolecular
 * @brief A class representing a molecular medium composition.
 * 
 * This class allows for the representation of molecular compositions, such as H2, by storing a list of nuclei IDs and the number of electrons.
 */
class MediumCompositionMolecular : public MediumComposition {
	protected:
		std::vector<int> nucleiIds;
		int nElectrons;

	public:
		MediumCompositionMolecular();
		MediumCompositionMolecular(int nElectrons);
		MediumCompositionMolecular(const std::vector<int>& ids, int nElectrons = 0);
		void add(int id);
		void setNucleiIds(const std::vector<int>& ids);
		void setNumberOfElectrons(int n);
		int getNumberOfElectrons() const;
		const std::vector<int>& getNucleiIds() const;
		bool isAdmixed() const override;
		bool isIonized() const override;
		bool isNeutral() const override;
		bool isMolecular() const override;
		unsigned int getNumberOfNucleons() const;
		double getCompositionWeight() const override;
		std::string getDescription() const override;
};


/**
 * @class MediumCompositionList
 * @brief A list of MediumComposition objects with associated weights.
 */
class MediumCompositionList : public MediumComposition {
	protected:
		std::vector<ref_ptr<MediumComposition>> compositions;
		std::vector<double> weights;

	public:
		MediumCompositionList();
		MediumCompositionList(const std::vector<ref_ptr<MediumComposition>>& comps, const std::vector<double>& weights = {});
		void add(ref_ptr<MediumComposition> composition, double weight = 1.0);
		const std::vector<ref_ptr<MediumComposition>>& getCompositions() const;
		const std::vector<double>& getWeights() const;
		bool isAdmixed() const override;
		bool isMolecular() const override;
		bool isIonized() const override;
		bool isNeutral() const override;
		unsigned int numberOfComponents() const override;
		double getCompositionWeight() const override;
		std::string getDescription() const override;
};

} // namespace crpropa

#endif // CRPROPA_MEDIUMCOMPOSITION_H
