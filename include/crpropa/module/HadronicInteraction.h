#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H


#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "crpropa/Module.h"
#include "crpropa/Sampler.h"
#include "crpropa/Source.h"
#include "crpropa/Units.h"
#include "crpropa/massDistribution/Density.h"
#include "kiss/logger.h"
#include "HepPID/ParticleIDMethods.hh"



namespace crpropa {


/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 * @class HadronicInteraction
 * @brief Interaction between nuclei.
 * 
 * This module simulates hadronic interactions of cosmic-ray nuclei with a target density field.
 * The interaction is sampled from a tabulated cross section file.
 * This closely follows the implementation from 
 * 	https://gitlab.ruhr-uni-bochum.de/doernjkj/hadronic-interaction-in-crpropa/ ,
 * by J. Dörner, with modifications by C. Prévotat to extend the energy range and allow heavier targets.
 */
class HadronicInteraction : public Module {
	public:
		class CrossSection {
			public:
				inline static const std::map<int, double> secondary_mc2 = {
					{-11, mass_electron * c_squared / GeV}, 
					{11, mass_electron * c_squared / GeV}, 
					{22, 0}, 
					{14, 0}, 
					{-14, 0}, 
					{12, 0}, 
					{-12, 0}, 
					{1000010010, mass_proton * c_squared / GeV}, 
					{-1000010010, mass_proton * c_squared / GeV}, 
					{1000000010, mass_neutron * c_squared / GeV}, 
					{-1000000010, mass_neutron * c_squared / GeV}
				}; // should be in GeV
				inline static const std::map<int, const std::string> restEnergiesNuclei = {
					{-1000000010, "an"},
					{-1000010010, "ap"},
					{1000020040, "He"},
					{1000060120, "C"},
					{1000130260, "Al"},
					{1000260560, "Fe"}
				}; // in J

				std::string filename;
				std::vector<double> energyPrimary; //< tabulated values of the kinetic energy of the primary
				std::vector<double> energySecondary; //< tabulated values of the energy of the secondary particle
				std::vector<std::vector<double>> sigma; //< integrated differential cross section (CDF); sigma[Tp][eps]
				double weightEnergyLoss; //< weight of the energy loss 
				// bool addSecondary; //< decide if secondaries are injected to the simulation
				int primaryId; // to know which table should be loaded ; need to know the CR and the target (p or He)
				int secondaryId; //< id of the secondary particle
				// int targetId; //< id of the target particles/medium
				ref_ptr<MediumComposition> targetMedium; //< target medium composition for the cross section

			public: 
				CrossSection();
				CrossSection(int primaryId, int secondaryId, ref_ptr<MediumComposition> target, double weightEnergyLoss = 1.0);
				void loadData();
				void setPrimaryId(int id);
				void setSecondaryId(int id);
				void setTargetMedium(ref_ptr<MediumComposition> target);
				void setWeightEnergyLoss(double weight);
				int getPrimaryId() const;
				int getSecondaryId() const;
				double getWeightEnergyLoss() const;
				ref_ptr<MediumComposition> getTargetMedium() const;

				static double totalInelasticCrossSectionProton(const double& energy);
				static double totalInelasticCrossSection(const double& T, const int& id, const int& At);
		};

	protected:
		void loadData(std::string file, CrossSection& cs);

	private:
		ref_ptr<Density> density; //< target field (proton number density)
		double limit; //< limit the next step to a fraction of the mean free path
		bool catastrophicLoss; //< deactivate the primary particle after the interaction
		bool havePhotons; //< if true, photons are produced and added as candidates
		bool haveElectrons; //< if true, electrons are produced and added as candidates
		bool haveNeutrinos; //< if true, neutrinos are produced and added as candidates
		bool haveNucleons; //< if true, nucleons are produced and added as candidates
		bool haveAntiNucleons; //< if true, anti-nucleons are produced 
		std::string interactionTag;
		std::vector<int> listOfPrimaries; //< list of primary IDs for which cross sections are loaded
		std::map<int, std::vector<CrossSection>> crossSectionList; // list of the included crosssections; dict : id of the CR
		ref_ptr<SamplerEvents> sampler;		// sampler for the interaction rate (thinning)

		inline static const std::map<int, double> restEnergiesNuclei = {
			{1000010010, mass_proton * c_squared }, 
			{1000020040, 4.002602 * c_squared * amu}, 
			{1000060120, 12.0106 * c_squared * amu}, 
			{1000130260, 26.98154 * c_squared * amu}, 
			{1000260560, 55.845 * c_squared * amu} 
			}; // in J

		inline static const std::map<int, double> minEnergiesNuclei = {
			{1000010010, 5 * GeV }, 
			{1000020040, 5 * GeV}, 
			{1000060120, 16 * GeV}, 
			{1000130260, 40 * GeV}, 
			{1000260560, 100 * GeV} 
			}; // in J

		inline static const std::vector<int> availableNuclei = {
			1000010010, // hydrogen
			1000020040, // helium
			1000060120, // carbon
			1000130260, // aluminium
			1000260560  // iron
		}; // list of available nuclei in the cross section files

	public:
		HadronicInteraction(ref_ptr<Density> density, bool havePhotons = false, bool haveElectrons = false, bool haveNeutrinos = false, bool haveNucleons = false, bool haveAntiNucleons = false, bool catastrophic = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), double limit = 0.1);
		void initData();
		void setHavePhotons(bool b);
		void setHaveElectrons(bool b);
		void setHaveNeutrinos(bool b);
		void setHaveNucleons(bool b);
		void setHaveAntiNucleons(bool b);
		void setLimit(double limit); 
		void setDensity(ref_ptr<Density> density);
		void setInteractionTag(std::string tag);
		void setCatastrophic(bool b);
		void setSampler(ref_ptr<SamplerEvents> sampler);
		void setListOfPrimaries(std::vector<int> primaries);
		bool getHavePhotons() const;
		bool getHaveElectrons() const;
		bool getHaveNeutrinos() const;
		bool getHaveNucleons() const;
		bool getHaveAntiNucleons() const;
		bool isCatastrophic() const;
		double getLimit() const;
		std::string getInteractionTag() const;
		void addCrossSection(CrossSection cs);
		bool isPrimaryImplemented(const int& id) const;
		void process(Candidate* candidate) const;
		void performInteraction(Candidate* candidate, const int& At) const;
		double getDensityAtPosition(const Vector3d& pos) const;

		// void printCrosssections() {
		// 	for(int i = 0; i < crosssectionList[1000010010].size(); i++) { // [1000010010] : print CS of protons
		// 		std::cout << i << "\t" << crosssectionList[1000010010][i].secondaryId << "\n";

		// 		std::vector< std::vector<double> > sigma = crosssectionList[1000010010][i].sigma;
		// 		for(int j = 0; j < 5; j++) {
		// 			std::cout << "\t";
		// 			for(int k = 0; k < 5; k++) 
		// 				std::cout << sigma[j][k] << "\t";
		// 			std::cout << "\n";
		// 		}
		// 	}
		// }
};

}



#endif // CRPROPA_HADRONICINTERACTION_H
