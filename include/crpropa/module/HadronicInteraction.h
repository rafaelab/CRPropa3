#pragma once

#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "crpropa/Module.h"
#include "crpropa/Sampler.h"
#include "crpropa/Source.h"
#include "crpropa/Units.h"
#include "crpropa/massDistribution/Density.h"
#include "crpropa/massDistribution/TargetMedium.h"
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
					{1000130270, "Al"},
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
				TargetMedium targetMedium; //< target medium composition for the cross section

			public: 
				CrossSection();
				CrossSection(int primaryId, int secondaryId, TargetMedium target, double weightEnergyLoss = 1.0);
				void loadData();
				void loadDataFromFile(const std::string& filename);
				void setPrimaryId(int id);
				void setSecondaryId(int id);
				void setTargetMedium(TargetMedium target);
				void setWeightEnergyLoss(double weight);
				int getPrimaryId() const;
				int getSecondaryId() const;
				double getWeightEnergyLoss() const;
				TargetMedium getTargetMedium() const;

				static double totalInelasticCrossSectionProton(const double& energy);
				static double totalInelasticCrossSection(const double& T, const int& id, const int& At);
		};

	private:
		std::vector<std::pair<ref_ptr<Density>, TargetMedium>> media; //< list of media with their target medium
		double limit; //< limit the next step to a fraction of the mean free path
		bool catastrophicLoss; //< deactivate the primary particle after the interaction
		bool havePhotons; //< if true, photons are produced and added as candidates
		bool haveElectrons; //< if true, electrons are produced and added as candidates
		bool haveNeutrinos; //< if true, neutrinos are produced and added as candidates
		bool haveNucleons; //< if true, nucleons are produced and added as candidates
		bool haveAntiNucleons; //< if true, anti-nucleons are produced 
		std::string interactionTag;
		std::vector<int> listOfPrimaries; //< list of primary IDs for which cross sections are loaded
		std::map<int, std::vector<CrossSection>> crossSectionList; // list of the included cross sections; keys: primary IDs, values: vector of cross sections for this primary
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
		HadronicInteraction(bool havePhotons = false, bool haveElectrons = false, bool haveNeutrinos = false, bool haveNucleons = false, bool haveAntiNucleons = false, bool catastrophic = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), double limit = 0.1);
		HadronicInteraction(ref_ptr<Density> density, TargetMedium target, bool havePhotons = false, bool haveElectrons = false, bool haveNeutrinos = false, bool haveNucleons = false, bool haveAntiNucleons = false, bool catastrophic = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), double limit = 0.1);
		HadronicInteraction(std::vector<ref_ptr<Density>> densities, std::vector<TargetMedium> targets, bool havePhotons = false, bool haveElectrons = false, bool haveNeutrinos = false, bool haveNucleons = false, bool haveAntiNucleons = false, bool catastrophic = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), double limit = 0.1);
		void initData();
		void setHavePhotons(bool b) noexcept;
		void setHaveElectrons(bool b) noexcept;
		void setHaveNeutrinos(bool b) noexcept;
		void setHaveNucleons(bool b) noexcept;
		void setHaveAntiNucleons(bool b) noexcept;
		void setCatastrophic(bool b) noexcept;
		void setLimit(double limit) noexcept; 
		void setInteractionTag(std::string tag);
		void setSampler(ref_ptr<SamplerEvents> sampler);
		void setListOfPrimaries(std::span<const int> primaries);
		bool getHavePhotons() const;
		bool getHaveElectrons() const;
		bool getHaveNeutrinos() const;
		bool getHaveNucleons() const;
		bool getHaveAntiNucleons() const;
		bool isCatastrophic() const;
		double getLimit() const;
		std::string getInteractionTag() const;
		std::vector<int> getListOfPrimaries() const;
		void addMedium(ref_ptr<Density> density, TargetMedium target);
		void addCrossSection(CrossSection cs);
		bool isPrimaryImplemented(const int& id) const;
		double getDensity(const Vector3d& pos, const double& z) const;
		void process(Candidate* candidate) const;
		void performInteraction(Candidate* candidate, const int& At) const;

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

/** @} 
 */

} // namespace crpropa



