#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H


#include <string>
#include <fstream>
#include <limits>

#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "crpropa/Module.h"
#include "crpropa/Sampler.h"
#include "crpropa/Source.h"
#include "crpropa/Units.h"
#include "crpropa/massDistribution/Density.h"
#include "kiss/logger.h"



namespace crpropa {


/**
 * \addtogroup EnergyLosses
 * @{
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

				int nucleusId; //< id of the secondary particle
				std::vector<double> Tp; //< tabulated values of the kinetic energy of the primary
				std::vector<double> eps; //< tabulated values of the energy of the secondary particle
				std::vector< std::vector<double>> sigma; //< integrated differential crossecton (CDF); sigma[Tp][eps]
				double weightEnergyLoss; //< weight of the energy loss 
				bool addSecondary; //< decide if secondaries are injected to the simulation
				int CR_Id; // to know which table should be loaded ; need to know the CR and the target (p or He)

			public: 
				CrossSection();
				CrossSection(std::string file);
				CrossSection(std::string file, int id, double lossFactor, bool addSecondary, int CR_Id);
				void loadData(std::string file);
		};

	protected:
		void loadData(std::string file, CrossSection& cs);

	private:
		ref_ptr<Density> density; //< target field (proton number density)
		double limit; //< limit the next step to a fraction of the mean free path
		bool deactivatePrimary; //< deactivate the primary particle after the interaction
		std::string interactionTag;
		std::map<int, std::vector<CrossSection>> crossSectionList; // list of the included crosssections; dict : id of the CR


		inline static const std::map<int, double> CR_rest_energy = {
			{1000010010, mass_proton * c_squared }, 
			{1000020040, 4.002602 * c_squared * amu}, 
			{1000060120, 12.0106 * c_squared * amu}, 
			{1000130260, 26.98154 * c_squared * amu}, 
			{1000260560, 55.845 * c_squared * amu} 
			}; // in J

		inline static const std::map<int, double> CR_min_energy = {
			{1000010010, 5 * GeV }, 
			{1000020040, 5 * GeV}, 
			{1000060120, 16 * GeV}, 
			{1000130260, 40 * GeV}, 
			{1000260560, 100 * GeV} 
			}; // in J

	public:
		HadronicInteraction(ref_ptr<Density> density);
		HadronicInteraction(ref_ptr<Density> density, std::string configFile, double limit = 0.1);


		void initData(std::string configFile);	
		void addCrossSection(CrossSection cs);

		void process(Candidate* candidate) const;
		void performInteraction(Candidate* candidate, int A_t) const;

		/* Total inelastic crosssection following the parametrisation from 
			@param Tp:	kinetic energy of the proton in [J]
		*/
		double totalInelasticCrossSection(double Tp, int id, int A_t) const;

		// /* allow the given nucleusId as a secondary
		// 	Search in the list of crosssections for the nucleus id and overwrite the information from the config file.
		// 	@param nucleusId:	nucleusId of the secondary
		// 	@param allow:	if true secondaries are allowed
		// */
		void allowSecondaryId(int nucleusId, bool allow = true);

		double getDensityAtPosition(Vector3d& pos) const;

		void setLimit(double limit); 
		double getLimit() const;

		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;

		void setDeactivatePrimary(bool b);
		bool getDeactivatePrimary() const;

		// void printCrosssections() {
		// 	for(int i = 0; i < crosssectionList[1000010010].size(); i++) { // [1000010010] : print CS of protons
		// 		std::cout << i << "\t" << crosssectionList[1000010010][i].nucleusId << "\n";

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
