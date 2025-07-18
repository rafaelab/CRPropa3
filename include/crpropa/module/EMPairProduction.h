#ifndef CRPROPA_EMPAIRPRODUCTION_H
#define CRPROPA_EMPAIRPRODUCTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Sampler.h"


namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */
 
/**
 @class EMPairProduction
 @brief Electron-pair production of photons with background photons.

 This module simulates electron-pair production of cosmic ray photons with background photons:
 gamma + gamma_b --> e+ + e- (Breit-Wheeler process).
 The resulting electron positron pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Optimised sampling strategies (e.g. thinning) are available. 
 Note that if sampling is used, the output must contain the column "weights", which should be included in the post-processing.
 */
class EMPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField; 	// target photon field
	ref_ptr<SamplerEvents> sampler;		// sampler for the interaction rate (thinning)
	bool haveElectrons;					// add secondary electrons to simulation
	double limit;						// limit the step to a fraction of the mean free path
	std::string interactionTag;


	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector< std::vector<double>> tabCDF;  //!< cumulative interaction rate

public:
	/** Constructor
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param sampler		    sampling object (see Sampler.h)
	 @param limit			step size limit as fraction of mean free path
	 */
	EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), double limit = 0.1);

	void setPhotonField(ref_ptr<PhotonField> photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);
	void setSampler(ref_ptr<SamplerEvents> sampler);
	void setInteractionTag(std::string tag);
	
	std::string getInteractionTag() const;

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void performInteraction(Candidate* candidate) const;
	void process(Candidate* candidate) const;
};
/** @}*/


} // namespace crpropa

#endif // CRPROPA_EMPAIRPRODUCTION_H
