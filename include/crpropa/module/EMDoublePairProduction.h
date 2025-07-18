#ifndef CRPROPA_EMDOUBLEPAIRPRODUCTION_H
#define CRPROPA_EMDOUBLEPAIRPRODUCTION_H

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Random.h"
#include "crpropa/Sampler.h"
#include "crpropa/Units.h"


namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class EMDoublePairProduction
 @brief Electron double pair production of photons with background photons.

 This module simulates electron double pair production of photons with background photons for several photon fields.
 The secondary electrons from this interaction are optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Optimised sampling strategies (e.g. thinning) are available. 
 Note that if sampling is used, the output must contain the column "weights", which should be included in the post-processing.
 */
class EMDoublePairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField;
	ref_ptr<SamplerEvents> sampler;
	bool haveElectrons;
	double limit;
	std::string interactionTag;

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]

public:
	/** Constructor
	 @param photonField		target photon field
	 @param haveElectrons	if true, add secondary electrons as candidates
	 @param sampler		    sampling object (see Sampler.h)
	 @param limit			step size limit as fraction of mean free path
	 */
	EMDoublePairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), double limit = 0.1);

	void setPhotonField(ref_ptr<PhotonField> photonField);
	void setHaveElectrons(bool haveElectrons);
	void setLimit(double limit);
	void setSampler(ref_ptr<SamplerEvents> sampler);
	void setInteractionTag(std::string tag);

	std::string getInteractionTag() const;

	void initRate(std::string filename);
	void process(Candidate* candidate) const;
	void performInteraction(Candidate* candidate) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_EMDOUBLEPAIRPRODUCTION_H
