#ifndef CRPROPA_EMPAIRPRODUCTION_H
#define CRPROPA_EMPAIRPRODUCTION_H

#include <fstream>
#include <cmath>

#include "crpropa/Common.h"
#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"


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
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class EMPairProduction: public Module {
private:
	ref_ptr<PhotonField> photonField; 	// target photon field
	bool haveElectrons;					// add secondary electrons to simulation
	double limit;						// limit the step to a fraction of the mean free path
	bool forwardApproximation;				// whether to correct the angle of the pairs (defaults to false)
	double thinning;					// factor of the thinning (0: no thinning, 1: maximum thinning)
	std::string interactionTag = "EMPP";

	// tabulated interaction rate 1/lambda(E)
	std::vector<double> tabEnergy;  //!< electron energy in [J]
	std::vector<double> tabRate;  //!< interaction rate in [1/m]
	
	// tabulated CDF(s_kin, E) = cumulative differential interaction rate
	std::vector<double> tabE;  //!< electron energy in [J]
	std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
	std::vector<std::vector<double>> tabCDF;  //!< cumulative interaction rate

public:
	/** Constructor
	 @param photonField		  		target photon field
	 @param haveElectrons	  		if true, add secondary electrons as candidates
	 @param thinning		  		weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
	 @param limit			  		step size limit as fraction of mean free path
	 @param forwardApproximation  	whether to assume that pairs are produced in the same direction as parent (defaults to true)
	 */
	EMPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, double thinning = 0, double limit = 0.1, bool forwardApproximation = true);

	// set the target photon field
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary electrons are added to the simulation
	void setHaveElectrons(bool haveElectrons);

	/** Limit the propagation step to a fraction of the mean free path
	 * @param limit fraction of the mean free path
	 */
	void setLimit(double limit);

	/** Apply thinning with a given thinning factor
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setThinning(double thinning);

	/** Determines whether the angle between the electron/positron will be considered or if it will be colinear with parent.
	 * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
	 */
	void setForwardApproximation(bool b);

	

	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */	
	void setInteractionTag(std::string tag);
	std::string getInteractionTag() const;

	void initRate(std::string filename);
	void initCumulativeRate(std::string filename);

	void performInteraction(Candidate *candidate) const;
	void process(Candidate *candidate) const;
};
/** @}*/


} // namespace crpropa

#endif // CRPROPA_EMPAIRPRODUCTION_H
