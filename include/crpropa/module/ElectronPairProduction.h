#ifndef CRPROPA_ELECTRONPAIRPRODUCTION_H
#define CRPROPA_ELECTRONPAIRPRODUCTION_H

#include <cstddef>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <typeinfo>

#include "crpropa/Module.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
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
 * @class ElectronPairProduction
 * @brief Bethe-Heitler pair production of charged nuclei with background photons.
 *
 * This module simulates electron-pair production as a continuous energy loss.
 * Several photon fields can be selected.
 * The production of secondary e+/e- pairs and photons can by activated.
 * By default, the module limits the step size to 10% of the energy loss length of the particle.
 *
 * The secondary electrons and positrons may be added to the simulation.
 * Moreover, optimisations on the sampling strategy of secondaries (acceptance/rejection => SamplerEvents; draw from distribution => SamplerDistribution) are available.
 */
class ElectronPairProduction: public Module {
	protected:
		ref_ptr<PhotonField> photonField;
		std::vector<double> tabLossRate; /*< tabulated energy loss rate in [J/m] for protons at z = 0 */
		std::vector<double> tabLorentzFactor; /*< tabulated Lorentz factor */
		std::vector<std::vector<double>> tabSpectrum; /*< electron/positron cdf(Ee|log10(gamma)) for log10(Ee/eV)=7-24 in 170 steps and log10(gamma)=6-13 in 70 steps and*/
		double limit; ///< fraction of energy loss length to limit the next step
		bool haveElectrons; /*< if true, secondary electrons will be added to the simulation */
		std::string interactionTag;
		ref_ptr<SamplerEvents> samplerEvents;
		ref_ptr<SamplerDistribution> samplerDistribution;
		int maximumSamples; ///< maximum number of samples of synchrotron photons (break condition; defaults to 00; 0 or <0 means no sampling)

	public:
		/**
		 * @brief Constructor for the ElectronPairProduction
		 * 
		 * @param photonField 	target photon field
		 * @param haveElectrons if true, secondary electrons will be added to the simulation
		 * @param samplingEvts 	sampler for events (acceptance/rejection of individual secondaries)
		 * @param samplingDist 	sampler for randomly drawing a few secondaries from the distribution
		 * @param limit 		step size limit as fraction of mean free path
		 */
		ElectronPairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons = false, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()), ref_ptr<SamplerDistribution> samplingDist = nullptr, int maximumSamples = 0, double limit = 0.1);

		// set a custom interaction tag to trace back this interaction
		void setInteractionTag(std::string tag);

		// set the target photon field
		void setPhotonField(ref_ptr<PhotonField> photonField);

		// decide if secondary electrons are added to the simulation
		void setHaveElectrons(bool haveElectrons);
		
		// limit the propagation step to a fraction of the mean free path
		void setLimit(double limit);

		// set the sampler for events (acceptance/rejection of individual secondaries)
		void setSamplerEvents(ref_ptr<SamplerEvents> sampler);

		// set the sampler for randomly drawing a few secondaries from the distribution (stored in a histogram)
		void setSamplerDistribution(ref_ptr<SamplerDistribution>sampler);

		// the distribution sampler will draw up to `nSamples` secondaries
		void setMaximumSamples(int nSamples);
		
		std::string getInteractionTag() const;
		
		void initRate(std::string filename);
		void initSpectrum(std::string filename);

		void process(Candidate* candidate) const;

		/**
		 * Calculates the energy loss length 1/beta = -E dx/dE in [m]
		 * @param	id		PDG particle ID
		 * @param lf		Lorentz factor
		 * @param z		redshift
		 *
		 * The energy loss length is tabulated for protons against CMB and IRB.
		 * Modification for nuclei and cosmological evolution of the photon background
		 * is considered with (cf. 10.1016/j.astropartphys.2012.07.010, eq. 3 and 5)
		 * beta_A,Z(E) = Z^2 / A * beta_p(E/A)
		 * beta(E,z) = (1+z)^3 beta((1+z)E).
		 */
		double lossLength(int id, double lf, double z = 0) const;
		
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_ELECTRONPAIRPRODUCTION_H
