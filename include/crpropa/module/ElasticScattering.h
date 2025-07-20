#ifndef CRPROPA_ELASTICSCATTERING_H
#define CRPROPA_ELASTICSCATTERING_H


#include <cmath>
#include <vector>

#include "crpropa/Common.h"
#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Random.h"
#include "crpropa/Sampler.h"
#include "crpropa/Units.h"


namespace crpropa {

/**
 @class ElasticScattering
 @brief Elastic scattering of background photons on cosmic-ray nuclei.
 */
class ElasticScattering: public Module {
	protected:
		static constexpr double lgmin = 6.;  // minimum log10(Lorentz-factor)
		static constexpr double lgmax = 14.; // maximum log10(Lorentz-factor)
		static constexpr size_t nlg = 201;   // number of Lorentz-factor steps
		static constexpr double epsmin = -15.4943;    // [log10(2 * eV) + 3]; log10 minimum photon background energy in nucleus rest frame for elastic scattering
		static constexpr double epsmax = -10.3743; // [log10(2 * eV) + 8.12]; log10 maximum photon background energy in nucleus rest frame for elastic scattering
		static constexpr size_t neps = 513; // number of photon background energies in nucleus rest frame

	private:
		ref_ptr<PhotonField> photonField;
		ref_ptr<SamplerEvents> sampler; // sampler for the interaction rate (thinning)
		std::vector<double> tabRate; // elastic scattering rate
		std::vector<std::vector<double>> tabCDF; // CDF as function of background photon energy
		std::string interactionTag;

	public:
		/** Constructor
		@param photonField		target photon field
		@param sampler		    sampling object (see Sampler.h)
		*/
		ElasticScattering(ref_ptr<PhotonField> photonField, ref_ptr<SamplerEvents> sampling = ref_ptr<SamplerEvents>(new SamplerEventsNull()));

		void setPhotonField(ref_ptr<PhotonField> photonField);
		void setInteractionTag(std::string tag);
		void setSampler(ref_ptr<SamplerEvents> sampler);

		std::string getInteractionTag() const;

		void initRate(std::string filename);
		void initCDF(std::string filename);
		void process(Candidate* candidate) const;

};

} // namespace crpropa

#endif // CRPROPA_ELASTICSCATTERING_H
