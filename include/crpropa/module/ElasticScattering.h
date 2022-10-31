#ifndef CRPROPA_ELASTICSCATTERING_H
#define CRPROPA_ELASTICSCATTERING_H


#include <cmath>
#include <limits>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "crpropa/Module.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
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
private:
	ref_ptr<PhotonField> photonField;
	ref_ptr<SamplerEvents> sampler;

	std::vector<double> tabRate; // elastic scattering rate
	std::vector<std::vector<double> > tabCDF; // CDF as function of background photon energy

	static const double lgmin; // minimum log10(Lorentz-factor)
	static const double lgmax; // maximum log10(Lorentz-factor)
	static const size_t nlg;   // number of Lorentz-factor steps
	static const double epsmin; // minimum log10(eps / J)
	static const double epsmax; // maximum log10(eps / J)
	static const size_t neps;   // number of eps steps

public:
	/** Constructor
	 @param photonField		target photon field
	 */
	ElasticScattering(ref_ptr<PhotonField> photonField, ref_ptr<SamplerEvents> sampler = NULL);
	void initRate(std::string filename);
	void initCDF(std::string filename);
	void setPhotonField(ref_ptr<PhotonField> photonField);
	void setSampler(ref_ptr<SamplerEvents> sampler);
	void process(Candidate *candidate) const;
};

} // namespace crpropa

#endif // CRPROPA_ELASTICSCATTERING_H
