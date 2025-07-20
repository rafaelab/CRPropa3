#include "crpropa/module/ElasticScattering.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {


ElasticScattering::ElasticScattering(ref_ptr<PhotonField> field, ref_ptr<SamplerEvents> sampling) {
	setInteractionTag("ES");
	setPhotonField(field);
	setSampler(sampling);
}

void ElasticScattering::setPhotonField(ref_ptr<PhotonField> field) {
	photonField = field;
	std::string fname = photonField->getFieldName();
	setDescription("ElasticScattering: " + fname);
	initRate(getDataPath("ElasticScattering/rate_" + fname.substr(0, 3) + ".txt"));
	initCDF(getDataPath("ElasticScattering/cdf_" + fname.substr(0, 3) + ".txt"));
}


void ElasticScattering::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

void ElasticScattering::setSampler(ref_ptr<SamplerEvents> s) {
	sampler = s;
}

std::string ElasticScattering::getInteractionTag() const {
	return interactionTag;
}

void ElasticScattering::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("ElasticScattering: could not open file " + filename);

	tabRate.clear();

	while (infile.good()) {
		if (infile.peek() == '#') {
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		double r;
		infile >> r;
		if (! infile)
			break;
		tabRate.push_back(r / Mpc);
	}

	infile.close();
}

void ElasticScattering::initCDF(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("ElasticScattering: could not open file " + filename);

	tabCDF.clear();
	std::string line;
	double a;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;

		std::stringstream lineStream(line);
		lineStream >> a;

		std::vector<double> cdf(neps);
		for (size_t i = 0; i < neps; i++) {
			lineStream >> a;
			cdf[i] = a;
		}
		tabCDF.push_back(cdf);
	}

	infile.close();
}

void ElasticScattering::process(Candidate* candidate) const {
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double E0 = candidate->current.getEnergy();

	if (not isNucleus(id))
		return;

	double lg = log10(candidate->current.getLorentzFactor() * (1 + z));
	if ((lg < lgmin) or (lg > lgmax))
		return;

	int A = massNumber(id);
	int Z = chargeNumber(id);
	int N = A - Z;

	double step = candidate->getCurrentStep();
	while (step > 0) {

		double rate = interpolateEquidistant(lg, lgmin, lgmax, tabRate);
		rate *= Z * N / double(A);  // TRK scaling
		rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);  // cosmological scaling

		// check for interaction
		Random& random = Random::instance();
		double randDist = -log(random.rand()) / rate;
		if (step < randDist)
			return;

		// draw random background photon energy from CDF
		size_t i = floor((lg - lgmin) / (lgmax - lgmin) * (nlg - 1)); // index of closest gamma tabulation point
		size_t j = random.randBin(tabCDF[i]) - 1; // index of next lower tabulated eps value
		double binWidth = (epsmax - epsmin) / (neps - 1); // logarithmic bin width
		double eps = pow(10, epsmin + (j + random.rand()) * binWidth);

		// boost to lab frame
		double cosTheta = 2 * random.rand() - 1;
		double E = eps * candidate->current.getLorentzFactor() * (1. - cosTheta);

		Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
		double w = sampler->computeWeight(id, E / (1 + z), E / E0, 0);
		if (w > 0)
			candidate->addSecondary(22, E, pos, w, interactionTag);

		// repeat with remaining step
		step -= randDist;
	}
}


} // namespace crpropa
