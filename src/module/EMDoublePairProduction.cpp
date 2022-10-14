#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

EMDoublePairProduction::EMDoublePairProduction(ref_ptr<PhotonField> photonField, bool haveElectrons, ref_ptr<Sampler> sampling, double limit) {
	setPhotonField(photonField);
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setSampler(sampling);
}

void EMDoublePairProduction::setPhotonField(ref_ptr<PhotonField> field) {
	photonField = field;
	std::string fname = photonField->getFieldName();
	setDescription("EMDoublePairProduction: " + fname);
	initRate(getDataPath("EMDoublePairProduction/rate_" + fname + ".txt"));
}

void EMDoublePairProduction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void EMDoublePairProduction::setLimit(double l) {
	limit = l;
}

void EMDoublePairProduction::setSampler(ref_ptr<Sampler> s) {
	if (s == NULL)
		sampler = new SamplerNull();
	else
		sampler = s;
}

void EMDoublePairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (!infile.good())
		throw std::runtime_error("EMDoublePairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabEnergy.clear();
	tabRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(pow(10, a) * eV);
				tabRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void EMDoublePairProduction::performInteraction(Candidate *candidate) const {
	// the photon is lost after the interaction
	candidate->setActive(false);

	if (not haveElectrons)
		return;

	// Use assumption of Lee 96 arXiv:9604098
	// Energy is equally shared between one e+e- pair, but take mass of second e+e- pair into account.
	// This approximation has been shown to be valid within -1.5%.
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double Ee = (E - 2 * mass_electron * c_squared) / 2;
	double f = Ee / E;

	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());

	// add secondaries
	double we = sampler->computeWeight( 11, Ee / (1 + z), 1 - f);
	double wp = sampler->computeWeight(-11, Ee / (1 + z), f);
	if (we > 0)
		candidate->addSecondary( 11, Ee / (1 + z), pos, we);
	if (wp > 0)
		candidate->addSecondary(-11, Ee / (1 + z), pos, wp);
}

void EMDoublePairProduction::process(Candidate *candidate) const {
	// check if photon
	if (candidate->current.getId() != 22)
		return;

	// scale the electron energy instead of background photons
	double z = candidate->getRedshift();
	double E = (1 + z) * candidate->current.getEnergy();

	// check if in tabulated energy range
	if (E < tabEnergy.front() or (E > tabEnergy.back()))
		return;

	// interaction rate
	double rate = interpolate(E, tabEnergy, tabRate);
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

	// check for interaction
	Random &random = Random::instance();
	double randDistance = -log(random.rand()) / rate;
	double step = candidate->getCurrentStep();
	if (step < randDistance) {
		candidate->limitNextStep(limit / rate);
		return;
	} else { // after performing interaction photon ceases to exist (hence return)
		performInteraction(candidate);
		return;
	}

}


} // namespace crpropa
