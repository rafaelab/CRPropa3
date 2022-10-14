#include <crpropa/Sampler.h>


namespace crpropa {


/****************************************************************************/

SamplerUniform::SamplerUniform(int pId, double s) {
	setSampling(s);
	setParticleId(pId);
}

void SamplerUniform::setSampling(double s) {
	sampling = s;
}

void SamplerUniform::setParticleId(int pId) {
	particleId = pId;
}

double SamplerUniform::getSampling() const {
	return sampling;
}

int SamplerUniform::getParticleId() const {
	return particleId;
}

double SamplerUniform::computeWeight(int id, double E, double f) const {
	if (id != particleId)
		return 0;

	if (sampling == 1.) {
		return 1;
	} else if (sampling == 0) {
		return 0;
	} else {
		Random &random = Random::instance();
		if (random.rand() < sampling) // accept and return weight
			return 1 / sampling;
		else // reject
			return 0;
	}
}


/****************************************************************************/

SamplerEnergyFraction::SamplerEnergyFraction(int pId, double s) {
	setSampling(s);
	setParticleId(pId);
}

void SamplerEnergyFraction::setSampling(double s) {
	sampling = s;
}

void SamplerEnergyFraction::setParticleId(int pId) {
	particleId = pId;
}

double SamplerEnergyFraction::getSampling() const {
	return sampling;
}

int SamplerEnergyFraction::getParticleId() const {
	return particleId;
}

double  SamplerEnergyFraction::computeWeight(int id, double E, double f) const {
	if (id != particleId)
		return 0;

	if (sampling == 1) {
		return 1;
	} else if (sampling == 0) {
		return 0;
	} else {
		Random &random = Random::instance();
		if (random.rand() < pow(f, 1 - sampling)) // accept and return weight
			return 1 / pow(f, 1 - sampling);
		else // reject
			return 0;
	}
}

/****************************************************************************/
SamplerNull::SamplerNull() {
}

double SamplerNull::computeWeight(int id, double E, double f) const {
	return 1.;
}

/****************************************************************************/

SamplerList::SamplerList() {
}

SamplerList::SamplerList(std::vector<ref_ptr<Sampler> > s) {
	for (size_t i = 0; i < s.size(); i++) {
		add(s[i]);
	}
}

void SamplerList::add(Sampler *sampler) {
	samplers.push_back(sampler);
}

double SamplerList::computeWeight(int id, double E, double f) const {
	double w = 1.;
	for (size_t i = 0; i < samplers.size(); i++) {
		w *= samplers[i]->computeWeight(id, E, f);
	}
	
	return w;
}


/****************************************************************************/

} // namespace crpropa