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

double SamplerUniform::computeWeight(int id, double E, double f, int counter) const {
	if (id != particleId)
		return 0;

	if (sampling >= 1.) {
		return 1;
	} else if (sampling <= 0) {
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

SamplerEnergyFractionPowerLaw::SamplerEnergyFractionPowerLaw(double idx, int pId, double s) {
	setSampling(s);
	setParticleId(pId);
	setIndex(idx);
}

void SamplerEnergyFractionPowerLaw::setSampling(double s) {
	sampling = s;
}

void SamplerEnergyFractionPowerLaw::setParticleId(int pId) {
	particleId = pId;
}

void SamplerEnergyFractionPowerLaw::setIndex(double idx) {
	index = idx;
}

double SamplerEnergyFractionPowerLaw::getSampling() const {
	return sampling;
}

int SamplerEnergyFractionPowerLaw::getParticleId() const {
	return particleId;
}

double SamplerEnergyFractionPowerLaw::getIndex() const {
	return index;
}

double  SamplerEnergyFractionPowerLaw::computeWeight(int id, double E, double f, int counter) const {
	if (id != particleId)
		return 0;

	if (sampling >= 1) {
		return 1;
	} else if (sampling < 0) {
		return 0;
	} else {
		Random &random = Random::instance();
		double r = 0;
		if (index < 0)
			r = pow(1 - f, (1 - sampling) * abs(index));
		else
			r = pow(f, (1 - sampling) * index);
		if (random.rand() < r) // accept and return weight
			return 1. / r;
		else // reject
			return 0;
	}
}


/****************************************************************************/
SamplerNull::SamplerNull() {
}

double SamplerNull::computeWeight(int id, double E, double f, int counter) const {
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

double SamplerList::computeWeight(int id, double E, double f, int counter) const {
	double w = 0;
	int k = 0;
	for (size_t i = 0; i < samplers.size(); i++) {
		double w0 = samplers[i]->computeWeight(id, E, f, counter);
		if (w0 > 0) {
			if (k == 0)
				w = 1;
			w *= w0;
			k++;
		}
	}
	
	return w;
}


/****************************************************************************/

} // namespace crpropa