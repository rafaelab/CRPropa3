#include <crpropa/Sampler.h>


namespace crpropa {


/****************************************************************************/

SamplerEventsEnergy::SamplerEventsEnergy(int pId, double s) {
	setSampling(s);
	setParticleId(pId);
}

void SamplerEventsEnergy::setSampling(double s) {
	sampling = s;
}

void SamplerEventsEnergy::setParticleId(int pId) {
	particleId = pId;
}

double SamplerEventsEnergy::getSampling() const {
	return sampling;
}

int SamplerEventsEnergy::getParticleId() const {
	return particleId;
}

double  SamplerEventsEnergy::computeWeight(int id, double E, double f, int counter) const {
	if (id != particleId)
		return 0;

	if (sampling >= 1) {
		return 1;
	} else if (sampling < 0) {
		return 0;
	} else {
		Random &random = Random::instance();
		double r = weightFunction(id, E, f, counter);
		if (random.rand() < r) // accept and return weight
			return 1. / r;
		else // reject
			return 0;
	}

}


/****************************************************************************/

SamplerEventsUniform::SamplerEventsUniform(int pId, double s) {
	setSampling(s);
	setParticleId(pId);
}

void SamplerEventsUniform::setSampling(double s) {
	sampling = s;
}

void SamplerEventsUniform::setParticleId(int pId) {
	particleId = pId;
}

double SamplerEventsUniform::getSampling() const {
	return sampling;
}

int SamplerEventsUniform::getParticleId() const {
	return particleId;
}

double SamplerEventsUniform::computeWeight(int id, double E, double f, int counter) const {
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

SamplerEventsEnergyFractionPowerLaw::SamplerEventsEnergyFractionPowerLaw(double idx, int pId, double s) : SamplerEventsEnergy(pId, s) {
	setIndex(idx);
}

void SamplerEventsEnergyFractionPowerLaw::setIndex(double idx) {
	index = idx;
}

double SamplerEventsEnergyFractionPowerLaw::getIndex() const {
	return index;
}

double SamplerEventsEnergyFractionPowerLaw::weightFunction(int id, double E, double f, int counter) const {
	return pow(f, (1 - sampling) * index);
}


/****************************************************************************/
SamplerEventsNull::SamplerEventsNull() {
}

double SamplerEventsNull::computeWeight(int id, double E, double f, int counter) const {
	return 1.;
}

/****************************************************************************/

SamplerEventsList::SamplerEventsList() {
}

SamplerEventsList::SamplerEventsList(std::vector<ref_ptr<SamplerEvents>> s) {
	for (size_t i = 0; i < s.size(); i++) {
		add(s[i]);
	}
}

void SamplerEventsList::add(SamplerEvents *SamplerEvents) {
	samplers.push_back(SamplerEvents);
}

double SamplerEventsList::computeWeight(int id, double E, double f, int counter) const {
	if (samplers.size() == 0)
		return 1.;

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

SamplerDistributionUniform::SamplerDistributionUniform(double vmin, double vmax, int nBins, std::string scale) {
	ref_ptr<Histogram1D> h = new Histogram1D(vmin, vmax, nBins, scale);	
	setSize(0);
	setDistribution(h);
}


void SamplerDistributionUniform::setDistribution(ref_ptr<Histogram1D> dist) {
	distribution = dist;
}

void SamplerDistributionUniform::setSize(int n) {
	datasetSize = n;
}


ref_ptr<Histogram1D> SamplerDistributionUniform::getDistribution() const {
	return distribution;
}

int SamplerDistributionUniform::getSize() const {
	return datasetSize;
}

std::vector<double> SamplerDistributionUniform::getSample(int nSamples) const {
	std::vector<double> sample;
	for (size_t i = 0; i < std::min(nSamples, distribution->getNumberOfBins()); i++) {
		sample.push_back(distribution->getSample());
	}

	return sample;
}

void SamplerDistributionUniform::transformToPDF() {
	distribution->transformToPDF();
}


void SamplerDistributionUniform::transformToCDF() {
	distribution->transformToCDF();
}

void SamplerDistributionUniform::append(const std::vector<double> &v) {
	for (size_t i = 0; i < v.size(); i++) {
		distribution->push(v[i]);
	}
	datasetSize = v.size();
}

void SamplerDistributionUniform::push(const double &v) {
	distribution->push(v);
	datasetSize++;
}

void SamplerDistributionUniform::clear() {
	distribution->clear();
	datasetSize = 0;
}



/****************************************************************************/

} // namespace crpropa