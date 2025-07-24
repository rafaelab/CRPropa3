#include "crpropa/module/SynchrotronRadiation.h"


namespace crpropa {

SynchrotronRadiation::SynchrotronRadiation(ref_ptr<MagneticField> field, bool havePhotons, ref_ptr<SamplerEvents> samplingEvts, ref_ptr<SamplerDistribution> samplingDist, int nSamples, double limit) {
	setInteractionTag("SYN");
	setField(field);
	setBrms(0);
	initSpectrum();
	setHavePhotons(havePhotons);
	setLimit(limit);
	setSecondaryThreshold(1e6 * eV);
	setMaximumSamples(nSamples);
	setSamplerEvents(samplingEvts);
	setSamplerDistribution(samplingDist);
}

SynchrotronRadiation::SynchrotronRadiation(double Brms, bool havePhotons, ref_ptr<SamplerEvents> samplingEvts, ref_ptr<SamplerDistribution> samplingDist, int nSamples, double limit) {
	setInteractionTag("SYN");
	setBrms(Brms);
	initSpectrum();
	setHavePhotons(havePhotons);
	setLimit(limit);
	setSecondaryThreshold(1e6 * eV);
	setMaximumSamples(nSamples);
	setSamplerEvents(samplingEvts);
	setSamplerDistribution(samplingDist);
}

void SynchrotronRadiation::setField(ref_ptr<MagneticField> f) {
	field = f;
}

ref_ptr<MagneticField> SynchrotronRadiation::getField() const {
	return field;
}

void SynchrotronRadiation::setBrms(double B) {
	Brms = B;
}

double SynchrotronRadiation::getBrms() const {
	return Brms;
}

void SynchrotronRadiation::setHavePhotons(bool photons) {
	havePhotons = photons;
}

bool SynchrotronRadiation::getHavePhotons() const {
	return havePhotons;
}


void SynchrotronRadiation::setLimit(double lim) {
	limit = lim;
}

double SynchrotronRadiation::getLimit() const {
	return limit;
}

void SynchrotronRadiation::setMaximumSamples(int nSamples) {
	maximumSamples = nSamples;
}

int SynchrotronRadiation::getMaximumSamples() const {
	return maximumSamples;
}


void SynchrotronRadiation::setSamplerEvents(ref_ptr<SamplerEvents> s) {
	samplerEvents = s;
}

ref_ptr<SamplerEvents> SynchrotronRadiation::getSamplerEvents() const {
	return samplerEvents;
}

void SynchrotronRadiation::setSamplerDistribution(ref_ptr<SamplerDistribution> s) {
	samplerDistribution = s;
}

ref_ptr<SamplerDistribution> SynchrotronRadiation::getSamplerDistribution() const {
	return samplerDistribution;
}

void SynchrotronRadiation::setSecondaryThreshold(double threshold) {
	secondaryThreshold = threshold;
}

double SynchrotronRadiation::getSecondaryThreshold() const {
	return secondaryThreshold;
}

void SynchrotronRadiation::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string SynchrotronRadiation::getInteractionTag() const {
	return interactionTag;
}


void SynchrotronRadiation::initSpectrum() {
	std::string filename = getDataPath("Synchrotron/spectrum.txt");
	std::ifstream infile(filename.c_str());

	if (! infile.good())
		throw std::runtime_error("SynchrotronRadiation: could not open file " + filename);

	// clear previously loaded interaction rates
	tabx.clear();
	tabCDF.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabx.push_back(pow(10, a));
				tabCDF.push_back(b);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	}
	infile.close();
}

void SynchrotronRadiation::process(Candidate* candidate) const {
	double charge = fabs(candidate->current.getCharge());
	if (charge == 0)
		return; // only charged particles

	// calculate gyroradius, evaluated at the current position
	double z = candidate->getRedshift();
	double B;
	if (field.valid()) {
		Vector3d Bvec = field->getField(candidate->current.getPosition(), z);
		B = Bvec.cross(candidate->current.getDirection()).getR();
	} else {
		B = sqrt(2. / 3) * Brms; // average perpendicular field component
	}
	B *= pow(1 + z, 2); // cosmological scaling
	double Rg = candidate->current.getMomentum().getR() / charge / B;

	// calculate energy loss
	double lf = candidate->current.getLorentzFactor();
	double dEdx = 1. / 6 / M_PI / epsilon0 * pow(lf * lf - 1, 2) * pow(charge / Rg, 2); // Jackson p. 770 (14.31)
	double step = candidate->getCurrentStep() / (1 + z); // step size in local frame
	double dE = step * dEdx;

	// apply energy loss and limit next step
	double E = candidate->current.getEnergy();
	candidate->current.setEnergy(E - dE);
	candidate->limitNextStep(limit * E / dEdx);

	// optionally add secondary photons
	if (havePhotons) {
		// check if photons with energies > 14 * Ecrit are possible
		double Ecrit = 3. / 4 * h_planck / M_PI * c_light * pow(lf, 3) / Rg;
		if (14 * Ecrit < secondaryThreshold)
			return;

		// draw photons up to the total energy loss
		// if maximumSamples is reached before that, compensate the total energy afterwards
		Random& random = Random::instance();
		double dE0 = dE;
		std::vector<double> energies;
		int counter = 0;
		while (dE > 0) {
			// draw random value between 0 and maximum of corresponding cdf
			// choose bin of s where cdf(x) = cdf_rand -> x_rand
			size_t i = random.randBin(tabCDF); // draw random bin (upper bin boundary returned)
			double binWidth = (tabx[i] - tabx[i - 1]);
			double x = tabx[i - 1] + random.rand() * binWidth; // draw random x uniformly distributed in bin
			double Ephoton = x * Ecrit;

			// if the remaining energy is not sufficient check for random accepting
			if (Ephoton > dE) {
				if (random.rand() > (dE / Ephoton))
					break; // not accepted
			}

			if (samplerDistribution == nullptr or maximumSamples <= 0) {
				double Es = Ephoton / 2;
				double w = samplerEvents->computeWeight(22, Es, Es / E, counter);
				Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
				if (w > 0)
					candidate->addSecondary(-11, Es, pos, w, interactionTag);
			} else {
				samplerDistribution->push(Ephoton / 2.);
			}

			if (maximumSamples > 0 and counter > maximumSamples)
				break;

			// energy loss
			dE -= Ephoton;

			// counter for sampling break condition;
			counter++;
		}

		if (samplerDistribution != nullptr and maximumSamples > 0) {
			samplerDistribution->transformToPDF();
			samplerDistribution->transformToCDF();
			std::vector<double> sampledPhotons = samplerDistribution->getSample(maximumSamples);
			double dEs = std::accumulate(sampledPhotons.begin(), sampledPhotons.end(), decltype(sampledPhotons)::value_type(0)) * 2;
			if (samplerDistribution->getSize() > 0) {
				double w = dE0 / dEs;
				for (size_t i = 0; i < sampledPhotons.size(); i++) {
					double Es = sampledPhotons[i];

					// this was not thoroughly tested yet, although it should work
					// w *= samplerEvents->computeWeight(-11, Es, Es / E0, i); 
					if (w > 0) {
						Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
						candidate->addSecondary(-11, Es, pos, w, interactionTag);
					}
				}
				samplerDistribution->clear();
			}
		}
	}
	
}

std::string SynchrotronRadiation::getDescription() const {
	std::stringstream s;
	s << "Synchrotron radiation";
	if (field.valid())
		s << " for specified magnetic field";
	else
		s << " for Brms = " << Brms / nG << " nG";
	if (havePhotons)
		s << ", synchrotron photons E > " << secondaryThreshold / eV << " eV";
	else
		s << ", no synchrotron photons";
	if (maximumSamples > 0)
		s << "maximum number of photon samples: " << maximumSamples;
	return s.str();
}


} // namespace crpropa
