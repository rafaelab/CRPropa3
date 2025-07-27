#include "crpropa/module/ElectronPairProduction.h"



namespace crpropa {

ElectronPairProduction::ElectronPairProduction(ref_ptr<PhotonField> photonField, bool electrons, ref_ptr<SamplerEvents> samplingEvts, ref_ptr<SamplerDistribution> samplingDist, int nSamples, double limit) {
	setInteractionTag("EPP");
	setPhotonField(photonField);
	setHaveElectrons(electrons);
	setLimit(limit);
	setSamplerEvents(samplingEvts);
	setSamplerDistribution(samplingDist);
	setMaximumSamples(nSamples); // default to no sampling
}

void ElectronPairProduction::setPhotonField(ref_ptr<PhotonField> field) {
	photonField = field;
	std::string fname = photonField->getFieldName();
	setDescription("ElectronPairProduction: " + fname);
	initRate(getDataPath("ElectronPairProduction/lossrate_" + fname + ".txt"));
	if (haveElectrons) { // Load secondary spectra only if electrons should be produced
		initSpectrum(getDataPath("ElectronPairProduction/spectrum_" + fname.substr(0,3) + ".txt"));
	}
}

void ElectronPairProduction::setHaveElectrons(bool electrons) {
	haveElectrons = electrons;
	if (haveElectrons) { // Load secondary spectra in case haveElectrons was changed to true
		std::string fname = photonField->getFieldName();
		initSpectrum(getDataPath("ElectronPairProduction/spectrum_" + fname.substr(0,3) + ".txt"));
	}
}

void ElectronPairProduction::setLimit(double lim) {
	limit = lim;
}

void ElectronPairProduction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

void ElectronPairProduction::setSamplerEvents(ref_ptr<SamplerEvents> s) {
	samplerEvents = s;
}

void ElectronPairProduction::setSamplerDistribution(ref_ptr<SamplerDistribution> s) {
	samplerDistribution = s;
}

void ElectronPairProduction::setMaximumSamples(int nSamples) {
	maximumSamples = nSamples;
}

std::string ElectronPairProduction::getInteractionTag() const {
	return interactionTag;
}

void ElectronPairProduction::initRate(std::string filename) {
	std::ifstream infile(filename.c_str());

	if (! infile.good())
		throw std::runtime_error("ElectronPairProduction: could not open file " + filename);

	// clear previously loaded interaction rates
	tabLorentzFactor.clear();
	tabLossRate.clear();

	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabLorentzFactor.push_back(pow(10, a));
				tabLossRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits <std::streamsize>::max(), '\n');
	}
	infile.close();
}

void ElectronPairProduction::initSpectrum(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (! infile.good())
		throw std::runtime_error("ElectronPairProduction: could not open file " + filename);

	double dNdE;
	tabSpectrum.resize(70);
	for (size_t i = 0; i < 70; i++) {
		tabSpectrum[i].resize(170);
		for (size_t j = 0; j < 170; j++) {
			infile >> dNdE;
			tabSpectrum[i][j] = dNdE * pow(10, (7 + 0.1 * j)); // read electron distribution pdf(Ee) ~ dN/dEe * Ee
		}
		for (size_t j = 1; j < 170; j++) {
			tabSpectrum[i][j] += tabSpectrum[i][j - 1]; // cdf(Ee), unnormalized
		}
	}
	infile.close();
}

double ElectronPairProduction::lossLength(int id, double lf, double z) const {
	double Z = chargeNumber(id);
	if (Z == 0)
		return std::numeric_limits<double>::max(); // no pair production on uncharged particles

	lf *= (1 + z);
	if (lf < tabLorentzFactor.front())
		return std::numeric_limits<double>::max(); // below energy threshold

	double rate;
	if (lf < tabLorentzFactor.back())
		rate = interpolate(lf, tabLorentzFactor, tabLossRate); // interpolation
	else
		rate = tabLossRate.back() * pow(lf / tabLorentzFactor.back(), -0.6); // extrapolation

	double A = nuclearMass(id) / mass_proton; // more accurate than massNumber(Id)
	rate *= Z * Z / A * pow_integer<3>(1 + z) * photonField->getRedshiftScaling(z);

	return 1. / rate;
}

void ElectronPairProduction::process(Candidate* c) const {
	int id = c->current.getId();
	if (not isNucleus(id))
		return; // only nuclei

	double lf = c->current.getLorentzFactor();
	double z = c->getRedshift();
	double E0 = c->current.getEnergy();
	double losslen = lossLength(id, lf, z);  // energy loss length
	if (losslen >= std::numeric_limits<double>::max())
		return;

	double step = c->getCurrentStep() / (1 + z); // step size in local frame
	double loss = step / losslen;  // relative energy loss

	if (haveElectrons) {
		double dE = c->current.getEnergy() * loss;  // energy loss
		int i = round((log10(lf) - 6.05) * 10);  // find closest cdf(Ee|log10(gamma))
		i = std::min(std::max(i, 0), 69);
		Random& random = Random::instance();

		double dE0 = dE;

		// draw pairs as long as their energy is smaller than the pair production energy loss
		int counter = 0;
		while (dE > 0) {
			// randomly draw an energy of the electron
			size_t j = random.randBin(tabSpectrum[i]);
			double Ee = pow(10, 6.95 + (j + random.rand()) * 0.1) * eV;

			 //  electron and positron in general don't have same lab frame energy, but averaged over many draws the result is consistent
			double Epair = 2 * Ee;

			// if the remaining energy is not sufficient check for random accepting
			if (Epair > dE)
				if (random.rand() > (dE / Epair))
					break; // not accepted

			// create pair and repeat with remaining energy
			dE -= Epair;

			if (samplerDistribution == nullptr or maximumSamples <= 0) {
				double Es = Epair / 2;
				double wp = samplerEvents->computeWeight(-11, Es, Es / E0, counter);
				double we = samplerEvents->computeWeight( 11, Es, Es / E0, counter);
				Vector3d pos = random.randomInterpolatedPosition(c->previous.getPosition(), c->current.getPosition());
				if (wp > 0)
					c->addSecondary(-11, Es, pos, wp, interactionTag);
				if (we > 0)
					c->addSecondary( 11, Es, pos, we, interactionTag);
			} else {
				samplerDistribution->push(Epair / 2.);
			}

			if (maximumSamples > 0 and counter > maximumSamples)
				break;

			counter++;
		}

		if (samplerDistribution != nullptr and maximumSamples > 0) {
			samplerDistribution->transformToPDF();
			samplerDistribution->transformToCDF();
			std::vector<double> sampledElectrons = samplerDistribution->getSample(maximumSamples);
			double dEs = std::accumulate(sampledElectrons.begin(), sampledElectrons.end(), decltype(sampledElectrons)::value_type(0)) * 2;
			if (samplerDistribution->getSize() > 0) {
				double w = dE0 / dEs / 2.;
				for (size_t i = 0; i < sampledElectrons.size(); i++) {
					double Es = sampledElectrons[i];
					
					// // this was not thoroughly tested yet, although it should work
					// double wp = w * samplerEvents->computeWeight(-11, Es, Es / E0, i);
					// double we = w * samplerEvents->computeWeight( 11, Es, Es / E0, i);

					double wp = 1;
					double we = 1;

					Vector3d pos = random.randomInterpolatedPosition(c->previous.getPosition(), c->current.getPosition());
					if (wp > 0)
						c->addSecondary(-11, Es, pos, wp, interactionTag);
					if (we > 0)
						c->addSecondary( 11, Es, pos, we, interactionTag);
				}
				samplerDistribution->clear();
			}
		}
	
	}



	c->current.setLorentzFactor(lf * (1 - loss));
	c->limitNextStep(limit * losslen);
}


} // namespace crpropa
