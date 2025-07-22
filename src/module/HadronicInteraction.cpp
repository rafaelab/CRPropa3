#include "crpropa/module/HadronicInteraction.h"



namespace crpropa {


HadronicInteraction::CrossSection::CrossSection() { 
}

HadronicInteraction::CrossSection::CrossSection(std::string file) { 
	loadData(file);
}

HadronicInteraction::CrossSection::CrossSection(std::string file, int id, double loss, bool add, int CR_Id) : nucleusId(id), weightEnergyLoss(loss), addSecondary(add), CR_Id(CR_Id) {
	loadData(file);
}

void HadronicInteraction::CrossSection::loadData(std::string file) { 
	std::ifstream infile(file.c_str());

	if (! infile.good()) 
		throw std::runtime_error("CrossSection::loadData could not open file " + file);

	// clear old data 
	sigma.clear();
	energyPrimary.clear();
	energySecondary.clear();

	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	// read secondary energy as the first row
	double e; 
	infile >> e; // ignore the first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> e;
		energySecondary.push_back(e + secondary_mc2.at(nucleusId)); // we would like this energy to become the total energy of the secondary here
	}	

	// read the following lines Tp, sigma values
	double a; // helper variable to store 
	while (infile.good()) {
		infile >> a;
		if (! infile) 
			break;	// end of file

		energyPrimary.push_back(a);

		std::vector<double> cdf;
		for (int i = 0; i < energySecondary.size(); i++){
			infile >> a;
			cdf.push_back(a);
		}
		sigma.push_back(cdf);
	}

	// finally close the file
	infile.close();
	
}

double HadronicInteraction::CrossSection::totalInelasticCrossSectionProton(const double& energy) {
	// Parametrisation from:
	//   Kafexhiu et al. PRD 90 (2014) 123014.
	// Note: only works for en >> Ethr
	double ethr = 2.797e8 * eV;
	double x = (energy - mass_proton * c_squared) / ethr;
	double sigma = (30.7 - 0.96 * log(x) + 0.18 * log(x) * log(x)) * pow(1. - pow(1. / x, 1.9), 3);
	sigma *= (milli * barn); 
	return sigma;
}

double HadronicInteraction::CrossSection::totalInelasticCrossSection(const double& T, const int& id, const int& At) { 
	double sigmaPP = totalInelasticCrossSectionProton(T); 
	double sigmaPP0 = totalInelasticCrossSectionProton(1 * TeV);

	// mass number of the projectile
	int Ap = massNumber(id); 

	double beta0;
	if (Ap == 1) {
		beta0 = 2.247 - 0.915 * (1 + pow(At, - 1 / 3));
	} else {
		beta0 = 1.581 - 0.876 * (pow(Ap, - 1 / 3) + pow(At, - 1 / 3));
	}

	// follows Sihver et al. PRC 47 (1993) 1225
	// also eq. 17 of Kafexhiu et al.
	double sigmaR = M_PI * pow_integer<2>(radius_electron) * pow_integer<2>((pow(Ap, 1. / 3) + pow(At, 1. / 3) - beta0 * (pow(Ap, -1. / 3) + pow(At, -1. / 3))));

	// eq. 18 of Kafexhiu et al.
	double G = 1 + std::log(std::max(1.0, sigmaPP / sigmaPP0)); 

	return sigmaR * G;
}

///////////////////////////////////////////////////////////////////////////////

HadronicInteraction::HadronicInteraction(ref_ptr<Density> density, bool havePhotons, bool haveElectrons, bool haveNeutrinos, bool haveAntiNucleons, bool deactivateAfter, ref_ptr<SamplerEvents> sampling, double limit) {
	setInteractionTag("HI");
	setDescription("HadronicInteraction");
	setDensity(density);
	setHavePhotons(havePhotons);
	setHaveElectrons(haveElectrons);
	setHaveNeutrinos(haveNeutrinos);
	setHaveAntiNucleons(haveAntiNucleons);
	setLimit(limit);
	setDeactivatePrimary(deactivateAfter);
	setSampler(sampling);
}

HadronicInteraction::HadronicInteraction(ref_ptr<Density> density, std::string configFile, bool deactivateAfter, ref_ptr<SamplerEvents> sampling, double limit) {
	setInteractionTag("HI");
	setDescription("HadronicInteraction");
	setDensity(density);
	setLimit(limit);
	setDeactivatePrimary(deactivateAfter);
	setSampler(sampling);
	if (configFile != "")
		initData(configFile);
}

void HadronicInteraction::loadData(std::string file, HadronicInteraction::CrossSection& cs) {
	// load data for a single cross section
	std::ifstream infile(file.c_str());

	if (! infile.good()) 
		throw std::runtime_error("HadronicInteraction: could not open file " + file);

	// clear old data 
	cs.sigma.clear();

	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	// read secondary energy as the first row
	double e; 
	infile >> e; // ignore the first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> e;
		try {
			cs.energySecondary.push_back(e + cs.secondary_mc2.at(cs.nucleusId));
		} catch (std::out_of_range& e) {
			KISS_LOG_ERROR << "HadronicInteraction::loadData: nucleusId " << cs.nucleusId << " not found in secondary_mc2 map." << std::endl;
			throw std::runtime_error("HadronicInteraction: nucleusId not found in secondary_mc2 map");
		}
	}	

	double a;
	while(infile.good()) {
		infile >> a;
		if(! infile) 
			break;
		
		cs.energyPrimary.push_back(a);
	
		std::vector<double> cdf;
		for (int i = 0; i < cs.energySecondary.size(); i++){
			infile >> a;
			cdf.push_back(a);
		}
		cs.sigma.push_back(cdf);
	}

	infile.close();
}

void HadronicInteraction::initData(std::string configFile) {
	std::ifstream infile;
	if ((configFile.at(0) == '/') or (configFile.at(0) == '.')) 
		infile.open(configFile.c_str());
	else 
		infile.open(getDataPath("hadronic_interaction/" + configFile).c_str()); // use path in data folder

	if (! infile.good()) 
		throw std::runtime_error("HadronicInteraction: could not open file " + configFile);

	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	while (infile.good()) {
		HadronicInteraction::CrossSection cs;
		infile >> cs.nucleusId; // rather the ID of the secondary
		infile >> cs.addSecondary;
		infile >> cs.weightEnergyLoss;
		infile >> cs.CR_Id;

		std::string path; 
		infile >> path;
		
		if(! infile) 
			break;	// end of file

		// decide if local or global pathes are used 
		std::string loadFile = path.substr(0,1); 
		if ((loadFile == "/") or (loadFile == "."))
			loadFile = path; // use absolute path
		else 
			loadFile = getDataPath("hadronic_interaction/" + path); // relative path in the data folder
		
		loadData(loadFile, cs);
		crossSectionList[cs.CR_Id].push_back(cs);
	}
}

void HadronicInteraction::addCrossSection(HadronicInteraction::CrossSection cs) { 
	// problem there : I would need the ID corresponding to the cs to add it at the right place to CSList ; maybe I can add the 
	crossSectionList.at(cs.CR_Id).push_back(cs);
}

void HadronicInteraction::setHavePhotons(bool b) {
	havePhotons = b;
}

void HadronicInteraction::setHaveElectrons(bool b) {
	haveElectrons = b;
}

void HadronicInteraction::setHaveNeutrinos(bool b) {
	haveNeutrinos = b;
}

void HadronicInteraction::setHaveAntiNucleons(bool b) {
	haveAntiNucleons = b;
}

void HadronicInteraction::setLimit(double lim) {
	limit = lim;
}

void HadronicInteraction::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

void HadronicInteraction::setDeactivatePrimary(bool b) {
	deactivatePrimary = b;
}

void HadronicInteraction::setDensity(ref_ptr<Density> d) {
	density = d;
}

void HadronicInteraction::setSampler(ref_ptr<SamplerEvents> s) {
	sampler = s;
}

double HadronicInteraction::getLimit() const {
	return limit;
}

bool HadronicInteraction::getDeactivatePrimary() const {
	return deactivatePrimary;
}

bool HadronicInteraction::getHavePhotons() const {
	return havePhotons;
}

bool HadronicInteraction::getHaveElectrons() const {
	return haveElectrons;
}

bool HadronicInteraction::getHaveNeutrinos() const {
	return haveNeutrinos;
}

bool HadronicInteraction::getHaveAntiNucleons() const {
	return haveAntiNucleons;
}

std::string HadronicInteraction::getInteractionTag() const {
	return interactionTag;
}

double HadronicInteraction::getDensityAtPosition(const Vector3d& pos) const {
	return density->getDensity(pos);
}

void HadronicInteraction::allowSecondaryId(int id, bool allow) {
	bool foundId = false;
	for(int i = 0; i < crossSectionList[id].size(); i++) {
		if(crossSectionList[id][i].nucleusId == id){
			foundId = true;
			crossSectionList[id][i].addSecondary = allow;
		}
	}

	if (! foundId) {
		KISS_LOG_WARNING << "HadronicInteraction:allowSecondaryId \tid " << id << " was not found in the list of cross section" << std::endl;  
		throw std::runtime_error("");
	}
}

bool HadronicInteraction::isPrimaryImplemented(const int& id) const {
	return crossSectionList.find(id) != crossSectionList.end();
}

void HadronicInteraction::process(Candidate* candidate) const {
	// candidate information 
	int id = candidate->current.getId();
	if (! isPrimaryImplemented(id))
		return;
	
	double z = candidate->getRedshift();
	double Eprimary = candidate->current.getEnergy() * (1 + z);
	double step = candidate->getCurrentStep() / (1 + z);
	Vector3d pos = candidate->current.getPosition();
	Vector3d dir = candidate->current.getDirection();


	Random& random = Random::instance();
	ref_ptr<MediumComposition> target = density->getTargetMedium();
	unsigned int nComponents = 1;
	if (target->isAdmixed())
		nComponents = target->numberOfComponents();

	for (size_t i = 0; i < nComponents; i++) {
		double At = target->getCompositionWeight();

		// decide on interaction 
		double n0 = density->getDensity(pos, z);
		double prob = CrossSection::totalInelasticCrossSection(Eprimary, id, At) * n0 * step;

		if (random.rand() < prob){
			performInteraction(candidate, At);
		} else {
			candidate->limitNextStep(limit * step / prob);
			double n1 = density->getDensity(pos + step * dir);
			double dNdX = fabs(n0 - n1) / n0 / step;
			candidate->limitNextStep(limit / dNdX);
		}
	}
}

void HadronicInteraction::performInteraction(Candidate* candidate, const int& At) const {
	
	bool activate = true;
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double Eprimary = candidate->current.getEnergy() * (1 + z);
	double K = Eprimary - CR_rest_energy.at(id); // kinetic energy of the primary in [J] ; Tp is the kinetic energy of the projectile
	double sigmaTot = CrossSection::totalInelasticCrossSection(K, id, At) / (milli * barn);
	double totalEnergyLoss = 0;
	double sampledEnergy = 0;

	std::vector<double> sigmaInterp;
	std::vector<double> dSigma{0, 0, 0};

	Random& random = Random::instance();

	// loop over all secondary species
	for (size_t i = 0; i < crossSectionList.at(id + At - 1).size(); i++) {
		sigmaInterp.clear();

		CrossSection cs = crossSectionList.at(id + At - 1)[i];
		int index = closestIndex(K / GeV, cs.energyPrimary); // then there is no interpolation there ; here Tp is total energy
		int l = 0;
		
		if (index == 0) {
			l = 1;
		} else if (index == cs.energyPrimary.size() - 1) {
			l = -1;
		}
		std::vector<double> subvector_Ep{cs.energyPrimary[index + l - 1], cs.energyPrimary[index + l], cs.energyPrimary[index + l + 1]};

		for (size_t k = 0; k < cs.sigma[index].size(); k++){
			dSigma[0] = cs.sigma[index + l - 1][k];
			dSigma[1] = cs.sigma[index + l][k];
			dSigma[2] = cs.sigma[index + l + 1][k];
			sigmaInterp.push_back(interpolate(K / GeV, subvector_Ep, dSigma));
		}
		
		// integrated differential crosssection #no interpolation # get the total integrated diff cross section
		double dSigma = sigmaInterp[cs.energySecondary.size() - 1];

		// number of produced secondaries (for a type of secondary)
		// draw a random bin from a (unnormalized) cumulative distribution function, 
		// int nSecondaries = randomRound(dsigma / sigmaTot);
		double sigmaRatio = dSigma / sigmaTot;
		int nSecondaries = static_cast<int>(std::floor(sigmaRatio));
		nSecondaries = sigmaRatio - nSecondaries;
		if (random.rand() > nSecondaries) {
			nSecondaries++;
		}
		

		// sample secondaries
		for(size_t j = 0; j < nSecondaries; j++) {
			 // draw a random bin from a (unnormalized) cumulative distribution function, without leading zero
			int bin = random.randBin(sigmaInterp);
			if (bin > 0) {
				 // from that, in the tables, E_sec should be total energy ; or I should add the mass energy ; also I should not take the middle of the bins
				sampledEnergy = cs.energySecondary[bin - 1] + random.rand() * (cs.energySecondary[bin] - cs.energySecondary[bin - 1]);
			} else {
				// first possible value for us is 1 GeV
				sampledEnergy = 1 + random.rand() * (cs.energySecondary[bin] - 1);
			}
			totalEnergyLoss += sampledEnergy * cs.weightEnergyLoss;

			// add secondary if allowed
			if(cs.addSecondary) {
				candidate->addSecondary(cs.nucleusId, sampledEnergy * GeV, 1., interactionTag);
				if ((cs.nucleusId == id) and ((id == 1000000010) or (id == 1000010010))){
					activate = false;
				}
			}
		}
	}

	// apply energy loss
	candidate->current.setEnergy(K - totalEnergyLoss * GeV);
	candidate->setActive(activate);
	if (deactivatePrimary)
		candidate->setActive(false);
		
}



} // namespace crpropa 
