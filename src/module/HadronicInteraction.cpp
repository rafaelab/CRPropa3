#include "crpropa/module/HadronicInteraction.h"


namespace crpropa {


HadronicInteraction::CrossSection::CrossSection() { 
}

HadronicInteraction::CrossSection::CrossSection(int primaryId, int secondaryId, ref_ptr<MediumComposition> target, double weightEnergyLoss) {
	setPrimaryId(primaryId);
	setSecondaryId(secondaryId);
	setTargetMedium(target);
	setWeightEnergyLoss(weightEnergyLoss);
	loadData();
}

void HadronicInteraction::CrossSection::loadData() { 

	std::string filename = getDataPath("HadronicInteraction/cumDiffSigma-prim_" + std::to_string(primaryId) + "-sec_" + std::to_string(secondaryId) + "-tgt_" + targetMedium->getLabel() + ".txt");

	std::filesystem::path filePath = filename;
	if (not std::filesystem::exists(filePath)) {
		throw std::runtime_error("CrossSection::loadData could not open file " + filePath.string());
	}
	std::ifstream infile(filePath);

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
		infile >> e; // kinetic energy of secondary
		energySecondary.push_back(e + secondary_mc2.at(secondaryId)); // we would like this energy to become the total energy of the secondary here
	}	

	// read the following lines Tp, sigma values
	double a; // helper variable to store 
	while (infile.good()) {
		infile >> a;
		if (not infile) 
			break;    // end of file

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

void HadronicInteraction::CrossSection::setPrimaryId(int id) {
	primaryId = id;
}

void HadronicInteraction::CrossSection::setSecondaryId(int id) {
	secondaryId = id;
}

void HadronicInteraction::CrossSection::setTargetMedium(ref_ptr<MediumComposition> target) {
	targetMedium = target;
}

void HadronicInteraction::CrossSection::setWeightEnergyLoss(double weight) {
	weightEnergyLoss = weight;
}

int HadronicInteraction::CrossSection::getPrimaryId() const {
	return primaryId;
}

int HadronicInteraction::CrossSection::getSecondaryId() const {
	return secondaryId;
}

double HadronicInteraction::CrossSection::getWeightEnergyLoss() const {
	return weightEnergyLoss;
}

ref_ptr<MediumComposition> HadronicInteraction::CrossSection::getTargetMedium() const {
	return targetMedium;
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

HadronicInteraction::HadronicInteraction(ref_ptr<Density> density, bool havePhotons, bool haveElectrons, bool haveNeutrinos, bool haveNucleons, bool haveAntiNucleons, bool catastrophic, ref_ptr<SamplerEvents> sampling, double limit) {
	setInteractionTag("HI");
	setDescription("HadronicInteraction");
	setDensity(density);
	setCatastrophic(catastrophic);
	setHavePhotons(havePhotons);
	setHaveElectrons(haveElectrons);
	setHaveNeutrinos(haveNeutrinos);
	setHaveNucleons(haveNucleons);
	setHaveAntiNucleons(haveAntiNucleons);
	setLimit(limit);
	setSampler(sampling);
	initData();
}

void HadronicInteraction::loadData(std::string file, HadronicInteraction::CrossSection& cs) {
	// load data for a single cross section
	std::ifstream infile(file.c_str());

	if (not infile.good()) 
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
			cs.energySecondary.push_back(e + cs.secondary_mc2.at(cs.secondaryId));
		} catch (std::out_of_range& e) {
			KISS_LOG_ERROR << "HadronicInteraction::loadData: secondaryId " << cs.secondaryId << " not found in secondary_mc2 map." << std::endl;
			throw std::runtime_error("HadronicInteraction: secondaryId not found in secondary_mc2 map");
		}
	}	

	double a;
	while(infile.good()) {
		infile >> a;
		if(not infile) 
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

void HadronicInteraction::initData() {

	ref_ptr<MediumComposition> target = density->getTargetMedium();
	std::string targetLabel = target->getLabel();

	std::cout << "HadronicInteraction: loading cross sections for target medium " << targetLabel << std::endl;
	if (listOfPrimaries.empty()) {
		listOfPrimaries = availableNuclei;
	}

	for (const auto& primaryId : listOfPrimaries) {

		if (haveNeutrinos) {
			std::vector<int> neutrinos = {-12, 12, -14, 14};
			for (const auto& neutrino : neutrinos) {
				CrossSection cs(primaryId, neutrino, target, 1.0);
				cs.loadData();
				addCrossSection(cs);
			}
		}

		if (haveElectrons) {
			std::vector<int> electrons = {-11, 11};
			for (const auto& electron : electrons) {
				CrossSection cs(primaryId, electron, target, 1.0);
				cs.loadData();
				addCrossSection(cs);
			}
		}

		if (havePhotons) {
			CrossSection cs(primaryId, 22, target, 1.0);
			cs.loadData();
			addCrossSection(cs);
		}
		
		if (haveNucleons) {
			std::vector<int> nucleons = {1000010010, 1000010010};
			for (const auto& nucleon : nucleons) {
				CrossSection cs(primaryId, nucleon, target, 1.0);
				cs.loadData();
				addCrossSection(cs);
			}
		}

		if (haveAntiNucleons) {
			std::vector<int> antiNucleons = {-1000010010, 1000010010};
			for (const auto& antiNucleon : antiNucleons) {
				CrossSection cs(primaryId, antiNucleon, target, 1.0);
				cs.loadData();
				addCrossSection(cs);
			}
		}

	}

	std::cout << "HadronicInteraction: loaded " << crossSectionList.size() << " primary IDs with their cross sections." << std::endl;

	// // Initialize data for hadronic interactions
	// for (const auto& [id, crossSections] : crossSectionList) {
	// 	for (const auto& cs : crossSections) {
	// 		// Perform any necessary initialization for each cross section
	// 	}
	// }   
}

// void HadronicInteraction::initData(std::string configFile) {
// 	std::filesystem::path configPath = configFile;
// 	if (not configPath.is_absolute()) {
// 		configPath = std::filesystem::absolute(getDataPath("hadronic_interaction/" + configFile));
// 	}

// 	if (not std::filesystem::exists(configPath)) {
// 		throw std::runtime_error("HadronicInteraction: could not open file " + configPath.string());
// 	}

// 	std::ifstream infile(configPath);

// 	while (infile.good()) {
// 		if (infile.peek() != '#') {
// 			HadronicInteraction::CrossSection cs;
// 			int primId;
// 			int secId;
// 			double w;
// 			bool 
			

// 			infile >> cs.secondaryId; 
// 			infile >> cs.addSecondary;
// 			infile >> cs.weightEnergyLoss;
// 			infile >> cs.primaryId;

// 			std::string path; 
// 			infile >> path;

// 			if (not infile) 
// 				break;    // end of file

// 			// decide if local or global paths are used 
// 			std::filesystem::path loadFile = path;
// 			if (not loadFile.is_absolute()) {
// 				loadFile = std::filesystem::absolute(getDataPath("hadronic_interaction/" + path));
// 			}

// 			loadData(loadFile.string(), cs);
// 			crossSectionList.try_emplace(cs.primaryId).first->second.push_back(cs);
// 		}

// 		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
// 	}

// 	infile.close();
// }

void HadronicInteraction::setListOfPrimaries(std::vector<int> primaries) {
	listOfPrimaries = primaries;
}

void HadronicInteraction::addCrossSection(HadronicInteraction::CrossSection cs) { 
	// problem there : I would need the ID corresponding to the cs to add it at the right place to CSList ; maybe I can add the 
	crossSectionList.try_emplace(cs.primaryId).first->second.push_back(cs);
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

void HadronicInteraction::setHaveNucleons(bool b) {
	haveNucleons = b;

	// if this option is true, catastrophic loss are activated
	if (haveNucleons) 
		setCatastrophic(true);
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

void HadronicInteraction::setCatastrophic(bool b) {
	catastrophicLoss = b;
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

bool HadronicInteraction::isCatastrophic() const {
	return catastrophicLoss;
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

bool HadronicInteraction::getHaveNucleons() const {
	return haveNucleons;
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

bool HadronicInteraction::isPrimaryImplemented(const int& id) const {
	return crossSectionList.find(id) != crossSectionList.end();
}

void HadronicInteraction::process(Candidate* candidate) const {
	int id = candidate->current.getId();
	if (not isPrimaryImplemented(id)) {
		#ifdef DEBUG
			KISS_LOG_ERROR << "HadronicInteraction:process \tid " << id << " is not implemented in the cross section list." << std::endl;
		#endif
		return;
	}

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
	
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double Eprimary = candidate->current.getEnergy() * (1 + z);
	double K = Eprimary - restEnergiesNuclei.at(id); // kinetic energy of the primary in [J] ; Tp is the kinetic energy of the projectile
	double sigmaTot = CrossSection::totalInelasticCrossSection(K, id, At) / (milli * barn);
	double totalEnergyLoss = 0;
	double sampledEnergy = 0;

	std::vector<double> sigmaInterp;
	std::vector<double> dSigma{0, 0, 0};

	Random& random = Random::instance();

// std::cout << "entering the loop" << std::endl;

	// loop over all secondary species
	for (size_t i = 0; i < crossSectionList.at(id + At - 1).size(); i++) {
		sigmaInterp.clear();

		CrossSection cs = crossSectionList.at(id + At - 1)[i];

		// std::cout << "i = " << i << " \t" << cs.primaryId << " \t" << cs.secondaryId << std::endl;

		int index = closestIndex(K / GeV, cs.energyPrimary);
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
		double sigmaRatio = dSigma / sigmaTot;
		int nSecondaries = static_cast<int>(std::floor(sigmaRatio));
		nSecondaries = sigmaRatio - nSecondaries;
		if (random.rand() > nSecondaries) {
			nSecondaries++;
		}

	// std::cout << "nSecondaries = " << nSecondaries << " for id = " << id << ", At = " << At << ", and secondaryId = " << cs.secondaryId << "\t\t sigmaRatio = " << sigmaRatio << std::endl;
		

		// sample secondaries
		for (int j = 0; j < nSecondaries; j++) {
			 // draw a random bin from a (unnormalized) cumulative distribution function, without leading zero
			int bin = random.randBin(sigmaInterp);
			if (bin > 0) {
				sampledEnergy = cs.energySecondary[bin - 1] + random.rand() * (cs.energySecondary[bin] - cs.energySecondary[bin - 1]);
			} else {
				sampledEnergy = 1 + random.rand() * (cs.energySecondary[bin] - 1);
			}
			totalEnergyLoss += sampledEnergy * cs.weightEnergyLoss;

			if (havePhotons and cs.secondaryId == 22) {
				candidate->addSecondary(cs.secondaryId, sampledEnergy * GeV, 1., interactionTag);
			} else if (haveElectrons and abs(cs.secondaryId) == 11) {
				candidate->addSecondary(cs.secondaryId, sampledEnergy * GeV, 1., interactionTag);
			} else if (haveNeutrinos and (abs(cs.secondaryId) == 12 or abs(cs.secondaryId) == 14 or abs(cs.secondaryId) == 16)) {
				candidate->addSecondary(cs.secondaryId, sampledEnergy * GeV, 1., interactionTag);
			} else if (haveAntiNucleons and (cs.secondaryId == -1000010010 or cs.secondaryId == -1000000010 or cs.secondaryId == -2212 or cs.secondaryId == -2112)) {
				candidate->addSecondary(cs.secondaryId, sampledEnergy * GeV, 1., interactionTag);
			} else {
				candidate->addSecondary(cs.secondaryId, sampledEnergy * GeV, 1., interactionTag);
			}

		}
	}

	// apply energy loss
	candidate->current.setEnergy(K - totalEnergyLoss * GeV);
	// if (catastrophicLoss)
	// 	candidate->setActive(false);
		
}



} // namespace crpropa
