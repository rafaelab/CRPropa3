#include "crpropa/module/HadronicInteraction.h"



namespace crpropa {


HadronicInteraction::HadronicInteraction(ref_ptr<Density> dens) : Module(), density(dens), limit(0.1) {
	setDescription("HI");
}

HadronicInteraction::HadronicInteraction(ref_ptr<Density> dens, std::string configFile, double limit) : Module(), density(dens), limit(limit) {
	setDescription("HI");
	if (configFile != "")
		initData(configFile);
}

void HadronicInteraction::addCrossSection(HadronicInteraction::CrossSection cs) { 
	// problem there : I would need the ID corresponding to the cs to add it at the right place to CSList ; maybe I can add the 
	crossSectionList.at(cs.CR_Id).push_back(cs);
}

HadronicInteraction::CrossSection::CrossSection() { 
}

HadronicInteraction::CrossSection::CrossSection(std::string file) { 
	loadData(file);
}

HadronicInteraction::CrossSection::CrossSection(std::string file, int id, double loss, bool add, int CR_Id) : nucleusId(id), weightEnergyLoss(loss), addSecondary(add), CR_Id(CR_Id) {
	loadData(file);
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

double HadronicInteraction::getLimit() const {
	return limit;
}

bool HadronicInteraction::getDeactivatePrimary() const {
	return deactivatePrimary;
}

std::string HadronicInteraction::getInteractionTag() const {
	return interactionTag;
}

double HadronicInteraction::totalInelasticCrossSection(double Tp, int id, int A_t) const { 
	// this is where I should add the modified formula for nucleus/nucleus interactions ; they only have p-p interactions

	int A_p = massNumber(id); // mass number of the projectile
	Tp = Tp / A_p; // in the formula, Tp is the projectile kinetic energy per nucleon
	double Tp_th = 0.2797 * GeV; // kinematic threshold // could replace this value by the real one, from the paper of Kafexhiu
	double r = Tp / Tp_th;
	double logr = std::log(r);

	double sigma_pp = (30.7 - 0.96 * logr + 0.18 * logr * logr) * pow_integer<3>(1 - pow(1 / r, 1.9)) * milli * barn; // warning, in this formula, Tp has to be the kinetic energy

	if ((id == 1000010010) and (A_t == 1))
		return std::max(0., sigma_pp); // clip with 0 

	double T_p0 = 1e3 * GeV; // in eV, fixed in the paper
	double r_p0 = T_p0 / Tp_th;
	double logr_p0 = std::log(r_p0);

	double sigma_Tp0 = (30.7 - 0.96 * logr_p0 + 0.18 * logr_p0 * logr_p0) * pow_integer<3>(1 - pow(1/r_p0, 1.9)) * milli * barn;
	
	double G = 1 + std::log(std::max(1.0, sigma_pp / sigma_Tp0)); // np.maximum : element-wise comparison
	double r0 = 1.36 * femto;
	double sigma_r0 = 3.14159265 * r0 * r0; // around 58.1 mb

	double beta0 = 2.247 - 0.915 * (1 + pow(A_t, -1/3));

	if (A_p != 1){
		beta0 = 1.581 - 0.876 * (pow(A_p, -1/3) + pow(A_t, -1/3)) ;
	}

	double sigma_r = sigma_r0 * pow_integer<2>((pow(A_p, 1/3) + pow(A_t, 1 / 3) - beta0 * ( pow(A_p, - 1 / 3) + pow(A_t, -1 / 3))));
	
	return sigma_r * G;
}

double HadronicInteraction::getDensityAtPosition(Vector3d& pos) const {
	return density->getDensity(pos);
}

void HadronicInteraction::CrossSection::loadData(std::string file) { // function of the class crosssection
	std::ifstream infile(file.c_str());

	if (! infile.good()) 
		throw std::runtime_error("CrossSection::loadData could not open file " + file);

	// clear old data 
	sigma.clear();
	Tp.clear();
	eps.clear();

	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	// read secondary energy as the first row
	double e; 
	infile >> e; // ignore the first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> e;
		eps.push_back(e + secondary_mc2.at(nucleusId)); // we would like this energy to become the total energy of the secondary here
	}	

	// read the following lines Tp, sigma values
	double a; // helper variable to store 
	while (infile.good()) {
		infile >> a;
		if (! infile) 
			break;	// end of file
		
		Tp.push_back(a);
	
		std::vector<double> cdf;
		for (int i = 0; i < eps.size(); i++){
			infile >> a;
			cdf.push_back(a);
		}
		sigma.push_back(cdf);
	}


	// finaly close the file
	infile.close();
	
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
		cs.eps.push_back(e + cs.secondary_mc2.at(cs.nucleusId));
	}	

	// read the following lines Tp, sigma values
	double a;
	while(infile.good()) {
		infile >> a;
		if(! infile) 
			break;
		
		cs.Tp.push_back(a);
	
		std::vector<double> cdf;
		for (int i = 0; i < cs.eps.size(); i++){
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
		
		if(!infile) 
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

void HadronicInteraction::allowSecondaryId(int id, bool allow) {
	bool foundId = false;
	for(int i = 0; i < crossSectionList[id].size(); i++) {
		if(crossSectionList[id][i].nucleusId == id){
			foundId = true;
			crossSectionList[id][i].addSecondary = allow;
		}
	}

	if (! foundId)
		KISS_LOG_WARNING << "HadronicInteraction::allowSecondaryId \tid " << id << " was not found in the list of cross section.\n";  
		throw;
}

void HadronicInteraction::process(Candidate* cand) const {
	int id = cand -> current.getId();

	// allow Iron, helium, C, Al and p
	if((id != 1000010010) and (id != 1000020040) and (id != 1000060120) and (id != 1000130260) and (id != 1000260560)) 
		return;

	// get the kinetic energy rather than total energy
	double Tp = cand->current.getEnergy(); 
	if (Tp < CR_min_energy.at(id))
		return;
	Tp = Tp - CR_rest_energy.at(id);

	double step = cand->getCurrentStep();
	Vector3d pos = cand->current.getPosition();
	Vector3d dir = cand->current.getDirection();

	// density of HI + HII + 2*H2
	double n0 = density->getNucleonDensity(pos); 

	//this is not a proba ; kind of strange there, maybe should check the effect of the limitNextStep (maybe this solves the issue)
	double prob_H = totalInelasticCrossSection(Tp, id, 1) * n0 * step; 

	int A_t = 1;
	Random& random = Random::instance();

	if (id == 1000010010 or id == 1000020040) {
		// He density ; this function doesn’t seem to exist yet
		double n_He = n0 * 0.1;
		// double n_He = density->getHeDensity(pos); 
		double prob_He = totalInelasticCrossSection(Tp, id, 4) * n_He * step; // here it should also be a 4
		double prob = prob_H + prob_He; // here Tp is the kinetic energy ; for targets : protons and Helium (ie A_t = 1 and 4) ; // 

		if (random.rand() < prob) {
			if (random.rand() * prob < prob_He) {
				A_t = 4;
			}
			performInteraction(cand, A_t);
		}
		else {
			cand->limitNextStep(limit * step / prob);
			double n1_H = density->getNucleonDensity(pos + step * dir);
			// double n1_He = density->getHeDensity(pos + step * dir); // same here for the He density
			double n1_He = 0.1 * n1_H;
			//double dNdX = std::max(fabs(n1_H - n0) / n0, fabs(n1_He - n_He) / n_He) / step;  // change the size of the step size to account for change of density
			double dNdX = fabs((n1_H - n0) * n_He + (n1_He - n_He) * n0) / (step * n0 * n_He);
			cand->limitNextStep(limit / dNdX);
		}
	} else {
		if (random.rand() < prob_H) {
			performInteraction(cand, A_t);
		} else {
			cand->limitNextStep(limit * step / prob_H);
			double n1_H = density->getNucleonDensity(pos + step * dir);
			double dNdX = fabs(n0 - n1_H) / n0 / step;  // change the size of the step size to account for change of density
			cand->limitNextStep(limit / dNdX);
		}
	}
}

void HadronicInteraction::performInteraction(Candidate* cand, int A_t) const {
	
	bool activate = true;
	int id = cand -> current.getId();
	double Tp = cand -> current.getEnergy(); // once again it’s total
	double totalEnergyLoss = 0; // in GeV
	double sigmaTot = totalInelasticCrossSection(Tp - CR_rest_energy.at(id), id, A_t) / (milli * barn); // the input of this function is kinetic energy of the primary
	double sampledEnergy = 0;

	std::vector<double> sigma_interp;
	std::vector<double> dsigma{0, 0, 0};

	Random &random = Random::instance();

	// loop over all secondary species
	for(size_t i = 0; i < crossSectionList.at(id + A_t - 1).size(); i++) {
		//std::cout << "begining 1st for" << std::endl;
		// crosssection
		sigma_interp.clear();

		CrossSection cs = crossSectionList.at(id + A_t - 1)[i];
		int index = closestIndex(Tp / GeV, cs.Tp); // then there is no interpolation there ; here Tp is total energy
		int l = 0;
		
		if (index == 0) {
			l = 1;
		} else if (index == cs.Tp.size() - 1) {
			l = -1;
		}
		std::vector<double> subvector_Ep{cs.Tp[index + l - 1], cs.Tp[index + l], cs.Tp[index + l + 1]};

		for (size_t k = 0; k < cs.sigma[index].size(); k++){
			dsigma[0] = cs.sigma[index + l - 1][k];
			dsigma[1] = cs.sigma[index + l][k];
			dsigma[2] = cs.sigma[index + l + 1][k];
			sigma_interp.push_back(interpolate(Tp / GeV, subvector_Ep, dsigma));
		}
		
		// integrated differential crosssection #no interpolation # get the total integrated diff cross section
		double dsigma = sigma_interp[cs.eps.size() - 1]; 

		 // number of produced secondaries (for a type of secondary)
		 // draw a random bin from a (unnormalized) cumulative distribution function, 
		// int nSecondaries = randomRound(dsigma / sigmaTot);
		double sigmaRatio = dsigma / sigmaTot;
		int nSecondaries = static_cast<int>(std::floor(sigmaRatio));
		nSecondaries = sigmaRatio - nSecondaries;
		if (random.rand() > nSecondaries) {
			nSecondaries++;
		}
		

		// sample secondaries
		for(size_t j = 0; j < nSecondaries; j++) {
			 // draw a random bin from a (unnormalized) cumulative distribution function, without leading zero
			int bin = random.randBin(sigma_interp);
			if (bin > 0) {
				 // from that, in the tables, E_sec should be total energy ; or I should add the mass energy ; also I should not take the middle of the bins
				sampledEnergy = cs.eps[bin - 1] + random.rand() * (cs.eps[bin] - cs.eps[bin - 1]);
			} else {
				// first possible value for us is 1 GeV
				sampledEnergy = 1 + random.rand() * (cs.eps[bin] - 1); 
			}
			totalEnergyLoss += sampledEnergy * cs.weightEnergyLoss;

			// add secondary if allowed
			if(cs.addSecondary) {
				cand->addSecondary(cs.nucleusId, sampledEnergy * GeV, 1., interactionTag);
				if ((cs.nucleusId == id) and ((id == 1000000010) or (id == 1000010010))){
					activate = false;
				}
			}
		}
	}
	//std::cout << "end for perform interaction" << std::endl;

	// apply energy loss
	cand->current.setEnergy(Tp - totalEnergyLoss * GeV);
	cand->setActive(activate);
	if (deactivatePrimary)
		cand->setActive(false);
}



} // namespace crpropa 
