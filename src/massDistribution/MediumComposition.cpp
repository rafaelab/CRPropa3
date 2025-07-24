#include "crpropa/massDistribution/MediumComposition.h"


namespace crpropa {


MediumCompositionElementary::MediumCompositionElementary(std::string label) {
}

MediumCompositionElementary::MediumCompositionElementary(int id, std::string label) {
	setParticleId(id);
}

void MediumCompositionElementary::setParticleId(int id) {
	if (isNucleus(id)) {
		KISS_LOG_WARNING << "MediumCompositionElementary is more appropriate for more fundamental particles. The Id provided is for a nucleus. Try using MediumCompositionAtomic. I will proceed anyway." << std::endl;
	}
	particleId = id;
}

int MediumCompositionElementary::getParticleId() const {
	return particleId;
}

bool MediumCompositionElementary::isAdmixed() const {
	return false;
}

bool MediumCompositionElementary::isMolecular() const {
	return false;
}

bool MediumCompositionElementary::isIonized() const {
	return false;
}

bool MediumCompositionElementary::isNeutral() const {
	return HepPID::charge(particleId) == 0;
}

unsigned int MediumCompositionElementary::getNumberOfNucleons() const {
	if (particleId == 2212 or particleId == 2112) { // p or n
		return 1;
	}
	if (isNucleus(particleId)) {
		return massNumber(particleId);
	}
	return 0; // no nucleons
}

double MediumCompositionElementary::getCompositionWeight() const {
	return 0;
}

std::string MediumCompositionElementary::getDescription() const {
	std::stringstream ss;
	ss << "MediumCompositionElementary with particle ID: " << particleId << std::endl;
	return ss.str();
}


///////////////////////////////////////////////////////////////////////////////

MediumCompositionAtomic::MediumCompositionAtomic(std::string label) {
}

MediumCompositionAtomic::MediumCompositionAtomic(int id, int nElectrons, std::string label) {
	setAtomId(id);
	setNumberOfElectrons(nElectrons);
}

void MediumCompositionAtomic::setAtomId(int id) {
	if (not isNucleus(id)) {
		KISS_LOG_WARNING << "MediumCompositionAtomic works for atomic nuclei only. You provided a particle that is not a nucleus." << std::endl;
		throw std::invalid_argument("MediumCompositionAtomic works for atomic nuclei only. You provided a particle that is not a nucleus.");
	}
	atomId = id;
}

void MediumCompositionAtomic::setNumberOfElectrons(int n) {
	nElectrons = n;
}

int MediumCompositionAtomic::getAtomId() const {
	return atomId;
}

int MediumCompositionAtomic::getNumberOfElectrons() const {
	return nElectrons;
}

bool MediumCompositionAtomic::isAdmixed() const {
	return false;
}

bool MediumCompositionAtomic::isMolecular() const {
	return false;
}

bool MediumCompositionAtomic::isIonized() const {
	return  nElectrons != chargeNumber(atomId);
}

bool MediumCompositionAtomic::isNeutral() const {
	return nElectrons == chargeNumber(atomId);
}

unsigned int MediumCompositionAtomic::getNumberOfNucleons() const {
	if (atomId == 2212 or atomId == 2112) { // p or n
		return 1;
	}
	if (isNucleus(atomId)) {
		return massNumber(atomId);
	}
	return 0; // no nucleons
}

double MediumCompositionAtomic::getCompositionWeight() const {
	return 1.0 * getNumberOfNucleons();
}

std::string MediumCompositionAtomic::getDescription() const {
	std::stringstream ss;
	ss << "MediumCompositionAtomic with atom ID: " << atomId << " and number of electrons: " << nElectrons << std::endl;
	return ss.str();
}


///////////////////////////////////////////////////////////////////////////////

MediumCompositionMolecular::MediumCompositionMolecular(std::string label) {
}

MediumCompositionMolecular::MediumCompositionMolecular(int nElectrons, std::string label) {
	atomIds.clear();
	setNumberOfElectrons(nElectrons);
}

MediumCompositionMolecular::MediumCompositionMolecular(const std::vector<int>& ids, int nElectrons, std::string label) {
	for (const auto& id : ids) {
		add(id);
	}
	setNumberOfElectrons(nElectrons);
}

void MediumCompositionMolecular::add(int id) {
	if (! isNucleus(id)) {
		KISS_LOG_WARNING << "MediumCompositionMolecular works for atomic nuclei only. You provided a particle that is not a nucleus." << std::endl;
		throw std::invalid_argument("MediumCompositionMolecular works for atomic nuclei only. You provided a particle that is not a nucleus.");
	}
	atomIds.push_back(id);
}

void MediumCompositionMolecular::setAtomIds(const std::vector<int>& ids) {
	atomIds = ids;
}

void MediumCompositionMolecular::setNumberOfElectrons(int n) {
	nElectrons = n;
}

int MediumCompositionMolecular::getNumberOfElectrons() const {
	return nElectrons;
}

const std::vector<int>& MediumCompositionMolecular::getAtomIds() const {
	return atomIds;
}

bool MediumCompositionMolecular::isAdmixed() const {
	return false;
}

bool MediumCompositionMolecular::isMolecular() const {
	return true;
}

bool MediumCompositionMolecular::isIonized() const {
	return nElectrons != 0;
}

bool MediumCompositionMolecular::isNeutral() const {
	return nElectrons == 0;
}

unsigned int MediumCompositionMolecular::getNumberOfNucleons() const {
	int nNucleons = 0;
	for (const auto& atom : atomIds)
		nNucleons += massNumber(atom);
	return nNucleons;
}

double MediumCompositionMolecular::getCompositionWeight() const {
	return 1.0 * getNumberOfNucleons();
}

std::string MediumCompositionMolecular::getDescription() const {
	std::stringstream ss;
	ss << "MediumCompositionMolecular with atom IDs: ";
	for (const auto& id : atomIds) {
		ss << id << " ";
	}
	ss << "and number of electrons: " << nElectrons << std::endl;
	return ss.str();
}


///////////////////////////////////////////////////////////////////////////////

MediumCompositionList::MediumCompositionList(std::string label) {
}

MediumCompositionList::MediumCompositionList(const std::vector<ref_ptr<MediumComposition>>& comps, const std::vector<double>& weights, std::string label) {
	if (comps.size() != weights.size()) {
		KISS_LOG_WARNING << "MediumCompositionList: Compositions and weights must have the same size." << std::endl;
		throw std::invalid_argument("MediumCompositionList: Compositions and weights must have the same size.");
	}

	for (size_t i = 0; i < comps.size(); ++i)
		add(comps[i], weights[i]);
}

void MediumCompositionList::add(ref_ptr<MediumComposition> composition, double weight) {
	compositions.push_back(composition);
	weights.push_back(weight);
}

const std::vector<ref_ptr<MediumComposition>>& MediumCompositionList::getCompositions() const {
	return compositions;
}

const std::vector<double>& MediumCompositionList::getWeights() const {
	return weights;
}

bool MediumCompositionList::isAdmixed() const {
	return numberOfComponents() > 1;
}

bool MediumCompositionList::isMolecular() const {
	return false;
}

bool MediumCompositionList::isIonized() const {
	return false;
}

bool MediumCompositionList::isNeutral() const {
	return false;
}

unsigned int MediumCompositionList::numberOfComponents() const {
	return compositions.size();
}

double MediumCompositionList::getCompositionWeight() const {
	double weight = 0.0;
	double sumWeights = accumulate(weights.begin(), weights.end(), 0.);
	for (size_t i = 0; i < numberOfComponents(); ++i) {
		weight += (weights[i] * compositions[i]->getCompositionWeight() / sumWeights);
	}

	return weight;
}

std::string MediumCompositionList::getDescription() const {
	std::stringstream ss;
	ss << "MediumCompositionList with " << compositions.size() << " components." << std::endl;
	for (const auto& comp : compositions)
		ss << comp->getDescription();

	return ss.str();
}


} // namespace crpropa
