#include "crpropa/ParticleID.h"


namespace crpropa {

int nucleusId(int a, int z) {
	if (z < 0)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with Z < 0, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < 1)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with A < 1, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	if (a < z)
		throw std::runtime_error(
				"crpropa::Nucleus: no nucleus with A < Z, A=" + kiss::str(a) + " Z="
						+ kiss::str(z));
	return 1000000000 + z * 10000 + a * 10;
}

int chargeNumber(int id) {
	return HepPID::Z(id);
}

int massNumber(int id) {
	if (id == 2112)
		return 1;
	return HepPID::A(id);
}

bool isNucleus(int id) {
	if (id == 2112)
		return true; // consider neutron as nucleus
	return HepPID::isNucleus(id);
}

bool isMeson(int id) {
	return HepPID::isMeson(id);
}

bool isBaryon(int id) {
	return HepPID::isBaryon(id);
}

bool isHadron(int id) {
	return HepPID::isHadron(id);
}

bool isLepton(int id) {
	return HepPID::isLepton(id);
}

bool hasUp(int id) {
	return HepPID::hasUp(id);
}

bool hasDown(int id) {
	return HepPID::hasDown(id);
}

bool hasCharm(int id) {
	return HepPID::hasCharm(id);
}

bool hasStrange(int id) {
	return HepPID::hasStrange(id);
}

bool hasTop(int id) {
	return HepPID::hasTop(id);
}

bool hasBottom(int id) {
	return HepPID::hasBottom(id);
}

std::string convertIdToName(int id) {
	// handle a few extra cases that HepPID doesn't like
	if (id == 1000000010) // neutron
		id = 2112;
	if (id == -1000000010) // anti-neutron
		id = -2112;
	if (id == -1000010010) // anti-proton
		id = -2212;
	return HepPID::particleName(id);
}

}
