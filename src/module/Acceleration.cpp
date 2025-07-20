#include "crpropa/module/Acceleration.h"


namespace crpropa {


AbstractAccelerationModule::AbstractAccelerationModule(double step) : Module(), stepLength(step) {
}


void AbstractAccelerationModule::add(StepLengthModifier* modifier) {
	modifiers.push_back(modifier);
}


void AbstractAccelerationModule::scatter(Candidate* candidate, const Vector3d& v) const {
	// particle momentum in lab frame
	const double E = candidate->current.getEnergy();
	const Vector3d p = candidate->current.getMomentum();

	// transform to rest frame of scatter center (p: prime)
	const double beta = v.getR() / c_light;
	const double gamma = 1. / sqrt(1 - beta * beta);
	const double Ep = gamma * (E - v.dot(p));
	const Vector3d pp = (p - v * E / (c_light * c_light)) * gamma;

	// scatter into random direction
	const Vector3d pp_new = Random::instance().randVector() * pp.getR();

	// transform back
	const double E_new = gamma * (Ep + v.dot(pp_new));
	const Vector3d p_new = (pp_new + v * Ep / (c_light * c_light)) * gamma;

	// update candidate properties
	candidate->current.setEnergy(E_new);
	candidate->current.setDirection(p_new / p_new.getR());
}


void AbstractAccelerationModule::process(Candidate* candidate) const {
	double currentStepLength = stepLength;
	for (auto m : modifiers) {
		currentStepLength = m->modify(currentStepLength, candidate);
	}

	double step = candidate->getCurrentStep();
	while (step > 0) {
		double randDistance = -1. * log(Random::instance().rand()) * currentStepLength;

		if (step < randDistance) {
			candidate->limitNextStep(0.1 * currentStepLength);
			return;
		}
		scatter(candidate, scatteringCenterVelocity(candidate));
		step -= randDistance;
	}
}


SecondOrderFermi::SecondOrderFermi(double scatterVelocity, double stepLength, unsigned int sizeOfPitchangleTable) : AbstractAccelerationModule(stepLength), scatterVelocity(scatterVelocity) {
	setDescription("SecondOrderFermi Acceleration");
	angle.resize(sizeOfPitchangleTable);
	angleCDF.resize(sizeOfPitchangleTable);

	// have a discretized table of beamed pitch angles
	for (size_t i = 0; i < sizeOfPitchangleTable; i++) {
		angle[i] = i * M_PI / (sizeOfPitchangleTable - 1);
		angleCDF[i] = (angle[i] + scatterVelocity / c_light * sin(angle[i])) / M_PI;
	}
}


Vector3d SecondOrderFermi::scatteringCenterVelocity(Candidate* candidate) const {
	size_t idx = closestIndex(Random::instance().rand(), angleCDF);
	Vector3d rv = Random::instance().randVector();
	Vector3d rotationAxis = candidate->current.getDirection().cross(rv);

	rv = candidate->current.getDirection().getRotated(rotationAxis, M_PI - angle[idx]);
	return rv * scatterVelocity;
}


DirectedFlowOfScatterCenters::DirectedFlowOfScatterCenters(const Vector3d& scatteringCenterVelocity) {
	setScatteringVelocity(scatteringCenterVelocity);
}

void DirectedFlowOfScatterCenters::setScatteringVelocity(Vector3d velocity) {
	scatterVelocity = velocity;
}

double DirectedFlowOfScatterCenters::modify(double step, Candidate* candidate) {
	double directionModifier = (-1. * scatterVelocity.dot(candidate->current.getDirection()) + c_light) / c_light;
	return step / directionModifier;
}


DirectedFlowScattering::DirectedFlowScattering(Vector3d scatteringCenterVelocity, double step) : AbstractAccelerationModule(step) {
	setScatteringVelocity(scatteringCenterVelocity);

	// In a directed field of scatter centers, the probability to encounter a scatter center depends on the direction of the candidate.
	StepLengthModifier* mod = new DirectedFlowOfScatterCenters(scatterVelocity);
	this->add(mod);
}

void DirectedFlowScattering::setScatteringVelocity(Vector3d velocity) {
	scatterVelocity = velocity;
}

Vector3d DirectedFlowScattering::scatteringCenterVelocity(Candidate* candidate) const { 
	// does not depend on candidate here.
	return scatterVelocity;
}


QuasiLinearTheory::QuasiLinearTheory(double referenceEnergy, double turbulenceIndex, double minimumRigidity) {
	setReferenceEnergy(referenceEnergy);
	setTurbulenceIndex(turbulenceIndex);
	setMinimumRigidity(minimumRigidity);
}

void QuasiLinearTheory::setReferenceEnergy(double energy) {
	referenceEnergy = energy;
}

void QuasiLinearTheory::setTurbulenceIndex(double index) {
	turbulenceIndex = index;
}

void QuasiLinearTheory::setMinimumRigidity(double rigidity) {
	minimumRigidity = rigidity;
}

double QuasiLinearTheory::modify(double steplength, Candidate* candidate) {
	if (candidate->current.getRigidity() < minimumRigidity) {
		return steplength * std::pow(minimumRigidity / (referenceEnergy / eV), 2. - turbulenceIndex);
	} else 	{
		return steplength * std::pow(candidate->current.getRigidity() / (referenceEnergy / eV), 2. - turbulenceIndex);
	}
}


ParticleSplitting::ParticleSplitting(Surface* surface, int crossingThreshold, int numberSplits, double minWeight, std::string counterId) {
	setSurface(surface);
	setCrossingThreshold(crossingThreshold);
	setNumberSplits(numberSplits);
	setMinimumWeight(minWeight);
	setCounterId(counterId);
}

void ParticleSplitting::setSurface(ref_ptr<Surface> surf) {
	surface = surf;
}

void ParticleSplitting::setCounterId(std::string id) {
	counterid = id;
}

void ParticleSplitting::setCrossingThreshold(int threshold) {
	crossingThreshold = threshold;
}

void ParticleSplitting::setNumberSplits(int splits) {
	numberSplits = splits;
}

void ParticleSplitting::setMinimumWeight(double weight) {
	minWeight = weight;
}

void ParticleSplitting::process(Candidate* candidate) const {
	const double currentDistance = surface->distance(candidate->current.getPosition());
	const double previousDistance = surface->distance(candidate->previous.getPosition());

	// candidate remains on the same side
	if (currentDistance * previousDistance > 0)
		return;

	if (candidate->getWeight() < minWeight)
		return;

	int nCrossings = 1;
	if (candidate->hasProperty(counterid))
		nCrossings = candidate->getProperty(counterid).toInt32() + 1;
	candidate->setProperty(counterid, nCrossings);

	if (nCrossings % crossingThreshold != 0)
		return;

	candidate->updateWeight(1. / numberSplits);

	for (size_t i = 1; i < numberSplits; i++) {
		// no recursive split as the weights of the secondaries created before the split are not affected
		ref_ptr<Candidate> newCandidate = candidate->clone(false);
		newCandidate->parent = candidate;
		uint64_t snr = Candidate::getNextSerialNumber();
		Candidate::setNextSerialNumber(snr + 1);
		newCandidate->setSerialNumber(snr);
		candidate->addSecondary(newCandidate);
	}
};

} // namespace crpropa
