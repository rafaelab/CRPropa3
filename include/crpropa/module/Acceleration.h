#ifndef CRPROPA_ACCELERATION_H
#define CRPROPA_ACCELERATION_H

#include <cmath>
#include <string>

#include "crpropa/Candidate.h"
#include "crpropa/Common.h"
#include "crpropa/Geometry.h"
#include "crpropa/Module.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"



namespace crpropa {

/** @addtogroup Acceleration
 *  @{
 */

/** @class StepLengthModifier
 * @brief Modifies the steplength of an acceleration module.
 */
class StepLengthModifier : public Referenced {
	public:
		/** Returns an update of the step length
		 * @param steplength 	modifies step length, e.g., based on scattering model
		 * @param candidate 	additional candidate properties are usually included in the calculation of the updated step length.
		 */
		virtual double modify(double stepLength, Candidate* candidate) = 0;
};


/**
 * @class AbstractAccelerationModule
 * @brief Core functionallity for acceleration by scattering with scattering centers moving in a velocity field.
 * @details The velocity field is implicity implemented in the derived classes for performance reasons. 
 *  Models for the dependence of the step length of the scatter process are set via modifiers.
 */
class AbstractAccelerationModule : public Module {
	double stepLength;
	std::vector<ref_ptr<StepLengthModifier>> modifiers;

	public:
	/// The parent's constructor need to be called on initialization!
	AbstractAccelerationModule(double step = 1. * parsec);

	// add a step length modifier to the model
	void add(StepLengthModifier* modifier);

	// update the candidate
	void process(Candidate* candidate) const;

	// Returns the velocity vector of the scatter centers in the rest frame of the candidate.
	// This needs to be implemented in inheriting classes.
	virtual Vector3d scatteringCenterVelocity(Candidate* candidate) const = 0;

	// Scatter the candidate with a center with given scatter center velocity into a random direction. 
	// Assumes that the candidate is ultra-relativistic (m = 0).
	void scatter(Candidate* candidate, const Vector3d& scatteringCenterVelocity) const;
};


/** @class SecondOrderFermi
 * @brief  Implements scattering with centers moving in isotropic directions.
 *   All scatter centers have the same velocity.
 */
class SecondOrderFermi : public AbstractAccelerationModule {
	private:
		double scatterVelocity;
		std::vector<double> angle;
		std::vector<double> angleCDF;

	public:
		/** Constructor
		@param scatterVelocity			velocity of scattering centers
		@param stepLength				average mean free path
		@param sizeOfPitchAngleTable	number of precalculated pitch angles
		*/
		SecondOrderFermi(double scatterVelocity = 0.1 * c_light, double stepLength = 1. * parsec, unsigned int sizeOfPitchAngleTable = 10000);
		virtual Vector3d scatteringCenterVelocity(Candidate* candidate) const;
};


/** @class DirectedFlowScattering
 * @brief Scattering in a directed flow of scatter centers.
 * @details Two of these region with different velocities can be used to create first order Fermi scenario.
 * Thanks to Aritra Ghosh, Groningen University, for first work in 2017 on the shock acceleration in CRPropa, leading to this module.
 */
class DirectedFlowScattering : public AbstractAccelerationModule {
	private:
		Vector3d scatterVelocity;

	public:
		/** Constructor
		 * @param scatteringCenterVelocity	velocity of scattering centers
		 * @param stepLength				average mean free path
		 */
		DirectedFlowScattering(Vector3d scatteringCenterVelocity, double stepLength = 1. * parsec);
		void setScatteringVelocity(Vector3d velocity);
		virtual Vector3d scatteringCenterVelocity(Candidate* candidate) const;
};


/** @class DirectedFlowOfScatterCenters
 * @brief In a directed flow, the step length depend on the direction of the particles as head-on collisions are more likely than tail-on collision. 
 * Therefore, propagation against the flow is harder.
 */
class DirectedFlowOfScatterCenters : public StepLengthModifier {
	private:
		Vector3d scatterVelocity;

	public:
		/** Constructor
		 * @param scatteringCenterVelocity	velocity of scattering centers
		 */
		DirectedFlowOfScatterCenters(const Vector3d& scatteringCenterVelocity);
		void setScatteringVelocity(Vector3d velocity);
		double modify(double steplength, Candidate* candidate);
};


/** @class QuasiLinearTheory
 * @brief Scales the steplength according to quasi linear theory.
 * @details Following quasi-linear theory [Schlickeiser1989], the mean free path \f$\lambda\f$ of a particle with energy \f$E\f$ and charge \f$Z\f$ in a field with turbulence spectrum \f$\frac{k}{k_{\min}}^{-q}\f$ is
 *   \f$ \lambda = {\left(\frac{B}{\delta B}\right)}^2 {\left(R_G\; k_{\min}\right)}^{1-q} R_G \equiv \lambda_0 {\left( \frac{E}{1 *   EeV}\frac{1}{Z} \right)}^{2-q} \f$,
 *  where \f$R_G = \frac{E}{B Z}\f$ is the gyro-radius of the particles.
 * This class implements the rigidity dependent scaling factor used to modify the base step length.
 * 
 * \b [Schlickeiser1989]
 *  "Cosmic-Ray Transport and Acceleration. II. Cosmic Rays in Moving Cold Media with Application to Diffusive Shock Wave Acceleration",
 *  R. Schlickeiser.
 *  The Astrophysical Journal 336 (1989) 264. 
 *  doi:10.1086/167010.
 */
class QuasiLinearTheory : public StepLengthModifier {
	private:
		double referenceEnergy;
		double turbulenceIndex;
		double minimumRigidity;

	public:
		/** Constructor
		 * @param referenceEnergy	reference energy - break of power spectrum
		 * @param turbulenceIndex	power law index of the isotropic magnetic turbulence power spectrum; default is set to Kolmogorov turbulence
		 * @param minimumRigidity	minimal rigidity
		 */
		QuasiLinearTheory(double referenceEnergy = 1. * EeV, double turbulenceIndex = 5. / 3, double minimumRigidity = 0);
		void setReferenceEnergy(double energy);
		void setTurbulenceIndex(double index);
		void setMinimumRigidity(double rigidity);
		double modify(double steplength, Candidate* candidate);
};


/** @class ParticleSplitting
 * @brief  Implements particle splitting, i.e. inverse thinning, to speed up the simulation.
 * @details After crossing a surface a given number of times, the particle is split to N particles with weight 1/N. 
 * This eases performance constraints in acceleration simulations due to the power law nature of many acceleration mechanisms.
 * This feature is thanks to Matthew Weiss, Penn State University for the first work on this in 2017.
 */
class ParticleSplitting : public Module {
	private:
		ref_ptr<Surface> surface;
		std::string counterid;
		int numberSplits;
		int crossingThreshold;
		double minWeight;

	public:
		/** Constructor
		@param surface             The surface to monitor
		@param crossingThreshold   Number of crossings after which a particle is split
		@param numberSplits        Number of particles the candidate is split into
		@param minWeight           Minimum weight to consider. Particles with a lower weight are not split again.
		@param counterid           An unique string to identify the particle property used for counting; useful if multiple splitting modules are present.
		*/
		ParticleSplitting(Surface* surface, int crossingThreshold = 50, int numberSplits = 5, double minWeight = 0.01, std::string counterid = "ParticleSplittingCounter");

		void setSurface(ref_ptr<Surface> surface);
		void setCounterId(std::string counterid);
		void setCrossingThreshold(int crossingThreshold);
		void setNumberSplits(int numberSplits);
		void setMinimumWeight(double minWeight);

		// update the candidate
		void process(Candidate* candidate) const;
};


/**  @} */ // end of group Acceleration

} // namespace crpropa

#endif
