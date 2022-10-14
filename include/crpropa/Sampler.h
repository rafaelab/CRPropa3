#ifndef CRPROPA_SAMPLER_H
#define CRPROPA_SAMPLER_H

#include <vector>

#include <crpropa/Common.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>



namespace crpropa {

 /**
 * \addtogroup Sampler
 * @{
 */
 

 /**
 @class Sampler
 @brief Abstract base class to handle all samplers.
 */
class Sampler: public Referenced {
	public:
		virtual double computeWeight(int id, double energy = 0, double energyFraction = 0) const = 0;
};



/**
 @class SamplerUniform
 @brief Throw away a fraction of the particles of a given type.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 Note that the energy fraction is implemented through setters and getters, but it will constantly be modified.
 This is because the structure of computeWeight should remain constant (or variadic arguments used). 
 This will be done in the future.
 */
class SamplerUniform: public Sampler {
	private:
		int particleId;
		double sampling;

	public:
		SamplerUniform(int particleId, double sampling);
		void setSampling(double sampling);
		void setParticleId(int particleId);
		double getSampling() const;
		int getParticleId() const;
		double computeWeight(int id, double energy = 0, double energyFraction = 0) const;
};


/**
 @class SamplerEnergyFraction
 @brief Throw away a fraction of the particles of a given type.
 The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 This is essentially the implementation of the thinning of the EM interactions. 
 In this case, the thinning is defined as: thinning = 1 - sampling
 */
class SamplerEnergyFraction: public Sampler {
	private:
		int particleId;
		double sampling;

	public:
		SamplerEnergyFraction(int particleId, double sampling);
		void setSampling(double sampling);
		void setParticleId(int particleId);
		double getSampling() const;
		int getParticleId() const;
		double computeWeight(int id, double energy = 0, double energyFraction = 0) const;
};

/**
 @class SamplerNull
 @brief Useful dummy class to be used as initializer.
 */
class SamplerNull : public Sampler {
	public:
		SamplerNull();
		double computeWeight(int id, double energy = 0, double energyFraction = 0) const;
};


/**
 @class SamplerList
 @brief List of objects of type "Sampler".
 */
class SamplerList : public Sampler {
	protected:
		std::vector< ref_ptr<Sampler> > samplers;
	public:
		SamplerList();
		SamplerList(std::vector<ref_ptr<Sampler>> samplers);
		void add(Sampler *sampler);
		inline void add(ref_ptr<Sampler> sampler) {
			add(sampler.get());
		}
		double computeWeight(int id, double energy = 0, double energyFraction = 0) const;
};

} // namespace crpropa

#endif