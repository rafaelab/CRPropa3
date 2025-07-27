#pragma once

#include <algorithm>
#include <vector>

#include "crpropa/Common.h"
#include "crpropa/Histogram.h"
#include "crpropa/Random.h"
#include "crpropa/Referenced.h"



namespace crpropa {

/**
 * \addtogroup Sampler
 * @{
 */

/**
 * @class SamplerEvents
 * @brief Abstract base class to handle individual events.
 */
class SamplerEvents: public Referenced {
	public:
		[[nodiscard]] virtual double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const = 0;
};


/**
 * @class SamplerEventsUniform
 * @brief Throw away a fraction of the particles of a given type.
 * The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 * Note that the energy fraction is implemented through setters and getters, but it will constantly be modified.
 * This is because the structure of computeWeight should remain constant (or variadic arguments used).
 * This will be done in the future.
 */
class SamplerEventsUniform: public SamplerEvents {
	protected:
		int particleId;
		double sampling;

	public:
		SamplerEventsUniform(int particleId, double sampling);
		void setSampling(double sampling);
		void setParticleId(int particleId);
		[[nodiscard]] double getSampling() const;
		[[nodiscard]] int getParticleId() const;
		[[nodiscard]] double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const override;
};


/**
 * @class SamplerEventsEnergy
 * @brief Throw away a fraction of the particles of a given type according to a given distribution.
 * This distribution is defined from the derived classes.
 * The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 * This is essentially the implementation of the thinning of the EM interactions for index=1.
 * In this case, the thinning is defined as: thinning = 1 - sampling.
 */
class SamplerEventsEnergy : public SamplerEvents {
	protected:
		int particleId;
		double sampling;

	public:
		SamplerEventsEnergy(int particleId, double sampling);
		void setSampling(double sampling);
		void setParticleId(int particleId);
		[[nodiscard]] double getSampling() const;
		[[nodiscard]] int getParticleId() const;
		[[nodiscard]] double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const override;
		[[nodiscard]] virtual double weightFunction(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const = 0;
};


/**
 * @class SamplerEventsEnergyFractionPowerLaw
 * @brief Throw away a fraction of the particles of a given type according to a power law associated with the energy fraction it takes.
 * The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 * This is essentially the implementation of the thinning of the EM interactions for index=1.
 * In this case, the thinning is defined as: thinning = 1 - sampling.

 */
class SamplerEventsEnergyFractionPowerLaw: public SamplerEventsEnergy {
	protected:
		int particleId;
		double sampling;
		double index;

	public:
		SamplerEventsEnergyFractionPowerLaw(double index, int particleId, double sampling);
		void setIndex(double index);
		[[nodiscard]] double getIndex() const;
		[[nodiscard]] double weightFunction(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const override;
};


// /**
//  @class SamplerEventsEnergyNormal
//  @brief Throw away a fraction of the particles of a given type according to a normal distribution.
//  */
// class SamplerEventsEnergyNormal : public SamplerEventsEnergy {
// 	protected:
// 		double mean;
// 		double standardDeviation;

// 	public:
// 		SamplerEventsEnergyNormal(int particleId, double sampling, double mean, double standardDeviation);
// 		void setStandardDeviation(double sigma);
// 		void setMean(double mean);
// 		double getMean() const;
// 		double getStandardDeviation() const;
// 		double weightFunction(int id, double energy = 0, double energyFraction = 0, int counter = 0) const;
// };


// /**
//  @class SamplerEventsEnergyLogNormal
//  @brief Throw away a fraction of the particles of a given type according to a log-normal distribution.
//  */
// class SamplerEventsEnergyLogNormal : public SamplerEventsEnergy {
// 	protected:
// 		double mean;
// 		double standardDeviation;

// 	public:
// 		SamplerEventsEnergyLogNormal(int particleId, double sampling, double mean, double standardDeviation);
// 		void setStandardDeviation(double sigma);
// 		void setMean(double mean);
// 		double getMean() const;
// 		double getStandardDeviation() const;
// 		double weightFunction(int id, double energy = 0, double energyFraction = 0, int counter = 0) const;
// };


/**
 * @class SamplerEventsEnergyFraction
 * @brief Throw away a fraction of the particles of a given type.
 * The sampling parameter controls whether all (sampling=1) or no (sampling=0) particles are accepted.
 * This is essentially the implementation of the thinning of the EM interactions. 
 * In this case, the thinning is defined as: thinning = 1 - sampling
 */
class SamplerEventsEnergyFraction: public SamplerEventsEnergyFractionPowerLaw {
	public:
		SamplerEventsEnergyFraction(int particleId, double sampling);
};


/**
 * @class SamplerEventsNull
 * @brief Useful dummy class to be used as initializer.
 */
class SamplerEventsNull : public SamplerEvents {
	public:
		SamplerEventsNull();
		[[nodiscard]] double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const override;
};


/**
 * @class SamplerEventsList
 * @brief List of objects of type "SamplerEvents".
 */
class SamplerEventsList : public SamplerEvents {
	protected:
		std::vector<ref_ptr<SamplerEvents>> samplers;
	public:
		SamplerEventsList();
		SamplerEventsList(std::vector<ref_ptr<SamplerEvents>> samplers);
		void add(SamplerEvents* samplers);
		[[nodiscard]] double computeWeight(const int& id, const double& energy = 0, const double& energyFraction = 0, const int& counter = 0) const override;
};


/**
 * @class SamplerEventsBinnedDistribution
 * @brief Builds a histogram and only after it is built samples from it.
 */
class SamplerDistribution : public Referenced {
	public:
		virtual std::vector<double> getSample(int nSamples) const = 0;
		virtual int getSize() const = 0;
		virtual ref_ptr<Histogram1D> getDistribution() const = 0;
		virtual void transformToPDF() = 0;
		virtual void transformToCDF() = 0;
		virtual void append(const std::vector<double>& v) = 0;
		virtual void push(const double& v) = 0;
		virtual void clear() = 0;
};

class SamplerDistributionUniform : public SamplerDistribution {
	protected:
		ref_ptr<Histogram1D> distribution;
		int datasetSize;

	public:
		SamplerDistributionUniform(double vmin, double vmax, int nBins, std::string scale = "lin");
		void setSize(int size);
		int getSize() const override;
		void setDistribution(ref_ptr<Histogram1D> dist);
		[[nodiscard]] ref_ptr<Histogram1D> getDistribution() const override;
		[[nodiscard]] std::vector<double> getSample(int nSamples) const override;
		void transformToPDF() override;
		void transformToCDF() override;
		void append(const std::vector<double>& v) override;
		void push(const double& v) override;
		void clear() override;
};

/** @} 
 */

} // namespace crpropa
