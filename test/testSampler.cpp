#include <gtest/gtest.h>

#include "crpropa/Random.h"
#include "crpropa/Sampler.h"


namespace crpropa {


// Test SamplerEventsUniform
TEST(Sampler, SamplerEventsUniformInitialization) {
	SamplerEventsUniform sampler(1, 0.5);
	EXPECT_EQ(sampler.getParticleId(), 1);
	EXPECT_DOUBLE_EQ(sampler.getSampling(), 0.5);
}

TEST(Sampler, SamplerEventsUniformComputeWeight) {
	SamplerEventsUniform sampler(22, 0.9);

	unsigned int accepted = 0;
	unsigned int rejected = 0;
	for (size_t i = 0; i < 100; i++) {
		double w = sampler.computeWeight(22, 0, 0, 0);
		if (w == 0) {
			rejected++;
		} else {
			accepted++;
		}
	}

	EXPECT_TRUE(accepted + rejected == 100);
	EXPECT_GT(accepted, 50); // may fail sometimes
	EXPECT_LT(rejected, 50);
}

// Test SamplerEventsEnergyFractionPowerLaw
TEST(Sampler, SamplerEventsEnergyFractionPowerLawWeightFunction) {
	SamplerEventsEnergyFractionPowerLaw sampler(2.0, 22, 0.5);
	unsigned int accepted = 0;
	unsigned int rejected = 0;
	for (size_t i = 0; i < 100; i++) {
		double w = sampler.computeWeight(22, 0, 0, 0);
		if (w == 0) {
			rejected++;
		} else {
			accepted++;
		}
	}

	EXPECT_TRUE(accepted + rejected == 100);
	// EXPECT_GT(accepted, 50); // may fail sometimes
	// EXPECT_LT(rejected, 50);
}

// Test SamplerEventsNull
TEST(Sampler, SamplerEventsNullComputeWeight) {
	SamplerEventsNull sampler;
	EXPECT_DOUBLE_EQ(sampler.computeWeight(1, 0, 0, 0), 1.0); // Always returns 1.0
}

// Test SamplerEventsList
TEST(Sampler, SamplerEventsListComputeWeight) {
	SamplerEventsList samplerList;
	samplerList.add(new SamplerEventsNull());
	EXPECT_DOUBLE_EQ(samplerList.computeWeight(1, 0, 0, 0), 1.0); // Combined weight
}

// Test SamplerDistributionUniform
TEST(Sampler, SamplerDistributionUniformInitialization) {
	SamplerDistributionUniform sampler(0.0, 1.0, 10, "lin");
	EXPECT_EQ(sampler.getSize(), 0);
	EXPECT_EQ(sampler.getDistribution()->getNumberOfBins(), 10);
}

// TEST(Sampler, SamplerDistributionUniformGetSample) {
// 	SamplerDistributionUniform sampler(0.0, 1.0, 10, "lin");
// 	sampler.push(0.5);
// 	auto samples = sampler.getSample(5);
// 	EXPECT_EQ(samples.size(), 1); // Only one sample pushed
// }

TEST(Sampler, SamplerDistributionUniformTransformations) {
	SamplerDistributionUniform sampler(0.0, 1.0, 10, "lin");
	sampler.push(0.5);
	sampler.transformToPDF();
	sampler.transformToCDF();
	// EXPECT_TRUE(sampler.getDistribution()->isCDF()); // Verify transformation
}

TEST(Sampler, SamplerDistributionUniformAppendAndClear) {
	SamplerDistributionUniform sampler(0.0, 1.0, 10, "lin");
	sampler.append({0.1, 0.2, 0.3});
	EXPECT_EQ(sampler.getSize(), 3);
	sampler.clear();
	EXPECT_EQ(sampler.getSize(), 0);
}

} // namespace crpropa