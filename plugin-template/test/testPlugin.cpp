#include "CRPropa.h"
#include "myPlugin.h"

#include "gtest/gtest.h"

// with this expression we do not need to write crpropa:: 
// for everything out of the crpropa namespace
using namespace crpropa;

// create a test over TEST(ModuleToTest, WhatToTest){}
TEST(MyModule, SimpleTest){

	ModuleList SIM;

	SIM.add(new SimplePropagation(1*pc, 1*pc));
	SIM.add(new MaximumTrajectoryLength(1000*pc));
	SIM.add(new myPlugin::MyModule());

	// define source to test this module and get candidate with correct
	// flag out of there
	Source source;
	// AddMyProperty is added in plugin.h but under the crpropa namespace:
	source.add(new AddMyProperty());
	ref_ptr<Candidate> Cand = source.getCandidate();
	EXPECT_EQ(Cand->getProperty("counter"), uint32_t(0));

	SIM.run(Cand);

	// in a test you can test your values like this:
	EXPECT_EQ(Cand->getProperty("counter"), uint32_t(1000));
	// but there are many more such EXPECT expressions
}

int main(int argc, char **argv) {
	// this runs all above defined tests
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}