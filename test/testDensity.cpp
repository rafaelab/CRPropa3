#include "crpropa/Grid.h"
#include "crpropa/Units.h"
#include "crpropa/massDistribution/TargetMedium.h"
#include "crpropa/massDistribution/ConstantDensity.h"
#include "crpropa/massDistribution/Cordes.h"
#include "crpropa/massDistribution/Ferriere.h"
#include "crpropa/massDistribution/Nakanishi.h"

#include "gtest/gtest.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace crpropa {

TEST(testConstantDensity, SimpleTest) {
	TargetMedium targetMedium("HI", "H", 1.0);
	ConstantDensity* density = new ConstantDensity(1.234, targetMedium);
	Vector3d position(1, 2, 3);
	EXPECT_DOUBLE_EQ(density->getDensity(position), 1.234);
}

TEST(testCordes, checkValueAtCertainPoints) {
	CordesHII n;
	Vector3d position(3.1 * kpc, 2.9 * kpc, -30 * pc);	//position for testing density
	EXPECT_NEAR(n.getDensity(position), 184500., 1);
}

TEST(testNakanishi, checkValueAtCertainPoints) {
	NakanishiHI nHI;
	NakanishiH2 nH2;

	Vector3d position(4 * kpc, -2.5 * kpc, -0.85 * kpc);
	EXPECT_NEAR(nHI.getDensity(position), 914., 1.);
	EXPECT_NEAR(nH2.getDensity(position), 0., 1.);

	position = Vector3d(50 * pc, 100 * pc, 10 * pc);
	EXPECT_NEAR(nHI.getDensity(position), 540867., 1.);
	EXPECT_NEAR(nH2.getDensity(position), 10335137., 1.);
}

TEST(testFerriere, checkValueAtCertainPoints) {
	FerriereHI nHI;
	FerriereHII nHII;
	FerriereH2 nH2;

	//testing density in inner Ring (R <= 3*kpc)
	Vector3d position1(60 * pc, -60 * pc, -20 * pc);	//testing position in region of CMZ
	Vector3d transf = FerriereBase::transformCoordinatesCMZ(position1);
	EXPECT_NEAR(transf.x, 5.9767 * pc, 1e-4 * pc);
	EXPECT_NEAR(transf.y, 12.8171 * pc, 1e-4 * pc);
	EXPECT_DOUBLE_EQ(transf.z, position1.z);	//no transformation in z component
	transf = FerriereBase::transformCoordinatesDisk(position1);
	EXPECT_NEAR(transf.x,  11.0660 * pc, 1e-4 * pc);
	EXPECT_NEAR(transf.y,  82.5860 * pc, 1e-4 * pc);
	EXPECT_NEAR(transf.z, -25.6338 * pc, 1e-4 * pc);

	// testing density
	EXPECT_NEAR(nHI.getDensity(position1), 6237723, 1.); 
	EXPECT_NEAR(nHII.getDensity(position1), 6243793, 1.);
	EXPECT_NEAR(nH2.getDensity(position1), 35484825, 1.);

	//testing position in region of the disc
	Vector3d position2(-500 * pc, -900 * pc, 35 * pc);
	EXPECT_NEAR(nHI.getDensity(position2), 5., 1.);
	EXPECT_NEAR(nHII.getDensity(position2), 48190, 1.);
	EXPECT_NEAR(nH2.getDensity(position2), 0., 1.);

	//testing the outer region R>3kpc and  3 kpc < R < R_sun
	Vector3d position3(5 * kpc, 4 * kpc, -29 * pc);
	EXPECT_NEAR(nHI.getDensity(position3), 540607, 1);
	EXPECT_NEAR(nHII.getDensity(position3), 66495, 1);
	EXPECT_NEAR(nH2.getDensity(position3), 2492685, 1);

	//testing the outer region R>3kpc and  R > R_sun
	Vector3d position4(10 * kpc, 2 * kpc, 50 * pc);
	EXPECT_NEAR(nHI.getDensity(position4), 431294, 1.);
	EXPECT_NEAR(nHII.getDensity(position4), 22109, 1.);
	EXPECT_NEAR(nH2.getDensity(position4), 54099, 1.);
}


TEST(testGridDensity, SimpleTest) {
	size_t Nx = 5;
	size_t Ny = 8;
	size_t Nz = 10;
	double spacing = 2.0;
	Vector3d origin(1., 2., 3.);

	ref_ptr<Grid1f> grid = new Grid1f(origin, Nx, Ny, Nz, spacing);

	// set some values for the grid
	grid->get(3, 2, 4) = 5;
	grid->get(3, 2, 5) = 12;
	grid->get(2, 3, 4) = 6;

	DensityGrid dens = DensityGrid(grid);

	// a point in the region where values are defined for the grid.
	Vector3d position = origin + Vector3d(2.2, 2.8, 4.1) * spacing; 
	double valueFromGrid =  grid->interpolate(position);
	double nTotal = dens.getDensity(position);

	// Check for values
	EXPECT_DOUBLE_EQ(valueFromGrid, nTotal);
	EXPECT_DOUBLE_EQ(nTotal, valueFromGrid);
}


} //namespace crpropa
