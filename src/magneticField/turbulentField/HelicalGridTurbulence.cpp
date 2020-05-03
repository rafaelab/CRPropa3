#include "crpropa/magneticField/turbulentField/HelicalGridTurbulence.h"
#include "crpropa/magneticField/turbulentField/GridTurbulence.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#ifdef CRPROPA_HAVE_FFTW3F
#include "fftw3.h"

namespace crpropa {

HelicalGridTurbulence::HelicalGridTurbulence(double Brms, double lMin, double lMax, double sindex, unsigned int seed, double H)
	: TurbulentField(Brms, sindex), lMin(lMin), lMax(lMax), seed(seed) {
	initGrid();
	initTurbulence(gridPtr, Brms, lMin, lMax, - sindex - 2, seed, H);
}

HelicalGridTurbulence::HelicalGridTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax, double alpha, unsigned int seed, double H) 
	: TurbulentField(Brms, sindex), gridPtr(grid), lMin(lMin), lMax(lMax), seed(seed) {
	initTurbulence(gridPtr, Brms, lMin, lMax, - sindex - 2, seed, H);
}

Vector3d HelicalGridTurbulence::getField(const Vector3d& pos) const {
	return gridPtr->interpolate(pos);
}

double HelicalGridTurbulence::getCorrelationLength() const {
	return turbulentCorrelationLength(lMin, lMax, -sindex - 2);
}

double HelicalGridTurbulence::turbulentCorrelationLength(double lMin, double lMax, double s) {
	double r = lMin / lMax;
	return lMax / 2 * (s - 1) / s * (1 - pow(r, s)) / (1 - pow(r, s - 1));
}

void HelicalGridTurbulence::initGrid() {
	gridPtr = new Grid3f(Vector3d(-boxSize/2), gridSize, spacing);
}

void HelicalGridTurbulence::initTurbulence(ref_ptr<Grid3f> grid, double Brms, double lMin, double lMax, double alpha, int seed, double H) {

	checkGridRequirementsTEMP(grid, lMin, lMax);

	Vector3d spacing = grid->getSpacing();
	size_t n = grid->getNx(); // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k;

	// only used if there is a helicity
	double Bktot, Bkplus, Bkminus, thetaplus, thetaminus;

	double kMin = spacing.x / lMax;
	double kMax = spacing.x / lMin;
	Vector3f b; // real b-field vector
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getR();

				// wave outside of turbulent range -> B(k) = 0
				if ((k < kMin) || (k > kMax)) {
					Bkx[i][0] = 0;
					Bkx[i][1] = 0;
					Bky[i][0] = 0;
					Bky[i][1] = 0;
					Bkz[i][0] = 0;
					Bkz[i][1] = 0;
					continue;
				}

				// construct an orthogonal base ek, e1, e2
				// (for helical fields together with the real transform the following convention
				// must be used: e1(-k) = e1(k), e2(-k) = - e2(k)
				if (ek.getAngleTo(n0) < 1e-3) { // ek parallel to (1,1,1)
					e1.setXYZ(-1, 1, 0);
					e2.setXYZ(1, 1, -2);
				} else { // ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.getR();
				e2 /= e2.getR();


				double Bkprefactor = mu0 / (4 * M_PI * pow(k, 3));
				Bktot = fabs(random.randNorm() * pow(k, alpha / 2));
				Bkplus  = Bkprefactor * sqrt((1 + H) / 2) * Bktot;
				Bkminus = Bkprefactor * sqrt((1 - H) / 2) * Bktot;
				thetaplus = 2 * M_PI * random.rand();
				thetaminus = 2 * M_PI * random.rand();
				double ctp = cos(thetaplus);
				double stp = sin(thetaplus);
				double ctm = cos(thetaminus);
				double stm = sin(thetaminus);

				Bkx[i][0] = ((Bkplus * ctp + Bkminus * ctm) * e1.x + (-Bkplus * stp + Bkminus * stm) * e2.x) / sqrt(2);
				Bkx[i][1] = ((Bkplus * stp + Bkminus * stm) * e1.x + ( Bkplus * ctp - Bkminus * ctm) * e2.x) / sqrt(2);
				Bky[i][0] = ((Bkplus * ctp + Bkminus * ctm) * e1.y + (-Bkplus * stp + Bkminus * stm) * e2.y) / sqrt(2);
				Bky[i][1] = ((Bkplus * stp + Bkminus * stm) * e1.y + ( Bkplus * ctp - Bkminus * ctm) * e2.y) / sqrt(2);
				Bkz[i][0] = ((Bkplus * ctp + Bkminus * ctm) * e1.z + (-Bkplus * stp + Bkminus * stm) * e2.z) / sqrt(2);
				Bkz[i][1] = ((Bkplus * stp + Bkminus * stm) * e1.z + ( Bkplus * ctp - Bkminus * ctm) * e2.z) / sqrt(2);

				Vector3f BkRe(Bkx[i][0], Bky[i][0], Bkz[i][0]);
				Vector3f BkIm(Bkx[i][1], Bky[i][1], Bkz[i][1]);
			} // for iz
		} // for iy
	} // for ix

	executeInverseFFTInplaceTEMP(grid, Bkx, Bky, Bkz);

	scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}

} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F
