import sys

try:
    import unittest
except:
    print("***********************************************************") 
    print("* WARNING!! Couldn't import python unittesting framework! *") 
    print("* No python tests have been executed                      *") 
    print("***********************************************************") 
    sys.exit(0)

try:
    import crpropa 
except Exception as e:
    print("*** CRPropa import failed")
    print(type(e), str(e))
    sys.exit(-1)

class TestMagneticLens(unittest.TestCase):
    def setUp(self):
        self.magneticLens = crpropa.MagneticLens(5)
        p = crpropa.Pixelization(5)
        m = crpropa.ModelMatrix(p.nPix(5), p.nPix(5), p.nPix(5))
        for i in range(p.nPix(5)):
            theta = 0.
            phi = 0.
            [phiout, thetaout] = p.pix2Direction(i, phi, theta)
            phiout += 0.4
            if phiout > 3.14:
                phiout -= 3.14
            j = p.direction2Pix(phiout, thetaout)
            m.setElement(i, j, 1)
        self.magneticLens.setLensPart(m, 18, 20)

    def test_transformCosmicRays(self):
        phi = 0
        theta = -0.01
        [r, phi, theta] = self.magneticLens.transformCosmicRay(19.0, phi,
                                                               theta)
        self.assertTrue(r)
        self.assertFalse(theta == -0.01)

    def test_transformModelVector_numpyArray(self):
        import numpy
        v = numpy.ones(self.magneticLens.getPixelization().nPix())
        self.magneticLens.transformModelVector(v, 2)

    def test_accessToPixelization(self):
        self.assertEquals(self.magneticLens.getPixelization().nPix(),
                          self.magneticLens.getPixelization().nPix(5))


class testPixelizationConsistency(unittest.TestCase):
    def testConsistencyWithHealpy(self):
        try:
            import healpy
        except ImportError:
            print "Consistency with healpy not tested as healpy is" \
                  "not available"
            return
        p = crpropa.Pixelization(4)
        from numpy import linspace
        from math import pi
        for theta in linspace(0, pi):
            for phi in linspace(-pi, pi):
                crpropaIdx = p.direction2Pix(phi, pi / 2 - theta)
                hpIdx = healpy.ang2pix(2 ** 4, theta, phi)
                self.assertEqual(crpropaIdx, hpIdx)


if __name__ == '__main__':
    unittest.main()
