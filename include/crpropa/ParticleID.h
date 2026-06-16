#ifndef CRPROPA_PARTICLE_ID_H
#define CRPROPA_PARTICLE_ID_H

#include <cstddef>
#include <vector>
#include <string>

namespace crpropa {
/**
 * \addtogroup PhysicsDefinitions
 * @{
 */

/**
 @file
 @brief 2012 Monte Carlo nuclear code scheme
 */

/** This implements the 2012 Monte Carlo nuclear code scheme.
 * Link: https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
 * Ion numbers are +/- 10LZZZAAAI.
 * AAA is A - total baryon number
 * ZZZ is Z - total charge
 * L is the total number of strange quarks.
 * I is the isomer number, with I=0 corresponding to the ground state.
 * 
 * Tip: To use electrons use 11 as ID directly 
 * @param A  Mass Number
 * @param Z  Charge Number
 */
int nucleusId(int A, int Z);
/** Returns the charge number of given ID
 * @param id  Particle ID after the 2012 Monte Carlo nuclear code scheme.
 */
int chargeNumber(int id);
/** Returns the mass number of given ID
 * @param id  Particle ID after the 2012 Monte Carlo nuclear code scheme.
 */
int massNumber(int id);
/** Checks whether given ID is nucleus or not
 * @param id  Particle ID after the 2012 Monte Carlo nuclear code scheme.
 */
bool isNucleus(int id);

/** Returns corresponding name of the particle with given ID
 * @param id  Particle ID after the 2012 Monte Carlo nuclear code scheme.
 */
std::string convertIdToName(int id); 

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PARTICLE_ID_H
