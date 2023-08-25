#ifndef CRPROPA_PARTICLE_ID_H
#define CRPROPA_PARTICLE_ID_H

#include <cstddef>
#include <vector>
#include <string>

#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"
#include "kiss/convert.h"



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
 * Ion numbers are +/- 10LZZZAAAI.
 * AAA is A - total baryon number
 * ZZZ is Z - total charge
 * L is the total number of strange quarks.
 * I is the isomer number, with I=0 corresponding to the ground state.
 */
int nucleusId(int a, int z);
int chargeNumber(int id);
int massNumber(int id);

bool isNucleus(int id);


/** Additional functionalities.
 *  This redundant implementation (as opposed to exposing the functions from the HepPID namespace) prevents compiler issues.
*/
bool isMeson(int id);
bool isBaryon(int id);
bool isHadron(int id);
bool isLepton(int id);
bool hasUp(int id);
bool hasDown(int id);
bool hasCharm(int id);
bool hasStrange(int id);
bool hasTop(int id);
bool hasBottom(int id);

/* Additional modules */
std::string convertIdToName(int id); 

/** @}*/
} // namespace crpropa

#endif // CRPROPA_PARTICLE_ID_H
