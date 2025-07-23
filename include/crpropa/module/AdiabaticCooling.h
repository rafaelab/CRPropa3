#ifndef CRPROPA_ADIABATICCOOLING_H
#define CRPROPA_ADIABATICCOOLING_H


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/advectionField/AdvectionField.h"
#include "kiss/logger.h"


namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
@class AdiabaticCooling
@brief Implements adiabatic cooling/heating due to advection.
*/

class AdiabaticCooling: public Module {
	protected:
		ref_ptr<AdvectionField> advectionField;
		double limit;

	public:
		/** Default constructor.
		 *  @param advectionField 	The advection field used for the adiabatic energy change
		 */
		AdiabaticCooling(ref_ptr<AdvectionField> advectionField);
		
		/** Constructor
		 *  @param advectionField 	the advection field used for the adiabatic energy change
		 *  @param limit 			maximum relative energy change allowed
		 */
		AdiabaticCooling(ref_ptr<AdvectionField> advectionField, double limit);
		
		void setLimit(double l);
		double getLimit() const;
		void process(Candidate *c) const;	
};
/** @}*/



}; // end namespace crpropa

#endif // CRPROPA_ADIABATICCOOLING_H
