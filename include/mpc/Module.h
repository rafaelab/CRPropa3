#ifndef MPC_MODULE_H_
#define MPC_MODULE_H_

#include "mpc/Candidate.h"
#include "mpc/Referenced.h"
#include "mpc/Common.h"

#include <string>

namespace mpc {

class Candidate;

/**
 @class Module
 @brief Abstract base class for modules
 */
class Module: public Referenced {
protected:
	std::string description;

public:
	Module();
	virtual ~Module() {
	}
	virtual std::string getDescription() const;
	virtual void setDescription(const std::string &description);
	virtual void process(Candidate *candidate) const = 0;
	inline void process(ref_ptr<Candidate> candidate) const {
		process(candidate.get());
	}
};

class SimulationVolumeDependentModule: public Module {
public:
	virtual void updateSimulationVolume(const Vector3 &origin, double size) = 0;
};

} // namespace mpc

#endif /* MPC_MODULE_H_ */
