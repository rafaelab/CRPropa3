#ifndef CRPROPA_MODULE_H
#define CRPROPA_MODULE_H


#include <string>
#include <typeinfo>

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"
#include "crpropa/Common.h"


namespace crpropa {

class Candidate;

/**
 * @class Module
 * @brief Abstract base class for modules
 */
class Module: public Referenced {
	protected:
		std::string description;
		
	public:
		Module();
		virtual ~Module() = default;
		void setDescription(const std::string& description);
		virtual std::string getDescription() const;
		virtual void process(Candidate* candidate) const = 0;
		inline void process(ref_ptr<Candidate> candidate) const {
			process(candidate.get());
		}
};


/**
 * @class AbstractCondition
 * @brief Abstract Module providing common features for conditional modules.
 */
class AbstractCondition: public Module {
	protected:
		ref_ptr<Module> rejectAction; 
		ref_ptr<Module> acceptAction;
		bool makeRejectedInactive;
		bool makeAcceptedInactive;
		std::string rejectFlagKey;
		std::string rejectFlagValue;
		std::string acceptFlagKey;
		std::string acceptFlagValue;

		void reject(Candidate* candidate) const;
		inline void reject(ref_ptr<Candidate> candidate) const {
			reject(candidate.get());
		}

		void accept(Candidate* candidate) const;
		inline void accept(ref_ptr<Candidate> candidate) const {
			accept(candidate.get());
		}

	public:
		AbstractCondition();
		void onReject(Module* rejectAction);
		void onAccept(Module* acceptAction);
		void setMakeRejectedInactive(bool makeInactive);
		void setMakeAcceptedInactive(bool makeInactive);
		void setRejectFlag(std::string key, std::string value);
		void setAcceptFlag(std::string key, std::string value);

		// return the reject flag (key & value), delimiter is the "&".
		std::string getRejectFlag();

		// return the accept flag (key & value), delimiter is the "&"
		std::string getAcceptFlag();
};

/**
 * @class Deactivation
 * @brief Direct deactivation of the candidate. Can be used for debuging.
 */
class Deactivation: public AbstractCondition {
	public: 
		void process(Candidate* candidate) const override { 
			reject(candidate); 
		}
};


} // namespace crpropa

#endif /* CRPROPA_MODULE_H */
