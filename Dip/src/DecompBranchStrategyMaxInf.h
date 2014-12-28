#ifndef DecompBranchStrategyMaxInf_h_
#define DecompBranchStrategyMaxInf_h_

#include "DecompBranchObject.h"
#include "DecompBranchStrategy.h"
#include "DecompAlgo.h"



/** This class implements maximum infeasibility branching. */
class DecompBranchStrategyMaxInf : public DecompBranchStrategy {

private:

   /** Illegal Assignment operator.*/
   DecompBranchStrategyMaxInf& operator=(const DecompBranchStrategyMaxInf& rhs);

public:

   /** MaxInf Constructor. */
   DecompBranchStrategyMaxInf() : DecompBranchStrategy() {
      type_ = static_cast<int>(DecompBranchingStrategyMaxInfeasibility);
   }

   /** MaxInf Constructor. */
   DecompBranchStrategyMaxInf(DecompAlgo* model) : DecompBranchStrategy(model) {
      type_ = static_cast<int>(DecompBranchingStrategyMaxInfeasibility);
   }

   /** Destructor. */
   virtual ~DecompBranchStrategyMaxInf() {}

   /** Copy constructor. */
   DecompBranchStrategyMaxInf(const DecompBranchStrategyMaxInf&);

   DecompBranchStrategyMaxInf(const DecompAlgo&);

   /** Clone a brancing strategy. */
   virtual DecompBranchStrategy* clone() const {
      return new DecompBranchStrategyMaxInf(*this);
   }

   /** Create a set of candidate branching objects. */
   virtual int createCandBranchObjects(int numPassesLeft, double ub);

   /** Compare branching object thisOne to bestSoFar. If thisOne is better
   	than bestObject, return branching direction(1 or -1), otherwise
   	return 0.
   	If bestSorFar is NULL, then always return branching direction(1 or -1).
   */
   virtual int betterBranchObject(DecompBranchObject* thisOne,
                                  DecompBranchObject* bestSoFar);
};

#endif
