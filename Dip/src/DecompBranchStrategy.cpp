#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>

#include "DecompBranchStrategy.h"
#include "AlpsDecompModel.h"

//#############################################################################


/* Compare N branch objects and identify bestObject_. Return index
   of best and sets way of branch bestObject_. */
DecompBranchObject*
DecompBranchStrategy::bestBranchObject()
{
   int i, betterDir;
   int bestDir = 0;
   int bestBrObjIndex = -1;

   if (numBranchObjects_ > 1) {
      //--------------------------------------------------
      // Clear best members.
      //--------------------------------------------------
      clearBest(algo_);
      //--------------------------------------------------
      // Select the best branching object.
      //--------------------------------------------------

      for (i = 0; i < numBranchObjects_; ++i) {
         betterDir = betterBranchObject(branchObjects_[i],
                                        bestBranchObject_);

         if (betterDir) {
            bestBrObjIndex = i;
            bestBranchObject_ = branchObjects_[i];
            bestDir = betterDir;
         }
      }

      if (bestBrObjIndex >= 0) {
         // Set branching direction.
         bestBranchObject_->setDirection(bestDir);
         // Need move this b obj to node later. Rest will be deleted.
         branchObjects_[bestBrObjIndex] = NULL;
      } else {
         bestBranchObject_ = NULL;
      }

      //--------------------------------------------------
      // Delete rest candidates.
      //--------------------------------------------------

      for (i = 0; i < numBranchObjects_; ++i) {
         if (branchObjects_[i]) {
            delete branchObjects_[i];
            branchObjects_[i] = NULL;
         }
      }
   } else {
      bestBranchObject_ = branchObjects_[0];
   }

   delete [] branchObjects_;
   branchObjects_ = NULL;
   numBranchObjects_ = 0;
   return bestBranchObject_;
}

//#############################################################################
