/**
 * diversifier.h
 * Purpose: Refresh the RNG seed.
 *
 *****************************************************************************/

#ifndef _DIVERSIFIER_H_
#define _DIVERSIFIER_H_

#include "mis_config.h"
#include "random_functions.h"

class diversifier {
    public:
        /**
         * Default Constructor.
         */
        diversifier() {};

        /** 
         * Default Destructor.
         */
        virtual ~diversifier() {};

        /**
         * Set a new seed for the given config
         *
         * @param config Config to update.
         */
        void diversify(MISConfig & config) {
            config.seed = random_functions::nextInt(0, 10000);
            random_functions::setSeed(config.seed);
        }
};

#endif
