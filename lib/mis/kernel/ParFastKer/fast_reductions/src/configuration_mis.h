/**
 * configuration.h
 * Purpose: Contains preset configurations for the evolutionary algorithms.
 *
 *****************************************************************************/

#ifndef _CONFIGURATION_MIS_H_
#define _CONFIGURATION_MIS_H_

#include "definitions.h"
#include "mis_config.h"
#include "interface/kaHIP_interface.h"

class configuration_mis {
    public:
        /**
         * Default Constructor.
         */
        configuration_mis() {} ;

        /**
         * Default Destructor.
         */
        virtual ~configuration_mis() {};

        /**
         * Set the standard configuration.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void standard( MISConfig & config );

        /**
         * Set the configuration for social network graphs.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void social( MISConfig & config );

        /**
         * Set the configuration for the experimental evaluation.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void full_standard( MISConfig & config );

        /**
         * Set the configuration for the experimental evaluation 
         * for social network graphs.
         * Use local search for combine operations.
         * Use ILS for improving offsprings.
         *
         * @param config Config to be initialized.
         */
        void full_social( MISConfig & config );
};

inline void configuration_mis::standard( MISConfig & mis_config ) {
    // Reductions
    mis_config.all_reductions                         = true;
    mis_config.num_reps                               = 5;
}

inline void configuration_mis::social( MISConfig & mis_config ) {
    standard(mis_config);
}

inline void configuration_mis::full_standard( MISConfig & mis_config ) {
    standard(mis_config);
}

inline void configuration_mis::full_social( MISConfig & mis_config ) {
    full_standard(mis_config);
}

#endif 
