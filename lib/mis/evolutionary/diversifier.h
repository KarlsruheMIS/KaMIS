/**
 * diversifier.h
 * Purpose: Refresh the RNG seed.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
