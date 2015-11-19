/**
 * operation_log.h
 * Purpose: Singleton deque to store the operations performed during ILS.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
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

#ifndef _OPERATION_LOG_H_
#define _OPERATION_LOG_H_

#include <deque>

#include "definitions.h"

class operation_log {
    public:
        /**
         * Default Constructor.
         */
        operation_log();

        /**
         * Default Destructor.
         */
        virtual ~operation_log();

        /**
         * Get the singleton instance.
         * 
         * @return Instance of the logger.
         */
        static operation_log *instance() {
            static operation_log inst;
            return &inst;
        };
        
        /**
         * Initialize the deque with the given size
         *
         * @param size The size of the deque.
         */
        void init(unsigned int size);
        
        /**
         * Log a insertion of the given node.
         *
         * @param x Node that was inserted.
         */
        void report_insert(NodeID x);

        /**
         * Log a removal of the given node.
         *
         * @param x Node that was removed.
         */
        void report_remove(NodeID x);

        /**
         * Look at a certain position of the deque.
         *
         * @param pos Position to look at.
         * @return Node at the position.
         */
        int peek(unsigned int pos);

        /**
         * Return the latest operation performed.
         *
         * @return Performed operation.
         */
        int unwind();

        /** 
         * Check if the deque is empty.
         *
         * @return True if empty.
         */
        bool is_empty();

        /**
         * Check if the deque is full.
         * 
         * @return True if full.
         */
        bool is_full();

        /**
         * Get the current size of the deque.
         *
         * @return Size of the deque.
         */
        unsigned int get_size();

        /**
         * Log is now active.
         */
        void activate();

        /**
         * Log is now inactive.
         */
        void deactivate();

        /**
         * Reset the log.
         */
        void reset();

        /**
         * Print the current status of the log.
         */
        void print();

    private:
        // Whether or not logging is active.
        bool active;
        // Underlying deque.
        std::deque<int> op_log;
};

#endif

