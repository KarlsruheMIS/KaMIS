/**
 * candidate_list.h
 * Purpose: Represents the list used to store the candidates for the local search.
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

#ifndef _CANDIDATE_LIST_H_
#define _CANDIDATE_LIST_H_

#include <vector>

#include "definitions.h"

class candidate_list {
    public:
        /**
         * Default Constructor.
         */
        candidate_list();

        /**
         * Default Destructor.
         */
        virtual ~candidate_list();

        /**
         * Initialize the candidate list.
         *
         * @param size Size of the list.
         */
        void init(unsigned int size);

        /**
         * Randomly permutate the current list.
         */
        void random_permute();

        /**
         * Swap the nodes at the given positions.
         *
         * @param i Position of the first node.
         * @param j Position of the second node.
         */
        void swap(unsigned int i, unsigned int j);

        /**
         * Remove a random element from the list.
         *
         * @return The removed candidate.
         */
        NodeID remove_random();
        
        /**
         * Get the node at the given position.
         *
         * @param pos Position of the node.
         * @return Node at the position.
         */
        NodeID pick(unsigned int pos);
        
        /** Pick and return a random node.
         *
         * @return Random node.
         */
        NodeID pick_random();

        /**
         * Insert the given node into the list.
         *
         * @param node Node to be inserted.
         */
        void insert(NodeID node);

        /** 
         * Remove the given node from the list.
         *
         * @param node Node to remove.
         */
        void remove(NodeID node);

        /**
         * Check if the list contains the given node.
         *
         * @param node Node that should be checked.
         * @return True if the node is in the list.
         */
        bool contains(NodeID node);

        /**
         * Check if the candidate list is empty.
         *
         * @return True if the list is empty.
         */
        bool is_empty();

        /**
         * Resets the list by removing all elements.
         */
        void reset();

        /**
         * Get the size of the candidate list.
         *
         * @return Number of candidates.
         */
        unsigned int get_size() { return count; };

        /**
         * Prints the candidate list.
         */
        void print();
        
    private:
        // The actual candidates.
        std::vector<NodeID> nodes;
        // Array containing the position of each node.
        std::vector<int> position;
        // Number of candidates.
        unsigned int count;

        /**
         * Removes the element at the given position.
         *
         * @return The node at the position.
         */
        NodeID remove_position(unsigned int position);
};

#endif

