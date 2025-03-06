/**
 * mis_permutation.h
 * Purpose: Data structure used for the local search algorithm.
 *          The structure itself is a permutation of all nodes divided into three blocks.
 *          The first block contains the solution nodes.
 *          The second block contains free nodes.
 *          The third block contains all remaining (non-free, non-solution) nodes.
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

#ifndef _MIS_PERMUTATION_H_
#define _MIS_PERMUTATION_H_

#include <vector>

#include "definitions.h"
#include "graph_access.h"
#include "data_structure/candidate_list.h"

class mis_permutation {
    friend class ils;
    public:
        /**
         * Default Constructor.
         */
        mis_permutation();

        /**
         * Default Destructor.
         */
        virtual ~mis_permutation();

        /**
         * Constructs and initializes the permutation for the given graph.
         * This is done by calculating the tightness of each node.
         * 
         * @param G Graph representation.
         */
        void construct(graph_access & G);

        /**
         * Returns the tightness of a node.
         *
         * @param node Node for which the tightness should be returned.
         * @return Node tightness.
         */
        int get_tightness(NodeID node) { return tightness[node]; };

        /**
         * Returns the position of a node within the permutation.
         *
         * @param node Node for which the position should be returned.
         * @return Node position in the permutation.
         */
        unsigned int get_position(NodeID node) { return position[node]; };

        /**
         * Get the solution node at a certain position.
         * Therefore, the position should be smaller than the solution size.
         *
         * @param i The node position.
         * @return Solution node at the position.
         */
        NodeID get_solution_node(unsigned int i);

        /**
         * Get the free node at a certain position.
         * Therefore, the position should be smaller than the free-blocks size.
         *
         * @param i The node position.
         * @return Free node at the position.
         */
        NodeID get_free_node(unsigned int i);

        /**
         * Get the non free, non solution node at a certain position.
         * 
         * @param i The node position.
         * @return Non free node at the position.
         */
        NodeID get_non_free_node(unsigned int i);

        /**
         * Get the non solution node at a certain position.
         *
         * @param i The node position.
         * @return Non solution node at the position.
         */
        NodeID get_non_solution_node(unsigned int i);

        /**
         * Get the current size of the solution.
         *
         * @return Solution size.
         */
        unsigned int get_solution_size();

		/**
		* Get the current weight of the solution.
		*
		* @return Solution weight.
		*/
		NodeWeight get_solution_weight();

        /**
         * Get the current amount of free nodes.
         *
         * @return Amount of free nodes.
         */
        unsigned int get_free_size();

        /**
         * Check whether or not the solution is maximal.
         * Maximal means there are no free nodes left.
         *
         * @return True if the solution is maximal.
         */
        bool is_maximal();

        /**
         * Add a node to the current solution.
         * Any neighbors that are invalid after the insertion will be removed
         * from the solution.
         *
         * @param node Node to be inserted.
         * @param G Graph representation.
         */
        void add_to_solution(NodeID node, graph_access & G);

        /**
         * Remove a node from the current solution.
         * Any neighbors that are changed by the removal will be moved 
         * to the according blocks.
         *
         * @param node Node to be removed.
         * @param G Graph representation.
         */
        void remove_from_solution(NodeID node, graph_access & G);

        /**
         * Check if the given node is in the solution.
         *
         * @return True if the node is in the solution.
         */
        bool is_solution_node(NodeID node);

        /**
         * Check if the given node is free.
         *
         * @return True if the node is free.
         */
        bool is_free_node(NodeID node);
        
        /**
         * Check if the given node is non-free and not in the solution.
         *
         * @return True if the node is non free and not in the solution.
         */
        bool is_non_solution_node(NodeID node);

        /**
         * Print infomation about the current state of the permutation.
         * The detailed version prints the nodes within each block.
         *
         * @param details True if the output should be more detailed.
         */
        void print(bool details);
        
        /**
         * Print the positions of all nodes.
         */
        void print_position();

        /**
         * Print the tightness of all nodes.
         */
        void print_tightness();

        /**
         * Check if the permutation is valid.
         * Validity means that all nodes are in the correct block
         * according to their tightness.
         *
         * @return True if the permutation is valid.
         */
        bool check_permutation();

        /**
         * Check if there are any inconsistencies.
         *
         * @return True if no inconsistencies were found.
         */
        bool check_consistency(graph_access & G);


	    size_t added_vertices = 0;

    private:
		// Pointer to graph_access in order to access ndoe weights
		graph_access* G_ptr;
        // The actual permutation.
        std::vector<NodeID> nodes;
        // Array containing the position of each node.
        std::vector<unsigned int> position;
        // Array containing the tightness of each node.
        std::vector<int> tightness;
        // Size of the solution.
        unsigned int solution_size;
        // Amount of free nodes.
        unsigned int free_size;
        // Amount of remaining nodes.
        unsigned int remaining_size;
        // Total size of the permutation.
        unsigned int total_size;
        // Number of inconsistencies.
        unsigned int inconsistencies;
        // List of onetight nodes.
        candidate_list onetight_all;

        /**
         * Calculate the tightness for a node.
         *
         * @param node Node for which the tightness should be calculated.
         * @param G Graph representation.
         * @return Tightness of the node.
         */
        int calculate_tightness(NodeID node, graph_access & G);

        /**
         * Get the size of the current permutation.
         *
         * @return Permutation size.
         */
        int get_size() { return solution_size + free_size + remaining_size; }

        /**
         * Move a node to the solution block.
         *
         * @param Node Node to be moved.
         * @param G Graph representation.
         */
        void move_to_solution(NodeID node, graph_access & G);

        /**
         * Move a node to the free block.
         *
         * @param Node Node to be moved.
         * @param G Graph representation.
         */
        void move_to_free(NodeID node, graph_access & G);

        /**
         * Move a node to the non-free, non-solution block.
         *
         * @param Node Node to be moved.
         * @param G Graph representation.
         */
        void move_to_non_free(NodeID node, graph_access & G);

        /**
         * Swap the position of the two nodes in the permutation.
         *
         * @param first_pos The position of the first node.
         * @param second_pos The position of the second node.
         */
        void swap_nodes(unsigned int first_pos, unsigned int second_pos);
};

#endif
