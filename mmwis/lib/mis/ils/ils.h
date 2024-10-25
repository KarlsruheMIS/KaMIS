/**
 * ils.h
 * Purpose: Perform the iterated local search (ILS) as described by Andrade et al.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _ILS_H_
#define _ILS_H_

#include <vector>

#include "timer.h"
#include "mmwis_config.h"
#include "local_search.h"
#include "graph_access.h"

namespace mmwis {

class ils {
       public:
        /**
         * Default Constructor.
         */
        ils(const mmwis::MISConfig &config);

        /**
         * Default Destructor.
         */
        virtual ~ils();

        /**
         * Main algorithm of the ILS.
         * More detailed information can be found in the original paper
         * of Andrade et al.
         *
         * @param config Config used for the ILS.
         * @param G Graph representation.
         * @param iteration_limit Maximum number of iterations.
         */
        void perform_ils(graph_access &G, unsigned int iteration_limit = 0, int offset = 0);

        /** 
         * Reset the ILS.
         * Clears the force-list and best solution.
         */
        void reset();

       private:
        // Config for ils
        mmwis::MISConfig config;
        // Array for storing the last time a node was forced in the solution.
        std::vector<NodeID> last_forced;
        // Main solution data structure.
        mis_permutation *perm;
        // List of candidates for the local search.
        maxNodeHeap *cand;
        // List of nodes that were forced in the solution.
        candidate_list *force_list;
        // List of onetight nodes.
        candidate_list *one;
        // Best solution found so far
        NodeID *best_solution;
        NodeID best_solution_size = 0;
        NodeWeight best_solution_weight = 0;
        NodeID best_solution_size_test = 0;
        NodeWeight best_solution_weight_test = 0;
        // Local search algorithm
        local_search local;
        // Timer for measuring the time taken for the ILS.
        timer t;

        // ILS config
        unsigned int plateau_down;
        unsigned int plateau_up;
        unsigned int plateau_best;
        unsigned int plateau_best_again;
        unsigned int pden_floor;
        unsigned int delta_penalty;
        bool limit_plateau;
        bool swap_on_failure;
        unsigned int used_iterations = 0;
        /**
         * Force the given node in the solution.
         *
         * @param config Config used for the ILS.
         * @param G Graph representation.
         * @param v Node to force into the solution.
         * @param force_list List storing the forced nodes.
         */
        void force(mmwis::MISConfig &config, graph_access &G, NodeID v, candidate_list *force_list = NULL,
                   bool with_candidates = true);

        /**
         * Undo all operations currently stored in the operation log.
         *
         * @param G Graph representation
         */
        void unwind(graph_access &G);
};

}

#endif
