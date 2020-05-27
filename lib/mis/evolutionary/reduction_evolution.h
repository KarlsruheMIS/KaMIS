/**
 * reduction_evolution.h
 * Purpose: Recursive evolutionary framework including reduction techniques.
 *
 *****************************************************************************/

#ifndef _REDUCTION_EVOLUTION_H_
#define _REDUCTION_EVOLUTION_H_

#include "timer.h"
#include "mis_config.h"
#include "population_mis.h"
#include "separator_pool.h"
#include "data_structure/graph_access.h"
#include "mis/kernel/branch_and_reduce_algorithm.h"
#include "../kernel/ParFastKer/fast_reductions/src/full_reductions.h"

// system includes
#include <memory> // unique_ptr

template <typename reducer>
class reduction_evolution {
    public:
        /**
         * Default Constructor.
         */
        reduction_evolution ();

        /**
         * Default Destructor.
         */
        virtual ~reduction_evolution ();

        /**
         * Perform the algorithm until the time limit is reached.
         * Uses reduction techniques in combination with an 
         * evolutionary approach to generate independent sets.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         * @param independent_set Current independent set.
         * @param best_nodes Best nodes from the upper level that should be removed.
         * @param recursive Recursive call.
         * @param is_size Current independent set size.
         */
        unsigned int perform_mis_search(MISConfig & mis_config, 
                                graph_access & G, 
                                std::vector<bool> & independent_set,
                                std::vector<NodeID> & best_nodes,
                                bool recursive = false, 
                                unsigned int is_size = 0);

        /**
         * Returns the best maximum independent set found so far.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         * @param out Individuum containing the best solution.
         * @return Size of the solution.
         */
        unsigned int collect_best_mis(MISConfig & mis_config, 
                                      graph_access & G, 
                                      individuum_mis & out);

    private:
        // Population used in the framework
        population_mis island;
        // Pool of node separators/partitions
        separator_pool *pool;

        unsigned int finer_max; 
        unsigned int coarser_max; 
        unsigned int is_base; 
        unsigned int pool_counter; 
        unsigned int reduction_counter; 

        /**
         * Performs a single round of the evolutionary algorithm.
         * Each round consists of a number of repetitions.
         * During each repetitions individuals are combined to generate a new offspring.
         * This offspring can then mutate and enter the population.
         * The details for the mutate and combine operators can be found
         * in the population description.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         */
        void perform_local_mis(MISConfig & mis_config, 
                               graph_access & G);

        /**
         * Set random seed and initialize timer.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         */
        void init(MISConfig & mis_config, 
                  graph_access & G);

        /**
         * Apply reductions to the given graph.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         * @param reduced Resulting reduced graph representation.
         * @param is_nodes Independent set nodes that were removed.
         * @param other_nodes Other nodes that were removed.
         * @param reverse_mapping Mapping from reduced to original graph.
         * @param best_nodes Set of nodes to extract before reductions.
         * @param recursive Recursive or initial call.
         * @param full_reducer Reducer for the complete set of reductions.
         *
         * @return size of the Resulting independent set.
         */
        unsigned int reduce(MISConfig & mis_config, 
                            graph_access & G, 
                            graph_access & reduced, 
                            std::vector<NodeID> & is_nodes,
                            std::vector<NodeID> & other_nodes,
                            std::vector<NodeID> & reverse_mapping,
                            std::vector<NodeID> & best_nodes,
                            bool recursive,
                            std::unique_ptr<reducer> & full_reducer);

        /**
         * Revert a set of reductions.
         * Adds removed independent set nodes into the original graph.
         *
         * @param mis_config Config for the framework.
         * @param independent_set Current best independent set.
         * @param reductions Set of applied reductions.
         * @param full_reducer Reducer used for creation of the reduced graph.
         *
         */
        void add_reductions(MISConfig & mis_config, 
                            std::vector<bool> & independent_set, 
                            std::unique_ptr<reducer> & full_reducer);

        /**
         * Create an initial population.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         */
        void fill_population(MISConfig & mis_config, 
                             graph_access & G);

        /**
         * Create a partition/separator pool.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         */
        void init_pool(MISConfig & mis_config, 
                       graph_access & G);

        void extract_nodes(MISConfig & mis_config, 
                           graph_access & G, 
                           std::vector<NodeID> & is_nodes,
                           std::vector<NodeID> & other_nodes,
                           std::unique_ptr<reducer> &full_reducer);

        /**
         * Extract a set of independent set nodes and 
         * their neighborhood from the graph.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         * @param nodes Nodes to be forcefully removed.
         */
        void get_extract_nodes(MISConfig & mis_config, 
                               graph_access & G, 
                               std::vector<NodeID> & nodes);

        /**
         * Map the independent set from the coarse to the finer level.
         * Mapping is performed according to the independent set sizes
         * on the different levels.
         *
         * @param mis_config Config for the framework.
         * @param coarser Reduced graph representation.
         * @param independent_set Current independent set.
         * @param coarser_is Independent set of the reduced graph.
         * @param reverse_mapping Mapping from reduced to original graph.
         */
        void perform_reverse_mapping(MISConfig & mis_config, 
                                     graph_access & coarser, 
                                     std::vector<bool> & independent_set,
                                     std::vector<bool> & coarser_is,
                                     std::vector<NodeID> & reverse_mapping);

        /**
         * Apply ILS to the current best independent set.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation of the current level.
         * @param independent_set Current independent set.
         */
        void build_final_solution(MISConfig & mis_config, 
                                  graph_access & G, 
                                  std::vector<bool> & independent_set);

        /**
         * Set number of local iterations based on the graph size.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         */
        void set_local_iterations(MISConfig & mis_config, 
                                  graph_access & G);

        /**
         * Performs a single round of the evolutionary algorithm.
         * Each round consists of a number of repetitions.
         * During each repetitions individuals are combined to generate a new offspring.
         * This offspring can then mutate and enter the population.
         * The details for the mutate and combine operators can be found
         * in the population description.
         *
         * @param mis_config Config for the framework.
         * @param G Graph representation.
         * @param independent_set Outputs the best independent set found.
         */
        void perform_evolutionary(MISConfig & mis_config,
                                  graph_access & G, 
                                  std::vector<bool> & independent_set);
};

#endif
