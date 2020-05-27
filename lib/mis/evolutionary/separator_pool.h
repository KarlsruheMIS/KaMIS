/**
 * separator_pool.h
 * Purpose: Create and maintain a number of node separators and according partitions.
 *
 *****************************************************************************/

#ifndef _SEPARATOR_POOL_H_
#define _SEPARATOR_POOL_H_

#include <vector>

#include "mis_config.h"
#include "definitions.h"
#include "population_mis.h"
#include "kaHIP_interface.h"
#include "data_structure/graph_access.h"

/**
 * Struct for storing a single partitioning.
 */
struct partition {
    int *partition_map;
    int k;
    unsigned int id;
};

/**
 * Struct for storing a node separator partitioning
 * and the according separator size.
 */
struct separator {
    int *separator_map;
    unsigned int separator_size;
    int k;
    unsigned int id;
};

class separator_pool {
    public:
        /**
         * Default Constructor.
         */
        separator_pool();

        /**
         * Default Destructor.
         */
        virtual ~separator_pool();

        /**
         * Initialize the pool.
         * Sets pool size and other parameters.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         */
        void init(MISConfig & config, graph_access & G);

        /**
         * Clear the pool of node separators
         */
        void clear_separators();

        /**
         * Clear the pool of k-separators
         */
        void clear_k_separators();

        /**
         * Clear the pool of partitionings.
         */
        void clear_partitions();

        /**
         * Clear the pool of k-partitionings.
         */
        void clear_k_partitions();

        /**
         * Clear the scores.
         *
         * @param config Config for the pool.
         */
        void clear_scores(MISConfig & config);

        /**
         * Generate a set number of node separators.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         */
        void generate_separators(MISConfig & config, graph_access & G);

        /**
         * Generate a set number of k-separators.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         */
        void generate_k_separators(MISConfig & config, graph_access & G);

        /**
         * Generate a set number of partitionings.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         */
        void generate_partitions(MISConfig & config, graph_access & G);

        /**
         * Generate a set number of k-partitionings.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         */
        void generate_k_partitions(MISConfig & config, graph_access & G);
        
        /**
         * Create a single (k-)separator.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param sep The resulting node separator.
         * @param number_of_blocks The number of blocks for the separator
         */
        void create_separator(MISConfig & config, graph_access & G, separator & sep, int number_of_blocks = 2);

        /**
         * Create a single (k-)partitioning.
         * Uses the KaHIP framework.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param part The resulting partitioning.
         * @param number_of_blocks The number of blocks for the partition.
         */
        void create_partition(MISConfig & config, graph_access & G, partition & part, int number_of_blocks = 2);

        /**
         * Apply the given node separator to the graph.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param sep Node separator to be applied.
         */
        void apply_separator(MISConfig & config, graph_access & G, separator & sep);

        /**
         * Apply the given (k-)partitioning to the graph.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param part Partitioning to be applied.
         */
        void apply_partition(MISConfig & config, graph_access & G, partition & part);

        /**
         * Calculate the scores for all k-partitions and individuals.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param pop The population of the evolutionary algorithm.
         */
        void calculate_partition_scores(MISConfig & config, graph_access & G, population_mis & pop);

        /**
         * Calculate the scores for all k-separators and individuals.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param pop The population of the evolutionary algorithm.
         */
        void calculate_separator_scores(MISConfig & config, graph_access & G, population_mis & pop);

        /**
         * Update all scores for a single individuum.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param ind The individuum to update.
         */
        void update_scores_for_individuum(MISConfig & config, graph_access & G, individuum_mis & ind);

        /**
         * Return a random separator from the pool.
         *
         * @param sep Random node separator.
         */
        void get_random_separator(separator & sep);

        /**
         * Return a random k-separator from the pool.
         *
         * @param sep Random node separator.
         */
        void get_random_k_separator(separator & sep);

        /**
         * Return a random partitioning from the pool.
         *
         * @param part Random partitioning.
         */
        void get_random_partition(partition & part);

        /**
         * Return a random k-partitioning from the pool.
         *
         * @param part Random partitioning.
         */
        void get_random_k_partition(partition & part);

        /**
         * Gather the best candidates for a given partition.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param pop The population used for getting individuals.
         * @param candidates The resulting candidate list.
         * @param part The partition for which the candidates should be found.
         */
        void get_best_candidates(MISConfig & config, graph_access & G, population_mis & pop, std::vector<individuum_mis> & candidates, partition & part);

        /**
         * Gather the best candidates for a given separator.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param pop The population used for getting individuals.
         * @param candidates The resulting candidate list.
         * @param sep The separator for which the candidates should be found.
         */
        void get_best_candidates(MISConfig & config, graph_access & G, population_mis & pop, std::vector<individuum_mis> & candidates, separator & sep);

        /** 
         * Renew the separator pool.
         * Removes and refreshes partitions/separators.
         * Calculates new scores if necessary.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param init Clear old partitions/separators?
         * @param scores Calculate scores?
         * @param pop Population of the evolutionary algorithm.
         */
        void renew_pool(MISConfig & config, graph_access & G, bool init, bool scores, population_mis & pop);

        /**
         * Print the sizes of all separators in the pool.
         */
        void print_separator_sizes();

        /**
         * Print the current score cache.
         *
         * @param config Config for the pool.
         */
        void print_scores(MISConfig & config);

        /**
         * Calculate the scores for all partitions/separators and the given population.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param pop Population of individuals.
         */
        void calculate_scores(MISConfig & config, graph_access & G, population_mis & pop);

        /**
         * Clear the stored graph data structures used for KaHIP calls.
         */
        void clear_graph();

    private:
        // Array for storing the node separators.
        std::vector<separator> internal_separators;
        // Array for storing the partitionings.
        std::vector<partition> internal_partitions;
        // Array for storing the k-partitionings.
        std::vector<partition> internal_k_partitions;
        // Array for storing the k-separators.
        std::vector<separator> internal_k_separators;

        // Number of node separators in the pool.
        unsigned int separators_size;
        // Number of k-separators in the pool.
        unsigned int k_separators_size;
        // Number of partitionings in the pool.
        unsigned int partitions_size;
        // Number of k-partitionings in the pool.
        unsigned int k_partitions_size;

        // Graph data structures used for the KaHIP-library calls.
        int *xadj;
        int *adjncy;
        // Score cache
        std::vector<unsigned int> scores;

        /**
         * Grow blocks for an applied separator.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param k Number of blocks.
         */
        void separator_bfs(MISConfig & config, graph_access & G, unsigned int k);

        /**
         * Calculate the score for a given individuum and partition/separator id.
         *
         * @param config Config for the pool.
         * @param G Graph representation.
         * @param ind Individuum for which the score should be calculated.
         * @param id Id of the separator/partition.
         */
        void calculate_score(MISConfig & config, graph_access & G, individuum_mis & ind, unsigned int id);
};

#endif

