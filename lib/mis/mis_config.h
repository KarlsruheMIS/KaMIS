/**
 * mis_config.h
 * Purpose: Configuration used for the evolutionary maximum independent set algorithms.
 *
 *****************************************************************************/

#ifndef _MIS_CONFIG_H_
#define _MIS_CONFIG_H_

#include <string>

#include "definitions.h"

// Configuration for the calculation of the MIS
struct MISConfig {
    // Name of the graph file.
    std::string graph_filename;
    // Name of the output file.
    std::string output_filename;
    // Seed for the RNG.
    int seed;
    // Size of the population used.
    unsigned int population_size;
    // Imbalance. Used for the KaHIP-interface calls.
    double imbalance;
    // Mode for the KaHIP-framework.
    unsigned int kahip_mode;
    // Use full kernelization (true) or FastKer (false)
    bool fullKernelization;
    // Time limit for the evolutionary algorithm
    double time_limit;
    // Number of repetitions in each round.
    unsigned int repetitions;
    // Insert a solution if no new solution has been inserted
    // for this amount of operations.
    unsigned int insert_threshold;
    // Update the pool of node separators.
    unsigned int pool_threshold;
    // Factor for timing based renewal of the separator.
    // Time after building > Factor * Time taken for building triggers renewal.
    double pool_renewal_factor;
    // Should the imbalance be randomized?
    bool randomize_imbalance;
    // Diversify the initial solution?
    bool diversify;
    // Tournament or random selection of individuals?
    bool enable_tournament_selection;
    // Use the vertex cover approach for the multiway combine operator.
    bool use_multiway_vc;
    // Number of blocks used in multiway operators.
    unsigned int multiway_blocks;
    // Number of individuals in a tournament
    unsigned int tournament_size;
    // Percentage for the mutation of a solution.
    int flip_coin;
    // Force parameter for the mutation.
    unsigned int force_k;
    // Number of candidates for forced insertion.
    unsigned int force_cand;
    // Number of initial node separators to be constructed.
    unsigned int number_of_separators;
    // Number of initial partitions to be constructed.
    unsigned int number_of_partitions;
    // Number of initial k-separators to be constructed.
    unsigned int number_of_k_separators;
    // Number of initial k-partitions to be constructed.
    unsigned int number_of_k_partitions;
    // Print result of each repetition.
    bool print_repetition;
    // Print the population after each round.
    bool print_population;
    // Write the log into a file
    bool print_log;
    // Write the inpendent set into a file
    bool write_graph;
    // Write the log into the console
    bool console_log;
    // Number of iterations for the ILS.
    unsigned int ils_iterations;
    // Use the Hopcroft-Karp algorithm to fix vertex cover candidates
    bool use_hopcroft;
    // Lower bound for the change rate between best individuals
    double best_limit;
    // Fraction of independent set nodes that should be removed before reduction
    double remove_fraction;
    // Optimize candidates for ILS.
    bool optimize_candidates;
    // Remove IS nodes from best individual before recursive reduction
    bool extract_best_nodes;
    // Apply all reductions to reduce the graph size
    bool all_reductions;
    // Threshold for reducing the graph
    unsigned int reduction_threshold;
    // Check graph sortedness
    bool check_sorted;
    // Use adaptive greedy starting solution
    bool start_greedy_adaptive;
};

#endif
