/**
 * mmwis_config.h
 * Purpose: Configuration used for the evolutionary maximum independent set algorithms.
 *
 *****************************************************************************/

#ifndef _MMWIS_CONFIG_H_
#define _MMWIS_CONFIG_H_

#include <string>
#include <cctype>
/* #include <algorithm> */

#include "definitions.h"

namespace mmwis {

// Configuration for the calculation of the MIS
struct MISConfig {
    enum Weight_Source {FILE, HYBRID, UNIFORM, GEOMETRIC, UNIT};
    enum Struction_Type {ORIGINAL, MODIFIED, EXTENDED, EXTENDED_REDUCED, NONE};
    enum Key_Type {RANDOM, DEGREE, INCREASE, APPROXIMATE_INCREASE};
    enum Backtrack_Type {IMMEDIATE_TIE_BREAKING, IMMEDIATE_EXCLUDE, END_MIN_KERNEL, NO_BACKTRACK};
    enum Vertex_Selection {degree, weight, degree_over_weight, weight_over_degree, hybrid, solution_participation};
    enum Reduction_Style {initial, time_ordering, weight_ordering, time_and_weight_ordering};
    enum StructionReduction_Style {DENSE, NORMAL};

    // Name of the graph file.
    std::string graph_filename;
    // Name of the output file.
    std::string output_filename;
    // Name of the kernel file.
    std::string kernel_filename;
    // disable simple reduction
    bool disable_fold1=false;
    bool disable_v_shape_min=false;
    bool disable_v_shape_mid=false;
    bool disable_v_shape_max=false;
    bool disable_triangle = false;
    bool disable_basic_se = false;
    bool disable_extended_se = false;
    bool disable_twin= false;
    bool disable_clique = false;
    bool disable_clique_neighborhood= false;
    bool disable_clique_neighborhood_fast= false;
    bool disable_generalized_fold= false;
    bool disable_generalized_neighborhood = false;
    bool disable_critical_set = false;
    bool disable_neighborhood= false;
    bool disable_heavy_set = false;
    // bound for number of nodes in heavy_set neighborhood graph =0 is disabled reduction completely
    int  heavy_set= 20;
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
    // Time limit for the total algorithm
    double time_limit;
    // Time limit for the evolutionary algorithm
    double evo_time_limit;
    // Time limit for the ils 
    double ils_time_limit;
    // Node Threshold when exact solver starts on reduced graph 
    double V_solve_exact;
    // time to take to solve exact when V_solve_exact threshold reached 
    double time_solve_exact;
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
    // Write the kernel into a file
    bool write_kernel;
    // Write the log into the console
    bool console_log;
    // Number of iterations for the ILS.
    unsigned int ils_iterations;
    // Use the max-flow algorithm to fix vertex cover candidates
    bool use_max_flow;
    // Lower bound for the change rate between best individuals
    double best_limit;
    // Fraction of independent set nodes added to the solution
    float fraction;
    // Optimize candidates for ILS.
    bool optimize_candidates;
    // Remove IS nodes from best individual before recursive reduction
    bool extract_best_nodes;
    // Apply all reductions to reduce the graph size
    bool all_reductions;
    // Threshold for reducing the graph
    unsigned int reduction_threshold;
    // Threshold for local search 
    unsigned int local_search_threshold;
    // Check graph sortedness
    bool check_sorted;
    // Use adaptive greedy starting solution
    bool start_greedy_adaptive;
	// Sort free nodes in local search
	bool sort_freenodes;
	// Perform reduction
	bool perform_reductions;
    // perform hils
    bool perform_hils;
    // Random weights or read weights from file
    Weight_Source weight_source;
    // Random weights or read weights from file
    Vertex_Selection vertex_selection;
    // Choose reduction order and amount for given graph type
    Reduction_Style reduction_style;
    // Choose struction reduction order
    StructionReduction_Style struction_reduction_style;

    // for struction
    Struction_Type struction_type;
    unsigned int struction_degree;
    unsigned int set_limit;
    int use_struction_initial_sol;
    double global_blow_up_factor;
    double phase_blow_up_factor;
    double key_reinsert_factor;
    unsigned int phase_blow_ups;
    unsigned int max_unimproving_phases;
    Key_Type key_type;
    Backtrack_Type backtrack_style;
    bool reduce_and_peel;
    bool plain_struction;
    bool disable_blow_up;

    void setBacktrackType(const std::string & back_type) {
        if (strCompare(back_type, "immediate_tie_breaking")) {
            backtrack_style = Backtrack_Type ::IMMEDIATE_TIE_BREAKING;
        } else if (strCompare(back_type, "immediate_exclude")) {
            backtrack_style = Backtrack_Type ::IMMEDIATE_EXCLUDE;
        } else if (strCompare(back_type, "end_min_kernel")) {
            backtrack_style = Backtrack_Type::END_MIN_KERNEL;
        } else if (strCompare(back_type, "no_backtrack")) {
            backtrack_style = Backtrack_Type ::NO_BACKTRACK;
        }
    }

    void setStructionType(const std::string & s_type) {
        if (strCompare(s_type, "original")) {
            struction_type = Struction_Type ::ORIGINAL;
        } else if (strCompare(s_type, "modified")) {
            struction_type = Struction_Type ::MODIFIED;
        } else if (strCompare(s_type, "extended")) {
            struction_type = Struction_Type::EXTENDED;
        } else if (strCompare(s_type, "extended_reduced")) {
            struction_type = Struction_Type ::EXTENDED_REDUCED;
        } else if (strCompare(s_type, "none")) {
            struction_type = Struction_Type ::NONE;
        }
    }

    void setKeyType(const std::string & k_type) {
        if (strCompare(k_type, "random")) {
            key_type = Key_Type ::RANDOM;
        } else if (strCompare(k_type, "degree")) {
            key_type = Key_Type ::DEGREE;
        } else if (strCompare(k_type, "increase")) {
            key_type = Key_Type::INCREASE;
        } else if (strCompare(k_type, "approximate_increase")) {
            key_type = Key_Type ::APPROXIMATE_INCREASE;
        }
    }

    void setReductionStyle(const std::string & redu_style) {
        if (strCompare(redu_style, "time"))                    { reduction_style = Reduction_Style::time_ordering;
        } else if (strCompare(redu_style, "weight"))           { reduction_style = Reduction_Style::weight_ordering;
        } else if (strCompare(redu_style, "time_and_weight"))  { reduction_style = Reduction_Style::time_and_weight_ordering;
        } else { reduction_style = Reduction_Style::initial;}
    }
 
    void setStructionReductionStyle(const std::string & redu_style) {
        if (strCompare(redu_style, "DENSE")) 
                struction_reduction_style = StructionReduction_Style::DENSE;
        else    struction_reduction_style = StructionReduction_Style::NORMAL;
    }

    void setWeightSource(const std::string & source) {
        if (strCompare(source, "file")) {
            weight_source = Weight_Source::FILE;
        } else if (strCompare(source, "hybrid")) {
            weight_source = Weight_Source::HYBRID;
        } else if (strCompare(source, "unit")) {
            weight_source = Weight_Source::UNIT;
        } else if (strCompare(source, "uniform")) {
            weight_source = Weight_Source::UNIFORM;
        } else if (strCompare(source, "geometric")) {
            weight_source = Weight_Source::GEOMETRIC;
        }
    }

    void setVertexSelection(const std::string & source) {
        if (strCompare(source, "degree")) {
            vertex_selection = Vertex_Selection::degree;
        } else if (strCompare(source, "weight")) {
            vertex_selection = Vertex_Selection::weight;
        } else if (strCompare(source, "hybrid")) {
            vertex_selection = Vertex_Selection::hybrid;
        } else if (strCompare(source, "weight_over_degree")) {
            vertex_selection = Vertex_Selection::weight_over_degree;
        } else { // solution_participation
            vertex_selection = Vertex_Selection::solution_participation;
        }
    }

    private:

    bool strCompare(const std::string & str1, const std::string & str2) {
        return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](unsigned char c1, unsigned char c2){ return std::toupper(c1) == std::toupper(c2); });
    }
};

}
#endif
