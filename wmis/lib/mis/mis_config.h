/**
 * mis_config.h
 * Purpose: Configuration used for the evolutionary maximum independent set algorithms.
 *
 *****************************************************************************/

#ifndef _MIS_CONFIG_H_
#define _MIS_CONFIG_H_

#include <string>
#include <cctype>
#include <algorithm>

#include "definitions.h"

// Configuration for the calculation of the MIS
struct MISConfig {
    enum Weight_Source {FILE, HYBRID, UNIFORM, GEOMETRIC};
    enum Reduction_Style {NORMAL, DENSE};

    // Name of the graph file.
    std::string graph_filename;
    // Name of the output file.
    std::string output_filename;
    // Seed for the RNG.
    int seed;
    // Time limit for the evolutionary algorithm
    double time_limit;
    // Write the log into a file
    bool print_log;
    // Write the inpendent set into a file
    bool write_graph;
    // Write the log into the console
    bool console_log;
    // Number of iterations for the ILS.
    unsigned int ils_iterations;
    // Lower bound for the change rate between best individuals
    double best_limit;
    // Number of candidates for forced insertion.
    unsigned int force_cand;
    // Optimize candidates for ILS.
    bool optimize_candidates;
    // Check graph sortedness
    bool check_sorted;
	// Sort free nodes in local search
	bool sort_freenodes;
	// Perform reduction
	bool perform_reductions;
    // Random weights or read weights from file
    Weight_Source weight_source;
    // Choose reduction order and amount for given graph type
    Reduction_Style reduction_style;


    void setReductionStyle(const std::string & redu_style) {
        if (strCompare(redu_style, "normal") || strCompare(redu_style, "sparse")) {
            reduction_style = Reduction_Style::NORMAL;
        } else if (strCompare(redu_style, "dense") || strCompare(redu_style, "osm")) {
            reduction_style = Reduction_Style::DENSE;
        }
    }

    void setWeightSource(const std::string & source) {
        if (strCompare(source, "file")) {
            weight_source = Weight_Source::FILE;
        } else if (strCompare(source, "hybrid")) {
            weight_source = Weight_Source::HYBRID;
        } else if (strCompare(source, "uniform")) {
            weight_source = Weight_Source::UNIFORM;
        } else if (strCompare(source, "geometric")) {
            weight_source = Weight_Source::GEOMETRIC;
        }
    }


    private:

    bool strCompare(const std::string & str1, const std::string & str2) {
        return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](unsigned char c1, unsigned char c2){ return std::toupper(c1) == std::toupper(c2); });
    }
};

#endif
