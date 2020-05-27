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
    // Directory containing partitions
    std::string partition_directory;
    // Name of the output file.
    std::string output_filename;
    // Number of repititions for benchmark
    int num_reps;
    // Write the log into a file
    bool print_log;
    // Write the inpendent set into a file
    bool write_graph;
    // Write the log into the console
    bool console_log;
    // Apply all reductions to reduce the graph size
    bool all_reductions;
    // Check graph sortedness
    bool check_sorted;
};

#endif
