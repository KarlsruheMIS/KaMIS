/**
 * mis_config.h
 * Purpose: Configuration used for the evolutionary maximum independent set algorithms.
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
