/**
 * reduction_evomis.cpp
 * Purpose: Main program for the evolutionary algorithm.
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

#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <argtable3.h>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>

#include "timer.h"
#include "mis_log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "mis_config.h"
#include "parse_parameters.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "strongly_connected_components.h"

bool is_IS(graph_access& G) {
    forall_nodes(G, node)
        if (G.getPartitionIndex(node) == 1)
            forall_out_edges(G, edge, node)
                NodeID neighbor = G.getEdgeTarget(edge);
                if (G.getPartitionIndex(neighbor) == 1) {
                    return false;
                }
            endfor

     endfor

    return true;
}

int main(int argn, char **argv) {
    struction::cout_handler::disable_cout();
    mis_log::instance()->restart_total_timer();

    ::mmwis::MISConfig mis_config;
    std::string graph_filepath;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
    if (ret_code)
        return 0;

    mis_config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string path = graph_filepath.substr(0,graph_filepath.find_last_of('/'));
    std::string name = mis_config.graph_filename.substr(0,mis_config.graph_filename.find_last_of('.'));

    // Read the graph
    graph_access G;
    graph_io::readGraphWeighted(G, graph_filepath);

    struction::branch_and_reduce_algorithm reducer(G, mis_config);

    graph_access &g = reducer.kernelize();

    mis_config.phase_blow_ups = 1;
    mis_config.struction_type = ::mmwis::MISConfig::Struction_Type::EXTENDED;

    struction::branch_and_reduce_algorithm tester(g, mis_config);
    tester.run_branch_reduce();
    struction::cout_handler::enable_cout();
    std::cout << name << "," << (tester.get_current_is_weight() + reducer.get_current_is_weight()) << std::endl;

    return 0;
}
