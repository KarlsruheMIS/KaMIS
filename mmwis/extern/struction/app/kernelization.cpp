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
#include "mmwis_config.h"
#include "parse_parameters.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "strongly_connected_components.h"

bool is_IS(graph_access& G) {
    forall_nodes(G, node) {
                if (G.getPartitionIndex(node) == 1) {
                    forall_out_edges(G, edge, node) {
                                NodeID neighbor = G.getEdgeTarget(edge);
                                if (G.getPartitionIndex(neighbor) == 1) {
                                    return false;
                                }
                            } endfor
                }
            } endfor

    return true;
}

void assign_weights(graph_access& G, const ::mmwis::MISConfig& mis_config) {
    constexpr NodeWeight MAX_WEIGHT = 200;

    if (mis_config.weight_source == ::mmwis::MISConfig::Weight_Source::HYBRID) {
        forall_nodes(G, node) {
                    G.setNodeWeight(node, (node + 1) % MAX_WEIGHT + 1);
                } endfor
    } else if (mis_config.weight_source == ::mmwis::MISConfig::Weight_Source::UNIFORM) {
        std::default_random_engine generator(mis_config.seed);
        std::uniform_int_distribution<NodeWeight> distribution(1,MAX_WEIGHT);

        forall_nodes(G, node) {
                    G.setNodeWeight(node, distribution(generator));
                } endfor
    } else if (mis_config.weight_source == ::mmwis::MISConfig::Weight_Source::GEOMETRIC) {
        std::default_random_engine generator(mis_config.seed);
        std::binomial_distribution<int> distribution(MAX_WEIGHT / 2);

        forall_nodes(G, node) {
                    G.setNodeWeight(node, distribution(generator));
                } endfor
    }
}


int main(int argn, char **argv) {
    mis_log::instance()->restart_total_timer();
    //mis_log::instance()->print_title();

    ::mmwis::MISConfig mis_config;
    std::string graph_filepath;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
    if (ret_code) {
        return 0;
    }

    mis_config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    mis_log::instance()->set_config(mis_config);

    std::string path = graph_filepath.substr(0,graph_filepath.find_last_of('/'));
    std::string name = mis_config.graph_filename.substr(0,mis_config.graph_filename.find_last_of('.'));
    std::cout << name << ",";

    // Read the graph
    graph_access G;
    graph_io::readGraphWeighted(G, graph_filepath);
    assign_weights(G, mis_config);

    std::string path_and_file = graph_filepath.substr(0,graph_filepath.find_last_of('-'));

    if (mis_config.weight_source == ::mmwis::MISConfig::Weight_Source::UNIFORM) {
        graph_io::writeGraphWeighted(G, path_and_file + "-uniform.graph");
    }

    mis_log::instance()->set_graph(G);

    struction::branch_and_reduce_algorithm reducer(G, mis_config);
    timer t;

    struction::cout_handler::enable_cout();
    graph_access &g = reducer.kernelize();
    double time = t.elapsed();

    std::vector<int> comp_map(g.number_of_nodes(), 0);
    size_t comp_count = strongly_connected_components().strong_components(g, comp_map);
    std::vector<size_t> comp_size(comp_count, 0);
    for (auto comp_id : comp_map)
        comp_size[comp_id]++;
    size_t max_component = comp_count != 0 ? *std::max_element(comp_size.begin(), comp_size.end()) : 0;

    struction::cout_handler::enable_cout();
    std::cout << G.number_of_nodes() << "," << G.number_of_edges() / 2 << ",";
    std::cout << g.number_of_nodes() << "," << g.number_of_edges() / 2 << ",";
    std::cout << comp_count << "," << max_component << ",";
    std::cout << time << "," << reducer.get_current_is_weight() << std::endl;

    if (mis_config.write_kernel)
        graph_io::writeGraphWeighted(g, mis_config.kernel_filename);
    return 0;
}
