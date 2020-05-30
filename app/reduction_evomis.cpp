/** 
 * reduction_evomis.cpp
 * Purpose: Main program for the evolutionary algorithm.
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <argtable3.h>

#include "timer.h"
#include "ils/ils.h"
#include "ils/local_search.h"
#include "mis_log.h"
#include "graph_io.h"
#include "reduction_evolution.h"
#include "mis_config.h"
#include "greedy_mis.h"
#include "parse_parameters.h"
#include "data_structure/graph_access.h"
#include "data_structure/mis_permutation.h"
#include "mis/kernel/ParFastKer/fast_reductions/src/full_reductions.h"

template<class reducer>
int run(MISConfig &mis_config, graph_access &G) {

    // Perform the evolutionary algorithm
    std::vector<bool> independent_set(G.number_of_nodes(), false);
    reduction_evolution<reducer> evo;
    std::vector<NodeID> best_nodes;
    evo.perform_mis_search(mis_config, G, independent_set, best_nodes);

    mis_log::instance()->print_results();
    if (mis_config.print_log) mis_log::instance()->write_log();
    if (mis_config.write_graph) graph_io::writeIndependentSet(G, mis_config.output_filename);

    std::cout <<  "checking solution ..."  << std::endl;
    int counter = 0;
    forall_nodes(G, node) {
            if( independent_set[node] ) {
                    counter++;
                    forall_out_edges(G, e, node) {
                            NodeID target = G.getEdgeTarget(e);
                            if(independent_set[target]) {
                                std::cout <<  "not an independent set!"  << std::endl;
                                exit(1);
                            }
                    } endfor
            }
    } endfor
    std::cout <<  "done ..."  << std::endl;
    std::cout << "Independent set has size " << counter << std::endl;
    return 0;
}
		

int main(int argn, char **argv) {
    mis_log::instance()->restart_total_timer();
    mis_log::instance()->print_title();
    
    MISConfig mis_config;
    std::string graph_filepath;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
    if (ret_code) {
        return 0;
    }
    mis_config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of( '/' ) +1);
    mis_log::instance()->set_config(mis_config);

    // Read the graph
    graph_access G;
    graph_io::readGraphWeighted(G, graph_filepath);
    mis_log::instance()->set_graph(G);
    
    // Print setup information
    mis_log::instance()->print_graph();
    mis_log::instance()->print_config();

    if(mis_config.fullKernelization) {
    	return run<branch_and_reduce_algorithm>(mis_config, G);
    } 
    else {
	// Might be a bit confusingly named, but this is FastKer
    	return run<full_reductions>(mis_config, G);
    }
}

