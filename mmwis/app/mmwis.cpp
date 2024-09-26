/** 
 * reduction_evomis.cpp
 * Purpose: Main program for the evolutionary algorithm.
 *
 *****************************************************************************/

#include <argtable3.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <random>

#include "timer.h"
#include "hils.h"
#include "mmwis_log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "reduction_evolution.h"
#include "mmwis_config.h"
#include "graph_io.h"
#include "parse_parameters.h"
#include "branch_and_reduce_algorithm.h"

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

std::vector<NodeID> reverse_mapping;
NodeWeight perform_reduction(std::unique_ptr<mmwis::branch_and_reduce_algorithm>& reducer, graph_access& G, graph_access& rG, const mmwis::MISConfig& config) {
	reducer = std::unique_ptr<mmwis::branch_and_reduce_algorithm>(new mmwis::branch_and_reduce_algorithm(G, config));
	reducer->reduce_graph();

	// Retrieve reduced graph
	reverse_mapping = std::vector<NodeID>(G.number_of_nodes(), 0);
	reducer->build_graph_access(rG, reverse_mapping);

	if (!is_IS(rG)) {
		std::cerr << "ERROR: reduced graph is not independent" << std::endl;
		exit(1);
	}

	NodeWeight is_weight = reducer->get_current_is_weight();

	return is_weight;
}


void assign_weights(graph_access& G, const mmwis::MISConfig& mis_config) {
	constexpr NodeWeight MAX_WEIGHT = 200;

	if (mis_config.weight_source == mmwis::MISConfig::Weight_Source::HYBRID) {
		forall_nodes(G, node) {
			G.setNodeWeight(node, (node + 1) % MAX_WEIGHT + 1);
		} endfor
	} else if (mis_config.weight_source == mmwis::MISConfig::Weight_Source::UNIFORM) {
		std::default_random_engine generator(mis_config.seed);
  		std::uniform_int_distribution<NodeWeight> distribution(1,MAX_WEIGHT);

		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	} else if (mis_config.weight_source == mmwis::MISConfig::Weight_Source::GEOMETRIC) {
		std::default_random_engine generator(mis_config.seed);
  		std::binomial_distribution<int> distribution(MAX_WEIGHT / 2);

		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	} else if (mis_config.weight_source == mmwis::MISConfig::Weight_Source::UNIT) {
		forall_nodes(G, node) {
			G.setNodeWeight(node, 1);
		} endfor

	}
}



template<class reducer>
int run(mmwis::MISConfig &mis_config, graph_access &G, NodeWeight weight_offset) {

    mmwis::mmwis_log::instance()->restart_total_timer();
    std::vector<int> forced_vertices(G.number_of_nodes(), -1);

    // Perform the evolutionary algorithm
    std::vector<bool> independent_set(G.number_of_nodes(), false);
    mmwis::reduction_evolution<mmwis::branch_and_reduce_algorithm> evo;
    std::vector<NodeID> best_nodes;
    std::vector<NodeID> worse_nodes;
    bool solved_exactly = false;
    evo.perform_mis_search(mis_config, G, independent_set, best_nodes, worse_nodes, solved_exactly, false);
    mmwis::mmwis_log::instance()->print_results();
    if (mis_config.print_log) mmwis::mmwis_log::instance()->write_log();
    /* if (mis_config.write_graph) graph_io::writeIndependentSet(G, mis_config.output_filename); */
    #ifdef NDEBUG
    if(is_IS(G)) {
        int solution_weight = weight_offset;
        forall_nodes(G, node) {
                if(independent_set[node]) {
                        solution_weight += G.getNodeWeight(node);
                }
        } endfor
        // std::cout << "MWIS_weight " << solution_weight << std::endl; 
    }
    #endif
    return 0;
}
		

int main(int argn, char **argv) {
    mmwis::mmwis_log::instance()->restart_overall_total_timer();
    mmwis::mmwis_log::instance()->print_memetic_title();
    
    mmwis::MISConfig mis_config;
    std::string graph_filepath;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
    if (ret_code) {
        return 0;
    }
    mis_config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of( '/' ) +1);
    mmwis::mmwis_log::instance()->set_config(mis_config);

    // Read the graph
    graph_access G;
    std::string comments;
    graph_io::readGraphWeighted(G, graph_filepath);
    assign_weights(G, mis_config);

    mmwis::mmwis_log::instance()->set_graph(G);

    NodeWeight weight_offset = 0;
    
    // Print setup information
    mmwis::mmwis_log::instance()->print_graph();
    mmwis::mmwis_log::instance()->print_config();

	if (mis_config.perform_reductions) {
        weight_offset = 0;
        run<mmwis::branch_and_reduce_algorithm>(mis_config, G, weight_offset);
        return 0;
	}

}

