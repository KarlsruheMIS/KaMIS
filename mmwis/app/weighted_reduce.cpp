/**
 * Purpose: Main program weighted kernelization
 *
 *****************************************************************************/

#include <cstddef>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <argtable3.h>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>

#include "definitions.h"
#include "timer.h"
#include "ils/ils.h"
#include "ils/local_search.h"
#include "mis_log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "mmwis_config.h"
#include "greedy_mis.h"
#include "parse_parameters.h"
#include "branch_and_reduce_algorithm.h"


void initial_is(graph_access& G) {
	std::vector<NodeID> nodes(G.number_of_nodes());
	for (size_t i = 0; i < nodes.size(); i++) {
		nodes[i] = i;
	}

	// sort in descending order by node weights
	std::sort(nodes.begin(), nodes.end(), [&G](NodeID lhs, NodeID rhs) {
		return G.getNodeWeight(lhs) > G.getNodeWeight(rhs);
	});

	for (NodeID n : nodes) {
		bool free_node = true;

		forall_out_edges(G, edge, n) {
			NodeID neighbor = G.getEdgeTarget(edge);
			if (G.getPartitionIndex(neighbor) == 1) {
				free_node = false;
				break;
			}
		} endfor

		if (free_node) G.setPartitionIndex(n, 1);
	}
}

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

NodeWeight perform_reduction(std::unique_ptr<branch_and_reduce_algorithm>& reducer, graph_access& G, graph_access& rG, const MISConfig& config) {
	reducer = std::unique_ptr<branch_and_reduce_algorithm>(new branch_and_reduce_algorithm(G, config));
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


void assign_weights(graph_access& G, const MISConfig& mis_config) {
	constexpr NodeWeight MAX_WEIGHT = 200;

	if (mis_config.weight_source == MISConfig::Weight_Source::HYBRID) {
		forall_nodes(G, node) {
			G.setNodeWeight(node, (node + 1) % MAX_WEIGHT + 1);
		} endfor
	} else if (mis_config.weight_source == MISConfig::Weight_Source::UNIFORM) {
		std::default_random_engine generator(mis_config.seed);
  		std::uniform_int_distribution<NodeWeight> distribution(1,MAX_WEIGHT);

		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	} else if (mis_config.weight_source == MISConfig::Weight_Source::GEOMETRIC) {
		std::default_random_engine generator(mis_config.seed);
  		std::binomial_distribution<int> distribution(MAX_WEIGHT / 2);

		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	
	} else if (mis_config.weight_source == MISConfig::Weight_Source::UNIT) {
		forall_nodes(G, node) {
			G.setNodeWeight(node, 1);
		} endfor
	}
}


int main(int argn, char **argv) {
    mis_log::instance()->restart_total_timer();
    //mis_log::instance()->print_title();

    MISConfig mis_config;
    std::string graph_filepath;

    // Parse the command line parameters;
    int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
    if (ret_code) {
        return 0;
    }

    mis_config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of( '/' ) +1);
    mis_log::instance()->set_config(mis_config);
    std::string name = mis_config.graph_filename.substr(0, mis_config.graph_filename.find_last_of('.'));
    std::string path = graph_filepath.substr(0,graph_filepath.find_last_of('/'));

    // Read the graph
    graph_access G;
    graph_io::readGraphWeighted(G, graph_filepath);
	assign_weights(G, mis_config);

    mis_log::instance()->set_graph(G);

   NodeWeight weight_offset = 0;
   std::unique_ptr<branch_and_reduce_algorithm> reducer;

	//if (mis_config.write_graph) {
		//// just reduce the graph and write it into a file
		//graph_access rG;

		//auto start = std::chrono::system_clock::now();
		//weight_offset = perform_reduction(reducer, G, rG, mis_config);
		//auto end = std::chrono::system_clock::now();

		//std::chrono::duration<float> reduction_time = end - start;

		//std::ofstream output_reduced(mis_config.output_filename);

		//output_reduced << "%reduction_time " << reduction_time.count() << "\n";
		//output_reduced << "%reduction_offset " << weight_offset << "\n";

		//graph_io::writeGraphNodeWeighted(rG, output_reduced);
		//return 0;
	//}

	//std::cout << "%nodes " << G.number_of_nodes() << std::endl;

		// recude graph and run local search
		graph_access rG;

		auto start = std::chrono::system_clock::now();
		weight_offset = perform_reduction(reducer, G, rG, mis_config);
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<float> reduction_time = end - start;


		std::cout << "graph " << name << "\n";
		std::cout << "number of nodes: " << G.number_of_nodes() << "\n";
		std::cout << "number of edges: " << G.number_of_edges() << "\n";
		std::cout << "reduced number of nodes " << rG.number_of_nodes() << "\n";
		std::cout << "reduced number of edges " << rG.number_of_edges() << "\n";
		std::cout << "%reduction_time " << reduction_time.count() << "\n";
		std::cout << "%reduction_offset " << weight_offset << std::endl;
		std::cout << "MIS_weight " << weight_offset << std::endl;


		reducer->reverse_reduction(G, rG, reverse_mapping);

		if (!is_IS(G)) {
			std::cerr << "ERROR: graph after inverse reduction is not independent" << std::endl;
            exit(1);
		}
		NodeWeight is_weight = 0;



   forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 1) {
					is_weight += G.getNodeWeight(node);
				}
			} endfor

		std::cout << "MIS_weight_offset " << is_weight << std::endl;
        std::cout << path << " " << name << std::endl;
        std::cout << "used ordering number: "<< mis_config.reduction_style<<  std::endl;
        if(mis_config.write_kernel) {
            graph_io::writeGraphWeighted(rG, mis_config.kernel_filename);
        }

    return 0;
}
