 /******************************************************************************
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
  *****************************************************************************/

#include "full_reductions.h"
#include <vector>
#include <memory>
#include "parallel_reductions.h"
#include <omp.h>
#include <algorithm>
#include <fstream>
#include "kaHIP_interface.h"

full_reductions::full_reductions(std::vector<std::vector<int>> &_adj, int _N)
: adj(_adj) 
{
}

void full_reductions::reduce_graph() {
    // std::cout << "CALL REDUCE" << std::endl;
    int size_after_force = adj.size();
    std::vector<char> removed(adj.size(), 0);
    // #pragma omp parallel for
    for(unsigned int vertex : forced_vertices) {
        removed[vertex] = 1;
        size_after_force--;
        for(int neighbor: adj[vertex]) {
            if(removed[neighbor] == 0) {
                size_after_force--;
            }
            removed[neighbor] = 1;
        }
    }
    // std::cout << "marked removed vertices" << std::endl;
    after_force_mapping.resize(adj.size());
    after_force_reverse_mapping.resize(size_after_force);
    int mapping_counter = 0;
    for(int vertex = 0; vertex < adj.size(); ++vertex) {
        if(removed[vertex] == 0) {
            after_force_mapping[vertex] = mapping_counter;
            after_force_reverse_mapping[mapping_counter] = vertex;
            mapping_counter++;
        }
    }
    if(mapping_counter != size_after_force) {
        std::cout << "Missmatching sizes!" << std::endl;
        exit(1);
    }
    // std::cout << "finished mapping" << std::endl;
    std::vector<ui> pstart(size_after_force + 1);
    std::vector<ui> edges = std::vector<ui>();

    // pstart[0] = 0;
    int adjacency_counter = 0;
    for(ui i = 0; i < size_after_force; ++i) {
        pstart[i] = adjacency_counter;
        for (int neighbor: adj[after_force_reverse_mapping[i]]) if(removed[neighbor] == 0){
                    edges.push_back(after_force_mapping[neighbor]);
                    adjacency_counter++;
            }

            pstart[i + 1] = edges.size();
    }
    pstart[size_after_force] = adjacency_counter;
    // std::cout << "finished building edgelist" << std::endl;

    LineartimeKernelizer = LinearTime::Graph();
    LineartimeKernelizer.inputEdgeList(pstart, edges);
    // std::cout << "finished adding edgelist to lineartime" << std::endl;
    // LineartimeKernelizer.inputGraphAdjList(adj);
    std::vector<std::vector<int>> LineartimeKernel;
    linearTimeOffset = LineartimeKernelizer.degree_two_kernal_and_remove_max_degree_without_contraction(LineartimeKernel);

    std::vector<int> partitions = std::vector<int>(LineartimeKernel.size());
    // int numPartitions = 32;//*std::max_element(partitions.begin(), partitions.end());
    // double imbalance = 0.03;
    // int edgecut = 0;
    // int m = 0;
    // for(auto neighborarray : LineartimeKernel) {
    //     m += neighborarray.size();
    // }
    // std::vector<int> xadj(LineartimeKernel.size() + 1);
    // std::vector<int> adjncy(m);
    // unsigned int adjncy_counter = 0;
    // for (unsigned int i = 0; i < LineartimeKernel.size(); ++i) {
    //     xadj[i] = adjncy_counter;
    //     for (int const neighbor : LineartimeKernel[i]){
    //             adjncy[adjncy_counter++] = neighbor;
    //         }
    //     std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
    // }
    // xadj[LineartimeKernel.size()] = adjncy_counter;
    // std::cout << "Call kaffpa" << std::endl;
    // int LinearTimeKernelSize = LineartimeKernel.size();
    // kaffpa( &LinearTimeKernelSize, 
    //         NULL,
    //         xadj.data(),
    //         NULL,
    //         adjncy.data(),
    //         &numPartitions, 
    //         &imbalance,
    //         false,
    //         1337,
    //         0,
    //         &edgecut,
    //         partitions.data());
	//std::cout << "Creating object" << std::endl;
	parallel_reducer = std::unique_ptr<parallel_reductions>(new parallel_reductions(LineartimeKernel, partitions));
	//std::cout << "Finished creating object" << std::endl;
	// std::cout << "Before call to parallel reduce_graph" << std::endl;
	parallel_reducer->reduce_graph_parallel(forced_vertices);
	// std::cout << "After call to parallel reduce_graph" << std::endl;
	//std::cout << "Kernel size after parallel run: " << parallel_reducer->size() << std::endl;
	// std::cout << "Before call to sequential reduce_graph" << std::endl;
    // parallel_reducers.back()->reduce_graph_sequential();
	// std::cout << "After call to sequential reduce_graph" << std::endl;
	// std::cout << "Kernel size after sequential run: " << parallel_reducers.back()->size() << std::endl;

    // std::vector<bool> independent_set(adj.size());
    // parallel_reducer->ExtendPartialSolution(independent_set);
    // exit(1);
}


size_t full_reductions::get_current_is_size_with_folds() {
    // std::cout << "get current size" << std::endl;
    // std::cout << parallel_reducer->get_current_is_size_with_folds() << std::endl;
    // std::cout << linearTimeOffset << std::endl;
    // std::cout << forced_vertices.size() << std::endl;
	return parallel_reducer->get_current_is_size_with_folds() + linearTimeOffset + forced_vertices.size();
}

size_t full_reductions::number_of_nodes_remaining() {
	return parallel_reducer->size();
}

void full_reductions::force_into_independent_set(std::vector<NodeID> &nodes_to_force){
    // parallel_reducer->force_into_independent_set(nodes_to_force);
    forced_vertices = std::vector<unsigned int>(nodes_to_force);
}

void full_reductions::convert_adj_lists(graph_access &G, std::vector<NodeID> &reverse_mapping) {
    parallel_reducer->getGraphAccess(G, reverse_mapping);
}
void full_reductions::extend_finer_is(std::vector<bool> &independent_set) {
    // for(int vertex = 0; vertex < adj.size(); ++vertex) {
    //     if(independent_set[vertex]) {
    //         for(int neighbor : adj[vertex]) {
    //             if(independent_set[neighbor]) {
    //                 std::cout << "Before: Not an independent set!" << std::endl;
    //                 exit(1);
    //             }
    //         }
    //     }
    // }
    // std::cout << "Before: Valid independent set!" << std::endl;
    parallel_reducer->ExtendPartialSolution(independent_set);
    LineartimeKernelizer.UndoReductions(independent_set);
    std::vector<bool> is_copy(independent_set);
    independent_set.resize(adj.size());
    for(int i = 0; i < adj.size(); ++i) {
        independent_set[i] = is_copy[after_force_mapping[i]];
    }
    for(int vertex : forced_vertices) {
        independent_set[vertex] = true;
        for(int neighbor: adj[vertex]) {
            independent_set[neighbor] = false;
        }
    }
    // for(int vertex = 0; vertex < adj.size(); ++vertex) {
    //     if(independent_set[vertex]) {
    //         for(int neighbor : adj[vertex]) {
    //             if(independent_set[neighbor]) {
    //                 std::cout << "After: Not an independent set!" << std::endl;
    //                 exit(1);
    //             }
    //         }
    //     }
    // }
    // std::cout << "After: Valid independent set!" << std::endl;
}
