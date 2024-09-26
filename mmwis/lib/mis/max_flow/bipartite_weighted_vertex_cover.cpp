/**
 * bipartite_weighted_vertex_cover.cpp
 * Purpose: Create a minimal vertex cover for a given graph.
 *****************************************************************************/

#include "bipartite_weighted_vertex_cover.h"
#include "definitions.h"
#include "flow_graph.h"
#include "push_relabel.h"
#include <limits>
#include <queue>

#include "graph_access.h"
#include "macros_assertions.h"
#include "random_functions.h"

bipartite_weighted_vertex_cover::bipartite_weighted_vertex_cover() {
    
}

bipartite_weighted_vertex_cover::~bipartite_weighted_vertex_cover() {

}


void bipartite_weighted_vertex_cover::max_flow_cover(graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, std::vector<NodeID> & vertex_cover) {

    // Create the bipartite graph and flow_network 
    graph_access bipartite;
    std::vector<NodeID> mapping;
    create_bipartite(G, lhs, rhs, bipartite, mapping);

    flow_graph flow_network;
    push_relabel max_flow_algorithm;
    NodeID source;
    NodeID sink;
    create_flow_network(bipartite, flow_network, source, sink);

    // Get a maximal flow for the bipartite flow network
    std::vector<NodeID> source_set;
    int max_flow = max_flow_algorithm.solve_max_flow_min_cut(flow_network, source, sink, true, source_set);
    int min_cut =0;

    std::vector<bool> nodes_reached(bipartite.number_of_nodes(), false);
    for(int i = 0; i < source_set.size(); i++ ) {
            if (source_set[i] == source) continue;
            nodes_reached[source_set[i]] = true;
    }
    forall_nodes(bipartite, node) {
            if (bipartite.getPartitionIndex(node) == 0 && !nodes_reached[node]) {
                    min_cut += bipartite.getNodeWeight(node);
                    vertex_cover.push_back(mapping[node]);
            } 
            if (bipartite.getPartitionIndex(node) == 1 && nodes_reached[node]) {
                    min_cut += bipartite.getNodeWeight(node);
                    vertex_cover.push_back(mapping[node]);
            } 

    } endfor
    
    ASSERT_EQ(max_flow, min_cut);

}

void bipartite_weighted_vertex_cover::create_bipartite(graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, graph_access & bipartite, std::vector<NodeID> & bipartite_mapping) {
    std::vector<NodeID> xadj;
    std::vector<NodeID> adjncy;
    std::vector<int> reverse_mapping(G.number_of_nodes(), -1);

    unsigned int lhs_size = 0;
    unsigned int rhs_size = 0;

    // Create the mapping
    forall_nodes(G, node) {
        if (lhs[node] == 0 && rhs[node] == 0) continue;
        else {
            if (lhs[node] == 1) lhs_size++;
            if (rhs[node] == 1) rhs_size++;
            reverse_mapping[node] = bipartite_mapping.size();
            bipartite_mapping.push_back(node);
        }
    } endfor

    // Extract bipartite graph
    unsigned int edges = 0;
    unsigned int nodes = 0;

    for (unsigned int i = 0; i < bipartite_mapping.size(); ++i) {
        NodeID original_node = bipartite_mapping[i];
        xadj.push_back(edges);
        nodes++;
        forall_out_edges(G, edge, original_node) {
            NodeID original_target = G.getEdgeTarget(edge);
            if (lhs[original_node] == 1 && rhs[original_target] == 1) {
                edges++;
                adjncy.push_back(reverse_mapping[original_target]);
            }
            else if (rhs[original_node] == 1 && lhs[original_target] == 1) {
                edges++;
                adjncy.push_back(reverse_mapping[original_target]);
            }
        } endfor
    }
    xadj.push_back(edges);

    int bi_size = nodes;
    int* bi_adjncy = new int[edges];
    int bi_xadj[bi_size + 1];
    std::vector<int> bi_vwgt(bi_size +1,0);
    std::vector<int> bi_adjwgt(edges, 1);
    // Create arrays
    for (unsigned int i = 0; i < adjncy.size(); ++i) {
        bi_adjncy[i] = adjncy[i];
    }
    for (unsigned int i = 0; i < xadj.size(); ++i) {
        bi_xadj[i] = xadj[i];
    }
    for (unsigned int i = 0; i < xadj.size(); ++i) {
        bi_vwgt[i] = G.getNodeWeight(i);
    }

    // Build the graph
    bipartite.build_from_metis_weighted(bi_size, bi_xadj, bi_adjncy, &bi_vwgt[0], &bi_adjwgt[0]);

    // Set the partition indices
    forall_nodes(bipartite, node) {
        NodeID original_node = bipartite_mapping[node];
        if (lhs[original_node] == 1) bipartite.setPartitionIndex(node, 0);
        else if (rhs[original_node] == 1) bipartite.setPartitionIndex(node, 1);
    } endfor

    ASSERT_EQ(bipartite.number_of_nodes(), lhs_size + rhs_size);
    delete[] bi_adjncy;
}

void bipartite_weighted_vertex_cover::create_flow_network(graph_access & bipartite, flow_graph & flow_network, NodeID & source, NodeID & sink) {
    int n = bipartite.number_of_nodes();
    std::vector<bool> included(bipartite.number_of_edges(),false);
    source = n;
    sink = n+1;
    flow_network.start_construction(n+2);
    forall_nodes(bipartite, node) {
        //edge to source and sink with nodeweight as capacity depending on partition index
        if (bipartite.getPartitionIndex(node) == 0) {
            flow_network.new_edge(source, node, bipartite.getNodeWeight(node));
        } else {
            flow_network.new_edge(node, sink, bipartite.getNodeWeight(node));
        }

        //edge between both nodes with inf capacity
        forall_out_edges(bipartite,e,node) {
                NodeID t = bipartite.getEdgeTarget(e);
                flow_network.new_edge(node,t,std::numeric_limits<long int>::max());
            } endfor
    } endfor

    flow_network.finish_construction();

}
