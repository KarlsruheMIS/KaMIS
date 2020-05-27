/**
 * bipartite_vertex_cover.cpp
 * Purpose: Create a minimal vertex cover for a given graph.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 *****************************************************************************/

#include "bipartite_vertex_cover.h"

#include <limits>
#include <queue>

#include "random_functions.h"

bipartite_vertex_cover::bipartite_vertex_cover() {
    
}

bipartite_vertex_cover::~bipartite_vertex_cover() {

}

bool bipartite_vertex_cover::hopcroft_bfs(graph_access & G, Matching & edge_matching, std::vector<int> & dist) {
    int INF = std::numeric_limits<int>::max();
    NodeID NIL = G.number_of_nodes();
    std::queue<NodeID> bfs_queue;

    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) {
            if (edge_matching[node] == NIL) {
                dist[node] = 0;
                bfs_queue.push(node);
            }
            else dist[node] = INF;
        }
    } endfor
    
    dist[NIL] = INF;

    while (!bfs_queue.empty()) {
        NodeID v = bfs_queue.front();
        bfs_queue.pop();
        if (dist[v] < dist[NIL]) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                // if (dist[edge_matching[u]] == dist[u])
                if (dist[edge_matching[u]] == INF) {
                    dist[edge_matching[u]] = dist[v] + 1;
                    bfs_queue.push(edge_matching[u]);
                }
            } endfor
        }
    }

    return dist[NIL] != INF;
}

bool bipartite_vertex_cover::hopcroft_dfs(graph_access & G, Matching & edge_matching, std::vector<int> & dist, std::vector<bool> & T, NodeID v) {
    int INF = std::numeric_limits<int>::max();
    NodeID NIL = G.number_of_nodes();

    if (v != NIL) {
        forall_out_edges(G, e, v) {
            NodeID u = G.getEdgeTarget(e);
            if (dist[edge_matching[u]] == dist[v] + 1) {
                if (hopcroft_dfs(G, edge_matching, dist, T, edge_matching[u])) {
                    edge_matching[u] = v;
                    edge_matching[v] = u;
                    return true;
                }
            }
        } endfor
        dist[v] = INF;
        return false;
    }
    return true;
}

void bipartite_vertex_cover::hopcroft_koenig(graph_access & G, Matching & edge_matching, std::vector<bool> & T) {
    NodeID NIL = G.number_of_nodes();

    std::vector<NodeID> L;
    std::vector<NodeID> R;
    forall_nodes(G, v) {
        if (G.getPartitionIndex(v) == 0) {
            if (edge_matching[v] == NIL) {
                T[v] = true;
                L.push_back(v);
            }
        }
    } endfor

    unsigned int nodes_added = 0;

    do {
        nodes_added = 0;
        for (unsigned int i = 0; i < L.size(); ++i) {
            NodeID v = L[i];
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                if (!T[u] && edge_matching[v] != u) {
                    T[u] = true;
                    R.push_back(u);
                    nodes_added++;
                }
            } endfor
        }
        L.clear();
        for (unsigned int i = 0; i < R.size(); ++i) {
            NodeID v = R[i];
            NodeID u = edge_matching[v];
            if (!T[u] && u != NIL) {
                T[u] = true;
                L.push_back(u);
                nodes_added++;
            }
        }
        R.clear();
    } while (nodes_added > 0);

}

int bipartite_vertex_cover::hopcroft(graph_access & G, Matching & edge_matching, std::vector<bool> & T) {
    int INF = std::numeric_limits<int>::max();
    NodeID NIL = G.number_of_nodes();
    std::vector<int> dist(G.number_of_nodes() + 1, INF);

    forall_nodes(G, v) {
        edge_matching[v] = NIL;
    } endfor
    

    NodePermutationMap permutation;
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);

    unsigned int matching = 0;
    while (hopcroft_bfs(G, edge_matching, dist)) {
        forall_nodes(G, node) {
            NodeID v = permutation[node];
            if (G.getPartitionIndex(v) == 0) {
                if (edge_matching[v] == NIL) {
                    if (hopcroft_dfs(G, edge_matching, dist, T, v)) {
                        matching++;
                    }
                }
            }
        } endfor
    }

    hopcroft_koenig(G, edge_matching, T);

    return matching;
}

void bipartite_vertex_cover::hopcroft_cover(graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, std::vector<NodeID> & vertex_cover) {
    // Create the bipartite graph
    graph_access bipartite;
    std::vector<NodeID> mapping;
    create_bipartite(G, lhs, rhs, bipartite, mapping);

    // Get a maximal matching for the bipartite graph
    Matching edge_matching(bipartite.number_of_nodes());
    std::vector<bool> T(bipartite.number_of_nodes(), false);
    // unsigned int matching_size = hopcroft(bipartite, edge_matching, T);

    forall_nodes(bipartite, node) {
        if (bipartite.getPartitionIndex(node) == 0 && !T[node]) {
            vertex_cover.push_back(mapping[node]);
            bipartite.setPartitionIndex(node, 1);
        }
        else if (bipartite.getPartitionIndex(node) == 1 && T[node]) {
            vertex_cover.push_back(mapping[node]);
            bipartite.setPartitionIndex(node, 1);
        }
        else {
            bipartite.setPartitionIndex(node, 0);
        }
    } endfor

    // ASSERT_EQ(matching_size, vertex_cover.size());
}

void bipartite_vertex_cover::create_bipartite(graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, graph_access & bipartite, std::vector<NodeID> & bipartite_mapping) {
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
    int bi_adjncy[edges];
    int bi_xadj[bi_size + 1];
    // Create arrays
    for (unsigned int i = 0; i < adjncy.size(); ++i) {
        bi_adjncy[i] = adjncy[i];
    }
    for (unsigned int i = 0; i < xadj.size(); ++i) {
        bi_xadj[i] = xadj[i];
    }

    // Build the graph
    bipartite.build_from_metis(bi_size, bi_xadj, bi_adjncy);

    // Set the partition indices
    forall_nodes(bipartite, node) {
        NodeID original_node = bipartite_mapping[node];
        if (lhs[original_node] == 1) bipartite.setPartitionIndex(node, 0);
        else if (rhs[original_node] == 1) bipartite.setPartitionIndex(node, 1);
    } endfor

    ASSERT_EQ(bipartite.number_of_nodes(), lhs_size + rhs_size);
}

