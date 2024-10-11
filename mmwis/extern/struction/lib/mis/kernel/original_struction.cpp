//
// Created by alex on 16.01.20.
//

#include "mwis_finder.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "original_struction.h"

using namespace struction;

typedef branch_and_reduce_algorithm::IS_status IS_status;

template<bool modified>
bool original_struction<modified>::reduce(branch_and_reduce_algorithm* br_alg, NodeID n, size_t max_nodes) {
    this->br_alg = br_alg;
    auto &status = br_alg->status;
    for (NodeID v : status.graph[n]) {
        if (status.weights[v] <= status.weights[n])
            return false;
    }

    auto &additional_nodes = br_alg->buffers[1];
    auto &layer_prefix_count = br_alg->buffers[2];
    auto &marks = br_alg->set_1;
    auto &neighbor_marks = br_alg->set_2;

    if (!find_additional_nodes( additional_nodes, status.graph[n], layer_prefix_count, max_nodes))
        return false;

    br_alg->set(n, IS_status::folded, false);

    const NodeID first_additional_node = status.graph.size();
    br_alg->push_nodes(additional_nodes.size());
    auto &neighbors = status.graph[n];
    neighbor_marks.clear();
    for (size_t i = 0; i < additional_nodes.size(); ++i) {
        status.weights[first_additional_node + i] = modified ? status.weights[additional_nodes[i]] : status.weights[n];
    }
    for (const auto neighbor : neighbors) {
        status.weights[neighbor] -= status.weights[n];
        ///TODO: Fix this?
        if (modified)
            neighbor_marks.add(neighbor);
    }


    for (size_t layer = 0; layer < neighbors.size(); ++layer) {
        for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
            // connect additional nodes with additional nodes of different layers
            for (size_t j = layer_prefix_count[layer + 1]; j < additional_nodes.size(); ++j) {
                status.graph.add_edge_undirected(first_additional_node + i, first_additional_node + j);
            }
            // connect additional nodes with additional nodes of same layer
            for (size_t j = i + 1; j < layer_prefix_count[layer + 1]; ++j) {
                if (status.graph.are_connected(additional_nodes[i], additional_nodes[j])) {
                    status.graph.add_edge_undirected(first_additional_node + i, first_additional_node + j);
                }
            }
        }

        // connect additional nodes with original nodes (adjacent to "i" = neighbors[layer])
        marks.clear();
        for (const auto neighbor : status.graph[neighbors[layer]]) {
            if (neighbor_marks.get(neighbor) || neighbor >= first_additional_node) continue;

            marks.add(neighbor);
            for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
                status.graph.add_edge_undirected(first_additional_node + i, neighbor);
            }
        }

        //connect additional nodes with original nodes (adjacent to "j" = additional_nodes[i])
        for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
            for (const auto neighbor : status.graph[additional_nodes[i]]) {
                if (neighbor_marks.get(neighbor) || marks.get(neighbor) || neighbor >= first_additional_node) continue;

                status.graph.add_edge_undirected(first_additional_node + i, neighbor);
            }
        }
    }

    if (modified) {
        // Form clique
        for (NodeID u : neighbors) {
            marks.clear();
            for (NodeID v : status.graph[u]) {
                marks.add(v);
            }
            for (NodeID w : neighbors) {
                if (w == u || marks.get(w)) continue;

                status.graph.add_edge_undirected(u, w);
            }
        }
        for (size_t layer = 0; layer < neighbors.size(); ++layer) {
            for (size_t layer_ = 0; layer_ < neighbors.size(); ++layer_) {
                if (layer_ == layer) continue;

                for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
                    status.graph.add_edge_undirected(first_additional_node + i, neighbors[layer_]);
                }
            }
        }
    }
    for (NodeID u = first_additional_node; u < status.graph.size(); ++u) {
        br_alg->add_next_level_node(u);
    }
    status.reduction_offset += status.weights[n];
}

template<bool modified>
bool original_struction<modified>::find_additional_nodes(sized_vector<NodeID> &additional_nodes, neighbor_list &neighbors,
                                               sized_vector<NodeID> &layer_prefix_count, size_t max_nodes) const {
    auto &status = br_alg->status;
    additional_nodes.clear();
    layer_prefix_count.set_size(neighbors.size() + 1);
    layer_prefix_count[0] = 0;
    for (size_t i = 0; i < neighbors.size(); ++i) {
        layer_prefix_count[i + 1] = layer_prefix_count[i];
        for (size_t j = i + 1; j < neighbors.size(); ++j) {
            if (!status.graph.are_connected(neighbors[i], neighbors[j])) {
                additional_nodes.push_back(neighbors[j]);
                ++layer_prefix_count[i + 1];
                if (additional_nodes.size() > max_nodes)
                    return false;
            }
        }
    }

    return true;
}

template class original_struction<true>;
template class original_struction<false>;
