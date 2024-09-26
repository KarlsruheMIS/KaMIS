//
// Created by alex on 14.12.19.
//

#include <stack>
#include <algorithm>
#include "graph_access.h"
#include "fast_set.h"
#include "maxNodeHeap.h"
#include "mwis_finder.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "extended_struction.h"

using namespace struction;

typedef branch_and_reduce_algorithm::IS_status IS_status;

template<bool reduced>
bool extended_struction<reduced>::reduce(branch_and_reduce_algorithm* br_alg, NodeID n, size_t max_nodes) {
    this->br_alg = br_alg;
    auto &status = br_alg->status;

    auto &additional_nodes = br_alg->buffers[1];
    auto &layer_prefix_count = br_alg->buffers[2];
    auto &tmp_marks = br_alg->set_1;
    if (!finder.findAllMWIS<reduced>(br_alg, status.graph[n], status.weights[n], max_nodes)) return false;

    std::vector<mwis> &sets = finder.get_sets();

    if (reduced && !findAdditionalNodes(n, sets, br_alg->buffers[3], layer_prefix_count, additional_nodes, max_nodes - sets.size())) {
        return false;
    }

    // Create set nodes
    const NodeID first_set_node = status.graph.size();
    br_alg->push_nodes(sets.size());
    for (size_t i = 0; i < sets.size(); ++i) {
        status.weights[first_set_node + i] = sets[i].weight - status.weights[n];
    }

    const NodeID first_additional_node = status.graph.size();
    //We have to obtain adjacencies between additional nodes before we remove N[n]
    if (reduced) {
        // Create additional nodes
        br_alg->push_nodes(additional_nodes.size());
        for (size_t i = 0; i < additional_nodes.size(); ++i) {
            status.weights[first_additional_node + i] = status.weights[additional_nodes[i]];
        }

        // connect additional nodes with additional nodes
        for (size_t layer = 0; layer < sets.size(); ++layer) {
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
        }
    }
    // Delete v and its neighbors
    br_alg->set(n, IS_status::folded, true);
    for (const auto u : status.graph[n])
        br_alg->set(u, IS_status::folded, false);


    // Connect set nodes to original nodes
    for (size_t i = 0; i < sets.size(); ++i) {
        tmp_marks.clear();
        for (const auto n : sets[i].nodes) {
            for (const auto neighbor : status.graph[n]) {
                if (status.node_status[neighbor] != IS_status::not_set || !tmp_marks.add(neighbor)) continue;
                status.graph.add_edge_undirected(first_set_node + i, neighbor);
            }
        }
    }

    if (reduced) {
        // connect additional nodes with original nodes
        for (size_t layer = 0; layer < sets.size(); ++layer) {
            tmp_marks.clear();
            for (const auto neighbor : status.graph[first_set_node + layer]) {
                for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
                    status.graph.add_edge_undirected(first_additional_node + i, neighbor);
                }
                tmp_marks.add(neighbor);
            }
            for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
                for (const auto neighbor : status.graph[additional_nodes[i]]) {
                    if (neighbor >= first_set_node || status.node_status[neighbor] != IS_status::not_set || tmp_marks.get(neighbor)) continue;
                    status.graph.add_edge_undirected(first_additional_node + i, neighbor);
                }
            }
        }

        // connect additional nodes with set nodes
        for (size_t layer = 0; layer < sets.size(); ++layer) {
            for (size_t layer_ = 0; layer_ < sets.size(); ++layer_) {
                if (layer == layer_) continue;
                for (size_t i = layer_prefix_count[layer]; i < layer_prefix_count[layer + 1]; ++i) {
                    status.graph.add_edge_undirected(first_additional_node + i, first_set_node + layer_);
                }
            }
        }
    }

    // form clique by all set nodes
    for (size_t i = 0; i < sets.size(); ++i)
        for (size_t j = i + 1; j < sets.size(); ++j)
            status.graph.add_edge_undirected(first_set_node + i, first_set_node + j);

    // Add all new nodes to modified list
    for (NodeID u = first_set_node; u < status.graph.size(); ++u) {
        br_alg->add_next_level_node(u);
    }

    status.reduction_offset += status.weights[n];
    restore_vec.emplace_back(n, sets.size());

    return true;
}

template<bool reduced_version>
void extended_struction<reduced_version>::restore(branch_and_reduce_algorithm* br_alg) {
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    // remove new nodes
    br_alg->pop_nodes(data.set_nodes);

    // restore neighbors of main node in reversed order
    auto &neighbors = status.graph[data.main];
    for (size_t i = neighbors.size(); i-- >= 1; ) {
        NodeID neighbor = neighbors[i];
        status.node_status[neighbor] = IS_status::not_set;
        status.graph.restore_node(neighbor);
    }
    status.remaining_nodes += status.graph[data.main].size();
    br_alg->unset(data.main);
    status.reduction_offset -= status.weights[data.main];
    restore_vec.pop_back();
}

template<bool reduced_version>
void extended_struction<reduced_version>::apply(branch_and_reduce_algorithm* br_alg) {
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    NodeID main = data.main;
    auto &neighbors = status.graph[main];
    for (NodeID n = status.n - data.set_nodes; n < status.n; ++n) {
        if (status.node_status[n] == IS_status::included) {
            NodeWeight set_weight = status.weights[n] + status.weights[main];
            restore(br_alg);
            includable_neighbors.clear();

            for (const auto neighbor : neighbors) {
                bool includable = true;
                for (NodeID v : status.graph[neighbor]) {
                    if (status.node_status[v] == IS_status::included) {
                        includable = false;
                        break;
                    }
                }
                if (includable)
                    includable_neighbors.push_back(neighbor);
            }
            finder.findAllMWIS(br_alg, includable_neighbors, set_weight - 1, 1);

            status.node_status[main] = IS_status::excluded;
            for (NodeID n : neighbors)
                status.node_status[n] = IS_status::excluded;
            for (NodeID n : finder.get_sets().back().nodes) {
                status.is_weight += status.weights[n];
                status.node_status[n] = IS_status::included;
            }


            status.is_weight -= status.weights[n];
            return;
        }
    }

    //No set node included --> include main node.
    restore(br_alg);
    status.node_status[main] = IS_status::included;
    status.is_weight += status.weights[main];
    for (NodeID n : neighbors)
        status.node_status[n] = IS_status::excluded;
}

template<bool reduced_version>
bool extended_struction<reduced_version>::findAdditionalNodes(NodeID n, const std::vector<mwis> &independent_sets, sized_vector<NodeID> &neighbor_ids,
                                                              sized_vector<NodeID> &layer_prefix_count, sized_vector<NodeID> &additional_nodes, size_t max_nodes) {
    auto &status = br_alg->status;
    auto &marks = br_alg->set_1;
    auto &neighbors = status.graph[n];

    additional_nodes.clear();
    layer_prefix_count.set_size(independent_sets.size() + 1);
    layer_prefix_count[0] = 0;

    neighbor_ids.set_size(status.n);
    for (size_t i = 0; i < neighbors.size(); ++i) {
        neighbor_ids[neighbors[i]] = i;
    }
    for (size_t i = 0; i < independent_sets.size(); ++i) {
        layer_prefix_count[i + 1] = layer_prefix_count[i];
        const auto &mwis_nodes = independent_sets[i].nodes;

        marks.clear();
        for (const NodeID n : mwis_nodes) {
            for (const NodeID v : status.graph[n])
                marks.add(v);
        }
        for (size_t j = neighbor_ids[mwis_nodes.back()] + 1; j < neighbors.size(); ++j) {
            if (marks.get(neighbors[j])) continue;

            ++layer_prefix_count[i + 1];
            additional_nodes.push_back(neighbors[j]);
            if (additional_nodes.size() > max_nodes)
                return false;
        }
    }

    return true;
}

template<bool reduced_version>
size_t extended_struction<reduced_version>::removed_vertices(branch_and_reduce_algorithm* br_alg, NodeID n) {
    return br_alg->deg(n) + 1;
}

template class extended_struction<true>;
template class extended_struction<false>;
