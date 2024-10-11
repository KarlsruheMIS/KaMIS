/**
 * struction_dynamic_graph.h
 * Purpose: Dynamic graph datastructure which allows hiding nodes and restoring
 * 			them in reverse order.
 * 
 *
 ******************************************************************************
 * Copyright (C) 2015-2018 Robert Williger
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

#ifndef _struction_DYNAMIC_GRAPH_H
#define _struction_DYNAMIC_GRAPH_H

#include <vector>
#include <algorithm>
#include "graph_access.h"
#include "fast_set.h"
#include "sized_vector.h"
namespace struction {
class struction_dynamic_graph {
public:
    using neighbor_list = std::vector<NodeID>;
    /*class neighbor_list {
        friend dynamic_graph;

    public:
        using iterator = std::vector<NodeID>::iterator;
        using const_iterator = std::vector<NodeID>::const_iterator;

        neighbor_list() = default;
        neighbor_list(size_t size) : neighbors(size) {};
        neighbor_list(const std::vector<NodeID>& neighbors) : neighbors(neighbors), counter(neighbors.size()) {}
        neighbor_list(std::vector<NodeID>&& neighbors) : neighbors(std::move(neighbors)), counter(this->neighbors.size()) {}

        NodeID &back() {return neighbors[counter - 1];}
        iterator begin() { return neighbors.begin(); }
        iterator end() { return neighbors.begin() + counter; }
        const_iterator begin() const { return neighbors.begin(); }
        const_iterator end() const { return neighbors.begin() + counter; }
        const_iterator cbegin() const { return neighbors.cbegin(); }
        const_iterator cend() const { return neighbors.cbegin() + counter; }

        size_t size() const noexcept { return counter; }
        void resize(size_t size) {neighbors.resize(size);}

        NodeID& operator[] (size_t index) { return neighbors[index]; }
        const NodeID& operator[] (size_t index) const { return neighbors[index]; }

    private:
        std::vector<NodeID> neighbors;
        size_t counter = 0;
    };*/

    struction_dynamic_graph(size_t nodes = 0) : graph(nodes), set(nodes) { graph.reserve(nodes); }

    struction_dynamic_graph(graph_access& G) : struction_dynamic_graph(G.number_of_nodes()) {
        neighbor_list* slot;

        forall_nodes(G, node)
            slot = &graph[node];
            slot->reserve(G.getNodeDegree(node));
            forall_out_edges(G, edge, node)
                slot->push_back(G.getEdgeTarget(edge));
            endfor
        endfor
    }

    struction_dynamic_graph(const std::vector<std::vector<NodeID>>& adj) : set(adj.size()) {
        graph.reserve(adj.size());

        for (const auto& vec : adj) {
            graph.push_back(vec);
        }
    }

    struction_dynamic_graph(std::vector<std::vector<NodeID>>&& adj) : set(adj.size()) {
        graph.reserve(adj.size());

        for (auto& vec : adj) {
            graph.push_back(std::move(vec));
        }
    }

    bool are_connected(NodeID u, NodeID v) const {
        for (const NodeID w : graph[u]) {
            if (w == v)
                return true;
        }
        return false;
    }

    void relink_directed(NodeID source, NodeID old_target, NodeID new_target) {
        auto& slot = graph[source];
        for (size_t pos = 0; pos < slot.size(); pos++) {
            if (slot[pos] == old_target) {
                slot[pos] = new_target;
                return;
            }
        }
    }

    void add_edge_directed(NodeID source, NodeID target) {
        graph[source].push_back(target);
    }

    void add_edge_undirected(NodeID node1, NodeID node2) {
        add_edge_directed(node1, node2);
        add_edge_directed(node2, node1);
    }

    // node itself still exists in the graph just the edges towards this node are hidden
    void hide_node(NodeID node) {
        for (auto neighbor : graph[node]) {
            hide_edge(neighbor, node);
        }
    }

    void hide_edge(NodeID source, NodeID target) {
        auto& slot = graph[source];
        for (size_t pos = slot.size(); pos-- > 0; ) {
            if (slot[pos] == target) {
                slot[pos] = slot.back();
                slot.pop_back();
                return;
            }
        }
    }

    template<typename P>
    void hide_edges(NodeID source, P pred, size_t max_edges = std::numeric_limits<size_t>::max()) {
        auto& slot = graph[source];
        for (size_t pos = slot.size(); pos-- > 0; ) {
            if (pred(slot[pos])) {
                slot[pos] = slot.back();
                slot.pop_back();
                if (--max_edges == 0)
                    return;
            }
        }
    }

    // restores last hidden node
    void restore_node(NodeID node) {
        for (auto neighbor : graph[node]) {
            graph[neighbor].push_back(node);
        }
    }

    // replaces the last target of the restored edge of node
    void replace_last_restored_edge(NodeID node, NodeID replacement) {
        graph[node].back() = replacement;
    }

    void push_nodes(size_t nodes) {
        graph.resize(size() + nodes);
        set.resize(size());
    }

    void pop_nodes(size_t nodes) {
        size_t new_size = size() - nodes;
        set.clear();
        for (NodeID n = size() - nodes; n < size(); ++n) {
            for (NodeID v : graph[n]) {
                if (v < new_size && set.add(v)) {
                    hide_edges(v, [new_size](NodeID w) {
                        return w >= new_size;
                    });
                }
            }
        }
        graph.resize(size() - nodes);
    }

    size_t size() const noexcept { return graph.size(); }

    neighbor_list& operator[] (NodeID node) { return graph[node]; }

    const neighbor_list& operator[] (NodeID node) const { return graph[node]; }

private:
    std::vector<neighbor_list> graph;
    fast_set set;
};
}

#endif 
