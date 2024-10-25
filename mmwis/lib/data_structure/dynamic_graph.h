/**
 * dynamic_graph.h
 * Purpose: Dynamic graph datastructure which allows hiding nodes and restoring
 * 			them in reverse order.
 * 
 *****************************************************************************/

#pragma once
#include <iostream>
#include <cstddef>
#include <type_traits>
#include <vector>
#include <algorithm>
#include "graph_access.h"

namespace mmwis {
class dynamic_graph {
public:
	class neighbor_list {
		friend dynamic_graph;

	public:
		using iterator = std::vector<NodeID>::iterator;
		using const_iterator = std::vector<NodeID>::const_iterator;

		neighbor_list() = default;
		neighbor_list(size_t size) : neighbors(size) {};
		neighbor_list(const std::vector<NodeID>& neighbors) : neighbors(neighbors), counter(neighbors.size()) {}
		neighbor_list(std::vector<NodeID>&& neighbors) : neighbors(std::move(neighbors)), counter(this->neighbors.size()) {}

		iterator begin() { return neighbors.begin(); }
		iterator end() { return neighbors.begin() + counter; }
		const_iterator begin() const { return neighbors.begin(); }
		const_iterator end() const { return neighbors.begin() + counter; }
		const_iterator cbegin() const { return neighbors.cbegin(); }
		const_iterator cend() const { return neighbors.cbegin() + counter; }

		size_t size() const noexcept { return counter; }
		size_t capacity() const noexcept { return neighbors.capacity(); }
		void resize(size_t size) {neighbors.resize(size);}
		void reserve(size_t size) {neighbors.reserve(size);}

		NodeID& operator[] (size_t index) { return neighbors[index]; }
		const NodeID& operator[] (size_t index) const { return neighbors[index]; }

	private:
		std::vector<NodeID> neighbors;
		size_t counter = 0;
	};

	dynamic_graph(size_t nodes = 0) : graph(nodes) { graph.reserve(nodes); }

	dynamic_graph(graph_access& G) : graph(G.number_of_nodes()) {
		neighbor_list* slot;

		forall_nodes(G, node) {
			slot = &graph[node];
			slot->resize(G.getNodeDegree(node));

			forall_out_edges(G, edge, node) {
				slot->neighbors[slot->counter++] = G.getEdgeTarget(edge);
			} endfor
		} endfor
	}

    //build only subgraph with given partitionIndex
	dynamic_graph(graph_access& G, sized_vector<NodeID> &nodes, std::vector<NodeID> &reverse_mapping, std::vector<NodeID> &partition_map) : graph(nodes.size()) {
	    reverse_mapping.resize(nodes.size(), 0);
		neighbor_list* slot;
        fast_set neighbors(G.number_of_nodes());

        NodeID partition_node = 0;
        for(NodeID node : nodes) {
            assert(partition_map[node] == G.number_of_nodes() && "Partition map can not be overwritten if set once.");
            partition_map[node] = partition_node;
            reverse_mapping[partition_node] = node;
            partition_node++;
        }


		for(NodeID node : nodes) {
            neighbors.add(node);
        }
        
        for(NodeID node : nodes) {
			slot = &graph[partition_map[node]];
            int degree = G.getNodeDegree(node);
            slot->resize(degree);

            int partition_degree = 0;
			forall_out_edges(G, edge, node) {
                NodeID target = G.getEdgeTarget(edge);
                NodeID partition_target = partition_map[G.getEdgeTarget(edge)];
                // if target was not mapped befor, it is not in the partition
                /* if (partition_target == G.number_of_nodes()) continue; */
                if (neighbors.get(target)) {
                    partition_degree++;
				    slot->neighbors[slot->counter++] = partition_target;
                }
			} endfor
            // in case of boundary nodes:
            if (partition_degree < degree) slot->resize(partition_degree);
            /* print_neighbors(partition_map[node]); */
		} 
	}

	dynamic_graph(const std::vector<std::vector<NodeID>>& adj) {
		graph.reserve(adj.size());

		for (const auto& vec : adj) {
			graph.push_back(vec);
		}
	}

	dynamic_graph(std::vector<std::vector<NodeID>>&& adj) {
		graph.reserve(adj.size());

		for (auto& vec : adj) {
			graph.push_back(std::move(vec));
		}
	}

	void add_edge_directed(NodeID source, NodeID target, bool not_hidden = false) {
		auto& slot = graph[source];
		slot.neighbors.push_back(target);
		slot.counter++;

        if (not_hidden) {
            for (size_t hidden_pos = slot.neighbors.size()-1; hidden_pos >= slot.counter; hidden_pos--) {
		        std::swap(slot.neighbors[hidden_pos-1], slot.neighbors[hidden_pos]);
            }
        }
	}

	void add_edge_undirected(NodeID node1, NodeID node2, bool not_hidden = false) {
		add_edge_directed(node1, node2, not_hidden);
		add_edge_directed(node2, node1, not_hidden);
	}

    void remove_edge_undirected(NodeID node1, NodeID node2) {
        remove_edge_directed(node1, node2);
        remove_edge_directed(node2, node1);
    }

    void remove_edge_directed (NodeID source, NodeID target) {
        hide_edge_to_last(source, target);
        if (graph[source].neighbors.back() == target) { graph[source].neighbors.pop_back(); 
        } 
        else {
        /* std::cout << "edge "<< source << "--" << target <<" was not removed." << std::endl; */
        /* std::cout << "last hidden node: "<< graph[source].neighbors.back() << "!=" << target << std::endl; */
        }

    }

	// node itself still exists in the graph just the edges towards this node are hidden
	void hide_node(NodeID node) {
		for (auto neighbor : graph[node]) {
			hide_edge(neighbor, node);
		}
	}

	void hide_edge_undirected(NodeID source, NodeID target) {
		hide_edge(source, target);
		hide_edge(target, source);
	}

	void hide_edge(NodeID source, NodeID target) {
		auto& slot = graph[source];
		for (size_t pos = 0; pos < slot.counter; pos++) {
			if (slot.neighbors[pos] == target) {
				std::swap(slot.neighbors[pos], slot.neighbors[--slot.counter]);
				return;
			}
		}
	}

	void hide_edge_to_last(NodeID source, NodeID target) {
		auto& slot = graph[source];
		for (size_t pos = 0; pos < slot.counter; pos++) {
			if (slot.neighbors[pos] == target) {
				std::swap(slot.neighbors[pos], slot.neighbors[--slot.counter]);
                for (size_t hidden_pos = slot.counter+1; hidden_pos < slot.neighbors.size(); hidden_pos++) {
				    std::swap(slot.neighbors[hidden_pos-1], slot.neighbors[hidden_pos]);
                }
				return;
			}
		}
 		for (size_t pos = slot.counter; pos < slot.neighbors.size(); pos++) {
			if (slot.neighbors[pos] == target) {
                for (size_t hidden_pos = pos+1; hidden_pos < slot.neighbors.size(); hidden_pos++) {
                    std::swap(slot.neighbors[hidden_pos-1], slot.neighbors[hidden_pos]);
                }
        std::cerr << "ERROR? edge "<< source << ", " << target <<" already hidden" << std::endl;
        return;
			}
		}
        std::cerr << "ERROR edge "<< source << ", " << target <<" to hide not found" << std::endl;
    }
	
	// restores last hidden node
	void restore_node(NodeID node) {
		// std::cout << "restore node "<< node << std::endl;
		for (auto neighbor : graph[node]) {
			// print_all_neighbors(neighbor);
			graph[neighbor].counter++;
			// print_all_neighbors(neighbor);
		}
	}

	// restores last hidden edge starting from node
	void restore_edge(NodeID node) {
		graph[node].counter++;
	}

	// replaces the last target of the restored edge of node
	void replace_last_restored_edge(NodeID node, NodeID replacement) {
		*(graph[node].end() - 1) = replacement;
	}

	// restores last hidden edge starting from node and replaces the target
	void restore_edge_and_replace(NodeID node, NodeID replacement) {
		graph[node].neighbors[graph[node].counter++] = replacement;
	}

	size_t size() const noexcept { return graph.size(); }

	neighbor_list& operator[] (NodeID node) { return graph[node]; }

	const neighbor_list& operator[] (NodeID node) const { return graph[node]; }

    void print_all_neighbors(NodeID node) {
		std::cout << "N("<<node<<"): ";
		for (size_t i = 0; i < graph[node].counter; i++) {
			std::cout << graph[node].neighbors[i] << " ";
		}
		std::cout << " | ";
		for (size_t i = graph[node].counter; i < graph[node].size(); i++) {
			std::cout << graph[node].neighbors[i] << " ";
		}
			std::cout << "\n";
	}

    void print_neighbors(NodeID node) {
        std::cout << "N("<<node<<"): ";
        for (size_t i = 0; i < graph[node].counter; i++) {
            std::cout << graph[node].neighbors[i] << " ";
        }
            std::cout << "\n";
    } 

private:
	std::vector<neighbor_list> graph;
};
}

