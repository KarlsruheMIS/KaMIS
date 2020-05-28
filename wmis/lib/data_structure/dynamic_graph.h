/**
 * dynamic_graph.h
 * Purpose: Dynamic graph datastructure which allows hiding nodes and restoring
 * 			them in reverse order.
 * 
 *****************************************************************************/

#ifndef _DYNAMIC_GRAPH_H
#define _DYNAMIC_GRAPH_H

#include <vector>
#include <algorithm>
#include "graph_access.h"

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
		void resize(size_t size) {neighbors.resize(size);}

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

	void add_edge_directed(NodeID source, NodeID target) {
		graph[source].neighbors.push_back(target);
		graph[source].counter++;
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
		for (size_t pos = 0; pos < slot.counter; pos++) {
			if (slot.neighbors[pos] == target) {
				std::swap(slot.neighbors[pos], slot.neighbors[--slot.counter]);
				return;
			}
		}
	}
	
	// restores last hidden node
	void restore_node(NodeID node) {
		for (auto neighbor : graph[node]) {
			graph[neighbor].counter++;
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

private:
	std::vector<neighbor_list> graph;
};

#endif 
