/******************************************************************************
* reductions.cpp
*****************************************************************************/


#include "reductions.h"
#include "branch_and_reduce_algorithm.h"
#include "data_structure/flow_graph.h"
#include "algorithms/push_relabel.h"

#include <utility>


typedef branch_and_reduce_algorithm::IS_status IS_status;

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			NodeWeight neighbor_weights = 0;

			for (NodeID u : status.graph[v]) {
				neighbor_weights += status.weights[u];
			}

			if (status.weights[v] >= neighbor_weights) {
				br_alg->set(v, IS_status::included);
			}
		}
	}

	//std::cout << "neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighborhood = br_alg->set_1;
	size_t oldn = status.remaining_nodes;

	NodeWeight neighbor_weights;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			neighbor_weights = 0;
			neighbors.clear();
			neighborhood.clear();

			for (NodeID neighbor : status.graph[v]) {
				neighbor_weights += status.weights[neighbor];
				neighbors.push_back(neighbor);
				neighborhood.add(neighbor);
			}

			if (status.weights[v] >= neighbor_weights) {
				br_alg->set(v, IS_status::included);
				continue;
			}

			std::sort(neighbors.begin(), neighbors.end(), [&status](const NodeID lhs, const NodeID rhs) { return status.weights[lhs] > status.weights[rhs]; });

			bool is_reducible = false;

			for (size_t i = 0; i < neighbors.size() && !is_reducible; i++) {
				NodeID neighbor1 = neighbors[i];

				if (!neighborhood.get(neighbor1))
					continue;

				for (NodeID neighbor2 : status.graph[neighbor1]) {
					if (neighbor2 != neighbor1 && neighborhood.get(neighbor2)) {
						// triangle [v, neighbor1, neighbor2] found
						neighbor_weights -= std::min(status.weights[neighbor1], status.weights[neighbor2]);
						neighborhood.remove(neighbor1);
						neighborhood.remove(neighbor2);

						if (status.weights[v] >= neighbor_weights) {
							is_reducible = true;
						}

						break;
					}
				}
			}

			if (is_reducible) {
				br_alg->set(v, IS_status::included);
			}
		}
	}

	//std::cout << "clique neighbor fast redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	this->br_alg = br_alg;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set && partition_into_cliques(v)) {
			br_alg->set(v, IS_status::included);
		}
	}

	//std::cout << "clique neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::partition_into_cliques(NodeID v) {
	auto& status = br_alg->status;
	auto& neighbors_vec = br_alg->buffers[0];
	auto& clique_neighbors_set = br_alg->set_1;

	target_weight = status.weights[v];
	neighbor_weights = 0;
	neighbors_vec.clear();

	for (NodeID neighbor : status.graph[v]) {
		neighbor_weights += status.weights[neighbor];
		neighbors_vec.push_back(neighbor);
	}

	if (neighbor_weights <= target_weight) {
		return true;
	}

	// partition neigbors of v into cliques
	NodeID max_neighbor;
	NodeWeight max_neighbor_weight;

	while (neighbors_vec.size() >= 2) {
		max_neighbor_weight = 0;
		clique_neighbors_set.clear();

		for (auto neighbor : neighbors_vec) {
			if (status.weights[neighbor] > max_neighbor_weight) {
				max_neighbor = neighbor;
				max_neighbor_weight = status.weights[neighbor];
			}
			clique_neighbors_set.add(neighbor);
		}

		clique_neighbors_set.remove(max_neighbor);
		if (expand_clique(max_neighbor, neighbors_vec, clique_neighbors_set))
			return true;
	}

	return false;
}

bool clique_neighborhood_reduction::expand_clique(NodeID max_neighbor, sized_vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set) {
	auto& status = br_alg->status;
	auto& temp_set = br_alg->set_2;
	auto& clique_set = br_alg->double_set;
	clique_set.clear();

	size_t local_max;
	NodeWeight local_max_weight;
	bool intersection_empty;

	while (true) {
		// intersect neighbors of clique with neighbors of max_neighbor
		intersection_empty = true;
		local_max_weight = 0;
		clique_set.add(max_neighbor);
		temp_set.clear();

		for (auto neighbor : status.graph[max_neighbor]) {
			if (clique_neighbors_set.get(neighbor)) {
				temp_set.add(neighbor);
				intersection_empty = false;

				if (status.weights[neighbor] > local_max_weight) {
					local_max = neighbor;
					local_max_weight = status.weights[neighbor];
				}
			}
		}

		if (intersection_empty)
			break;

		// add local_max to current clique
		neighbor_weights -= local_max_weight;
		if (neighbor_weights <= target_weight)
			return true;

		std::swap(clique_neighbors_set, temp_set);
		clique_neighbors_set.remove(local_max);
		max_neighbor = local_max;
	}

	auto& reamining_neighbors = br_alg->buffers[1];
	reamining_neighbors.clear();

	// adjust neigbors_vec
	for (auto neighbor : neighbors_vec) {
		if (!clique_set.get(neighbor))
			reamining_neighbors.push_back(neighbor);
	}

	std::swap(reamining_neighbors, neighbors_vec);
	return false;
}

bool critical_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t n = br_alg->status.n;
	size_t oldn = status.remaining_nodes;

	// build bipartite flow graph
	// node '+ n' shows that we refer to the node in the rest[1] partition
	flow_graph fg;
	fg.start_construction(2 * n + 2);

	const NodeID source = 2 * n;
	const NodeID sink = source + 1;

	for (NodeID node = 0; node < n; node++) {
		if (status.node_status[node] == IS_status::not_set) {
			// add source and target edges
			fg.new_edge(source, node, status.weights[node]);
			fg.new_edge(node + n, sink, status.weights[node]);

			// add edges between node and its neighbors in the other partition
			for (NodeID neighbor : status.graph[node]) {
				// each outgoing edge has enough capacity to support full flow of the single incoming edge into 'node'
				fg.new_edge(node, neighbor + n, status.weights[node]);
			}
		}
	}

	fg.finish_construction();


	// solve max-flow problem
	push_relabel flow_solver;
	std::vector<NodeID> dummy_vec;
	flow_solver.solve_max_flow_min_cut(fg, source, sink, false, dummy_vec);

	auto& max_cs_set = br_alg->double_set;

	max_cs_set.clear();
	// (source, node) edges where flow < capacity indicate that node is in the maximum critical set
	forall_out_edges(fg, edge, source) {
		NodeID node = fg.getEdgeTarget(source, edge);
		if (fg.getEdgeFlow(source, edge) < fg.getEdgeCapacity(source, edge)) {
			max_cs_set.add(node);
		}
	} endfor

		// isolated nodes in the maximum critical set form the maximum independent critical set
		for (NodeID node = 0; node < n; node++) {
			if (status.node_status[node] == IS_status::not_set && max_cs_set.get(node)) {
				bool isolated = true;

				for (NodeID neighbor : status.graph[node]) {
					if (max_cs_set.get(neighbor)) {
						isolated = false;
						break;
					}
				}

				//TODO: check if continue works! (was break)
				if (!isolated) continue;

				// found isolated node
				br_alg->set(node, IS_status::included);
			}
		}

	//std::cout << "cs redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool fold2_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	std::array<NodeID, 2> neighbors;
	size_t neighbor_count;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set && br_alg->deg(v) == 2) {
			neighbor_count = 0;
			bool skip = false;

			for (NodeID u : status.graph[v]) {
				if (status.weights[u] > status.weights[v]) {
					skip = true;
					break;
				}

				neighbors[neighbor_count++] = u;
			}

			if (skip) continue;

			if (status.weights[v] >= status.weights[neighbors[0]] + status.weights[neighbors[1]]) {
				// v is always best choice of the three vertices
				br_alg->set(v, IS_status::included);
				continue;
			}

			bool is_clique = false;
			for (NodeID u : status.graph[neighbors[0]]) if (u == neighbors[1]) {
				// clique of size 3
				if (status.weights[v] >= status.weights[u] && status.weights[v] >= status.weights[neighbors[0]]) {
					br_alg->set(v, IS_status::included);
				}
				is_clique = true;
				break;
			}

			if (is_clique) continue;

			fold(br_alg, { v,{ neighbors[0], neighbors[1] } });
		}
	}

	//std::cout << "fold2 redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void fold2_reduction::fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;

	br_alg->set(nodes.rest[1], IS_status::folded, false);
	br_alg->set(nodes.rest[0], IS_status::folded, true);

	restore_vec.push_back({ nodes, status.weights[nodes.main], status.graph[nodes.main],{} });

	status.reduction_offset += status.weights[nodes.main];
	status.weights[nodes.main] = status.weights[nodes.rest[0]] + status.weights[nodes.rest[1]] - status.weights[nodes.main];

	std::vector<NodeID> new_neighbors;
	neighbors.clear();
	neighbors.add(nodes.main);

	for (size_t i = 0; i < 2; i++) {
		for (auto neighbor : status.graph[nodes.rest[i]]) {
			if (neighbors.add(neighbor)) {
				new_neighbors.push_back(neighbor);
				status.graph.restore_edge_and_replace(neighbor, nodes.main);
				restore_vec.back().node_vecs[i].push_back(neighbor);
			}
		}
	}

	status.graph[nodes.main] = dynamic_graph::neighbor_list(std::move(new_neighbors));
	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.main);
	br_alg->add_next_level_neighborhood(nodes.main);
}

void fold2_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	// is "restored" in following loop
	status.graph.hide_node(data.nodes.main);
	status.graph[data.nodes.main] = std::move(data.main_neighbor_list);

	for (size_t i = 0; i < 2; i++) {
		br_alg->unset(data.nodes.rest[i]);

		for (auto neighbor : data.node_vecs[i]) {
			status.graph.replace_last_restored_edge(neighbor, data.nodes.rest[i]);
		}
	}

	status.weights[data.nodes.main] = data.main_weight;
	status.reduction_offset -= data.main_weight;

	restore_vec.pop_back();
}

void fold2_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
	auto main_status = status.node_status[nodes.main];
	restore(br_alg);

	if (main_status == IS_status::included) {
		status.node_status[nodes.main] = IS_status::excluded;
		status.node_status[nodes.rest[0]] = IS_status::included;
		status.node_status[nodes.rest[1]] = IS_status::included;

		status.is_weight += status.weights[nodes.rest[0]] + status.weights[nodes.rest[1]];
	} else {
		status.node_status[nodes.main] = IS_status::included;
		status.node_status[nodes.rest[0]] = IS_status::excluded;
		status.node_status[nodes.rest[1]] = IS_status::excluded;

		status.is_weight += status.weights[nodes.main];
	}
}

bool clique_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& set_1 = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	auto& isolated = br_alg->buffers[1];
	std::vector<NodeID> non_isolated;

	size_t oldn = status.remaining_nodes;

	for (size_t node_idx = 0; node_idx < marker.current_size(); node_idx++) {
		NodeID node = marker.current_vertex(node_idx);

		if (status.node_status[node] == IS_status::not_set) {
			neighbors.clear();
			set_1.clear();
			set_1.add(node);

			// find potential clique
			for (NodeID neighbor : status.graph[node]) {
				neighbors.push_back(neighbor);
				set_1.add(neighbor);
			}

			// check if clique
			isolated.clear();
			isolated.push_back(node);
			non_isolated.clear();

			size_t max_isolated_idx = 0;
			weighted_node max_isolated{ node, status.weights[node] };
			weighted_node max_non_isolated{ 0, 0 };

			bool is_clique = true;

			for (auto neighbor : neighbors) {
				size_t count = 0;
				bool is_isolated = true;

				for (NodeID neighbor_2nd : status.graph[neighbor]) {
					if (set_1.get(neighbor_2nd)) count++;
					else is_isolated = false;
				}

				if (is_isolated) {
					isolated.push_back(neighbor);
					if (status.weights[neighbor] > max_isolated.weight) {
						max_isolated = { neighbor, status.weights[neighbor] };
						max_isolated_idx = isolated.size() - 1;
					}
				}
				else {
					non_isolated.push_back(neighbor);
					if (status.weights[neighbor] > max_non_isolated.weight) {
						max_non_isolated = { neighbor, status.weights[neighbor] };
					}
				}

				is_clique = count == neighbors.size();
				if (!is_clique) break;
			}

			if (!is_clique) continue;

			// one of "isolated" members has highest weight of clique: Add to IS
			// also handles completely isolated cliques
			if (max_isolated.weight >= max_non_isolated.weight) {
				br_alg->set(max_isolated.node, IS_status::included);
				continue;
			}


			// remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -> we can always pick "max_isolated" over them
			isolated[max_isolated_idx] = isolated.back();
			isolated.pop_back();

			for (auto neighbor : isolated) {
				br_alg->set(neighbor, IS_status::excluded);
			}

			for (size_t i = 0; i < non_isolated.size(); i++) {
				NodeID neighbor = non_isolated[i];
				if (status.weights[neighbor] <= max_isolated.weight) {
					br_alg->set(neighbor, IS_status::excluded);
					non_isolated[i] = non_isolated.back();
					non_isolated.pop_back();
					i--;
				}
			}

			fold(br_alg, std::move(max_isolated), std::move(non_isolated));
		}
	}

	//std::cout << "clique redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void clique_reduction::fold(branch_and_reduce_algorithm* br_alg, const weighted_node& isolated, std::vector<NodeID>&& non_isolated) {
	auto& status = br_alg->status;

	br_alg->set(isolated.node, IS_status::folded);
	status.reduction_offset += isolated.weight;

	for (auto node : non_isolated) {
		status.weights[node] -= isolated.weight;
		br_alg->add_next_level_neighborhood(node);
	}

	status.folded_queue.push_back(get_reduction_type());
	br_alg->add_next_level_neighborhood(non_isolated);

	restore_vec.emplace_back(isolated, std::move(non_isolated));
}

void clique_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.isolated.node);
	status.reduction_offset -= data.isolated.weight;

	for (auto node : data.non_isolated) {
		status.weights[node] += data.isolated.weight;
	}

	restore_vec.pop_back();
}

void clique_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto isolated = restore_vec.back().isolated.node;

	bool set_isolated = true;

	for (auto node : restore_vec.back().non_isolated) {
		if (status.node_status[node] == IS_status::included) {
			set_isolated = false;
			break;
		}
	}

	status.is_weight += restore_vec.back().isolated.weight;

	restore(br_alg);

	if (set_isolated) {
		status.node_status[isolated] = IS_status::included;
	}
	else {
		status.node_status[isolated] = IS_status::excluded;
	}
}

bool twin_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& twin_candidates_set = br_alg->set_1;
	auto& tmp_set = br_alg->set_2;
	size_t oldn = status.remaining_nodes;

	NodeID twin;
	NodeWeight neighbors_weight;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			neighbors.clear();
			neighbors_weight = 0;

			for (NodeID neighbor : status.graph[v]) {
				neighbors.push_back(neighbor);
				neighbors_weight += status.weights[neighbor];
			}

			if (status.weights[v] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				continue;
			}

			twin_candidates_set.clear();
			bool candidates_empty = true;

			for (NodeID neighbor : status.graph[neighbors[0]]) {
				if (neighbor != v && br_alg->deg(neighbor) == neighbors.size()) {
					twin_candidates_set.add(neighbor);
					candidates_empty = false;
					twin = neighbor;
				}
			}

			for (size_t i = 1; i < neighbors.size() && !candidates_empty; i++) {
				NodeID neighbor = neighbors[i];
				tmp_set.clear();
				candidates_empty = true;

				for (NodeID candidate : status.graph[neighbor]) {
					if (twin_candidates_set.get(candidate)) {
						tmp_set.add(candidate);
						candidates_empty = false;
						twin = candidate;
					}
				}

				std::swap(twin_candidates_set, tmp_set);
			}

			if (candidates_empty)
				continue;

			if (status.weights[v] + status.weights[twin] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				br_alg->set(twin, IS_status::included);
			} else {
				fold(br_alg, v, twin);
			}
		}
	}

	//std::cout << "twin redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void twin_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin) {
	auto& status = br_alg->status;

	restore_vec.push_back({ main, twin });

	br_alg->set(twin, IS_status::folded, true);
	status.weights[main] += status.weights[twin];

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(main);
	br_alg->add_next_level_neighborhood(main);
}

void twin_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.twin);
	status.weights[data.main] -= status.weights[data.twin];

	restore_vec.pop_back();
}

void twin_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto main = restore_vec.back().main;
	auto twin = restore_vec.back().twin;

	restore(br_alg);

	if (status.node_status[main] == IS_status::included) {
		status.node_status[twin] = IS_status::included;
	} else {
		status.node_status[twin] = IS_status::excluded;
	}
}

bool domination_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			NodeWeight neighbors_weight = 0;
			size_t neighbors_count = 0;
			neighbors.clear();

			for (NodeID neighbor : status.graph[v]) {
				neighbors.add(neighbor);
				neighbors_weight += status.weights[neighbor];
				neighbors_count++;
			}

			if (status.weights[v] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				continue;
			}

			neighbors.add(v);
			bool is_subset;

			for (NodeID neighbor : status.graph[v]) {
				if (br_alg->deg(neighbor) > neighbors_count)
					continue;

				is_subset = true;

				for (NodeID neighbor2 : status.graph[neighbor]) {
					if (!neighbors.get(neighbor2)) {
						is_subset = false;
						break;
					}
				}

				if (is_subset && status.weights[neighbor] >= status.weights[v]) {
					br_alg->set(v, IS_status::excluded);
					break;
				}
			}
		}
	}

	//std::cout << "domination redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool generalized_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;
	auto config = br_alg->config;

	cout_handler::disable_cout();

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			NodeWeight neighbors_weight = 0;
			NodeWeight max_neighbor_weight = 0;

			for (NodeID neighbor : status.graph[v]) {
				neighbors_weight += status.weights[neighbor];

				if (status.weights[neighbor] > max_neighbor_weight)
					max_neighbor_weight = status.weights[neighbor];
			}

			if (status.weights[v] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				continue;
			}
			if (status.weights[v] < max_neighbor_weight) {
				//MWIS in N(v) >= max_neighbor_weight > w(v)
				continue;
			}

			if (status.graph[v].size() > branch_and_reduce_algorithm::REDU_RECURSION_LIMIT) {
				continue;
			}

			//std::cerr << "%gnr try neighborhood of size " << status.graph[v].size() << " start node: " << v << std::endl;

			// compute MWIS in N(v)
			config.time_limit = status.graph[v].size() / 10.0;
			br_alg->build_induced_neighborhood_subgraph(neighborhood_graph, v);
			branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

			if (!neighborhood_br_alg.run_branch_reduce()) {
				std::cerr << "%generalized_neighborhood_reduction br_call time out" << std::endl;
				continue;
			}

			if (status.weights[v] >= neighborhood_br_alg.get_current_is_weight())
				br_alg->set(v, IS_status::included);
		}
	}

	cout_handler::enable_cout();

	//std::cout << "generalized_neighborhood_reduction -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool generalized_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighbors_set = br_alg->set_1;
	auto& MWIS_set = br_alg->set_2;
	auto& reverse_mapping = br_alg->buffers[1];
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;
	auto config = br_alg->config;

	cout_handler::disable_cout();

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			neighbors.clear();
			neighbors_set.clear();
			NodeWeight neighbors_weight = 0;
			NodeWeight max_neighbor_weight = 0;

			for (NodeID neighbor : status.graph[v]) {
				neighbors_weight += status.weights[neighbor];
				neighbors.push_back(neighbor);
				neighbors_set.add(neighbor);

				if (status.weights[neighbor] > max_neighbor_weight)
					max_neighbor_weight = status.weights[neighbor];
			}

			if (status.weights[v] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				continue;
			}

			if (status.graph[v].size() > branch_and_reduce_algorithm::REDU_RECURSION_LIMIT) {
				continue;
			}

			//std::cerr << "%gfr try neighborhood of size " << status.graph[v].size() << " start node: " << v << std::endl;

			// compute MWIS in N(v)
			config.time_limit = status.graph[v].size() / 10.0;

			br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);
			branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

			if (!neighborhood_br_alg.run_branch_reduce()) {
				std::cerr << "%generalized_fold_reduction br_call time out" << std::endl;
				continue;
			}

			NodeWeight MWIS_weight = neighborhood_br_alg.get_current_is_weight();
			NodeWeight min_MWIS_neighbor_weight = std::numeric_limits<NodeWeight>::max();

			if (status.weights[v] >= MWIS_weight) {
				// same as in generalized_neighborhood_reduction
				br_alg->set(v, IS_status::included);
				continue;
			}

			neighborhood_br_alg.apply_branch_reduce_solution(neighborhood_graph);
			MWIS_set.clear();

			forall_nodes(neighborhood_graph, node) {
				if (neighborhood_graph.getPartitionIndex(node) == 1) {
					const NodeID neighbor = neighbors[node];
					MWIS_set.add(neighbor);

					if (status.weights[neighbor] < min_MWIS_neighbor_weight)
						min_MWIS_neighbor_weight = status.weights[neighbor];
				}
			} endfor

			if (status.weights[v] < MWIS_weight - min_MWIS_neighbor_weight) {
				// multiple IS exist that have bigger weight than v
				continue;
			}

			bool check_failed = false;

			// check that no other IS in N(v) exists with weight greater than v
			for (const NodeID neighbor : status.graph[v]) {
				if (!MWIS_set.get(neighbor))
					continue;

				neighbors.remove(std::find(neighbors.begin(), neighbors.end(), neighbor));
				neighbors_set.remove(neighbor);

				br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);
				branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

				if (!neighborhood_br_alg.run_branch_reduce()) {
					std::cerr << "%generalized_fold_reduction br_call loop time out" << std::endl;
					check_failed = true;
				}
				else if (neighborhood_br_alg.get_current_is_weight() >= status.weights[v]) {
					check_failed = true;
				}

				neighbors.push_back(neighbor);
				neighbors_set.add(neighbor);

				if (check_failed)
					break;
			}

			if (!check_failed) {
				fold(br_alg, v, MWIS_set, MWIS_weight);
				continue;
			}

			auto& neighborhood_intersection_set = MWIS_set;
			bool remove_node;

			// we can't fold but we can possibly remove some neighbors of v
			do {
				for (const NodeID node : status.graph[v]) {
					neighborhood_intersection_set.clear();

					for (const NodeID neighbor : status.graph[node]) {
						if (neighbors_set.get(neighbor)) {
							neighborhood_intersection_set.add(neighbor);
						}
					}

					// "force" node into an IS (= remove node and its neighbors from N(v) and compute MWIS in remaining N(v))

					neighbors.remove(std::find(neighbors.begin(), neighbors.end(), node));
					neighbors_set.remove(node);

					for (const NodeID neighbor : status.graph[node]) {
						if (neighborhood_intersection_set.get(neighbor)) {
							neighbors.remove(std::find(neighbors.begin(), neighbors.end(), neighbor));
							neighbors_set.remove(neighbor);
						}
					}

					br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);

					config.time_limit = neighbors.size() / 10.0;
					branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

					if (!neighborhood_br_alg.run_branch_reduce()) {
						std::cerr << "%generalized_fold_reduction br_call loop time out" << std::endl;
						remove_node = false;
					}
					else {
						// if the weight of every MWIS in N(v) which contains "node" is smaller than w(v) then we can remove "node"
						remove_node = neighborhood_br_alg.get_current_is_weight() + status.weights[node] <= status.weights[v];
					}

					for (const NodeID neighbor : status.graph[node]) {
						if (neighborhood_intersection_set.get(neighbor)) {
							neighbors.push_back(neighbor);
							neighbors_set.add(neighbor);
						}
					}

					if (remove_node) {
						br_alg->set(node, IS_status::excluded);
						break; // break and restart loop because set(..) modifies the range which we currently iterate
					}

					neighbors.push_back(node);
					neighbors_set.add(node);
				}
			} while (remove_node);
		}
	}

	cout_handler::enable_cout();

	//if (oldn - status.remaining_nodes != 0)
		//std::cout << "%generalized_fold_reduction improvement: " << (oldn - status.remaining_nodes) << std::endl;

	//std::cout << "generalized_fold_reduction -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void generalized_fold_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main_node, fast_set& MWIS_set, NodeWeight MWIS_weight) {
	auto& status = br_alg->status;

	restore_vec.emplace_back();
	restore_data& data = restore_vec.back();
	data.main_weight = status.weights[main_node];
	data.MWIS_weight = MWIS_weight;

	auto& nodes = data.nodes;
	nodes.main = main_node;

	// temporary copy for iteration
	data.main_neighbor_list = status.graph[main_node];

	for (auto neighbor : data.main_neighbor_list) {
		if (MWIS_set.get(neighbor))
			nodes.MWIS.push_back(neighbor);
		else
			br_alg->set(neighbor, IS_status::excluded);
	}

	// reverse order because of later "restore_edge_and_replace"
	for (int i = nodes.MWIS.size() - 1; i >= 1; i--) {
		br_alg->set(nodes.MWIS[i], IS_status::folded, false);
	}

	br_alg->set(nodes.MWIS[0], IS_status::folded, true);

	data.main_neighbor_list = status.graph[main_node];

	// "move" weight into redu offset
	status.reduction_offset += data.main_weight;

	status.weights[nodes.main] = MWIS_weight - data.main_weight;

	std::vector<NodeID> new_neighbors;
	auto& neighbors = MWIS_set;
	neighbors.clear();
	neighbors.add(main_node);

	for (NodeID MWIS_node : nodes.MWIS) {
		std::vector<NodeID> node_vec;

		for (auto neighbor : status.graph[MWIS_node]) {
			if (neighbors.add(neighbor)) {
				new_neighbors.push_back(neighbor);
				status.graph.restore_edge_and_replace(neighbor, nodes.main);
				node_vec.push_back(neighbor);
			}
		}

		data.MWIS_node_vecs.push_back(std::move(node_vec));
	}

	status.graph[nodes.main] = dynamic_graph::neighbor_list(std::move(new_neighbors));
	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.main);
	br_alg->add_next_level_neighborhood(nodes.main);
}

void generalized_fold_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	// is "restored" in following loops
	status.graph.hide_node(data.nodes.main);
	status.graph[data.nodes.main] = std::move(data.main_neighbor_list);

	for (size_t i = 0; i < data.nodes.MWIS.size(); i++) {
		br_alg->unset(data.nodes.MWIS[i]);

		for (auto neighbor : data.MWIS_node_vecs[i]) {
			status.graph.replace_last_restored_edge(neighbor, data.nodes.MWIS[i]);
		}
	}

	status.weights[data.nodes.main] = data.main_weight;
	status.reduction_offset -= data.main_weight;

	restore_vec.pop_back();
}

void generalized_fold_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
	auto MWIS_weight = restore_vec.back().MWIS_weight;
	auto main_status = status.node_status[nodes.main];
	restore(br_alg);

	if (main_status == IS_status::included) {
		status.node_status[nodes.main] = IS_status::excluded;

		for (auto node : nodes.MWIS) {
			status.node_status[node] = IS_status::included;
		}

		status.is_weight += MWIS_weight;
	} else {
		status.node_status[nodes.main] = IS_status::included;

		for (auto node : nodes.MWIS) {
			status.node_status[node] = IS_status::excluded;
		}

		status.is_weight += status.weights[nodes.main];
	}
}
