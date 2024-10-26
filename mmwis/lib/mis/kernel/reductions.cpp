/******************************************************************************
* reductions.cpp
*****************************************************************************/

#include "reductions.h"
#include "branch_and_reduce_algorithm.h"
#include "data_structure/flow_graph.h"
#include "push_relabel.h"
#include "definitions.h"

#include <cstddef>
#include <utility>

#include <iostream>

using namespace mmwis;

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_neighborhood) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			NodeWeight neighbor_weights = br_alg->get_unset_neighbor_weight_sum(v, status);

			if (status.weights[v] >= neighbor_weights) 
				br_alg->set(v, IS_status::included);
		}
	}

	// std::cout << "neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.neighborhood_reduced_nodes+= (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_clique_neighborhood_fast) return false;

	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighborhood = br_alg->set_1;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

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

				if (!neighborhood.get(neighbor1)) continue;

				for (NodeID neighbor2 : status.graph[neighbor1]) {
					if (neighbor2 != neighbor1 && neighborhood.get(neighbor2)) {

						// triangle [v, neighbor1, neighbor2] found
						neighbor_weights -= std::min(status.weights[neighbor1], status.weights[neighbor2]);
						neighborhood.remove(neighbor1);
						neighborhood.remove(neighbor2);

						if (status.weights[v] >= neighbor_weights) 
							is_reducible = true;

						break;
					}
				}
			}

			if (is_reducible) br_alg->set(v, IS_status::included);
		}
	}

	// std::cout << "clique neighbor fast redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.clique_neighborhood_fast_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_clique_neighborhood) return false;

	this->br_alg = br_alg;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set && partition_into_cliques(v)) 
			br_alg->set(v, IS_status::included);
	}

	// std::cout << "clique neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.clique_neighborhood_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::partition_into_cliques(NodeID v) {
	auto& status = br_alg->status;
	auto& neighbors_vec = br_alg->buffers[0];
	auto& clique_neighbors_set = br_alg->double_set;

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

		if (intersection_empty) break;

		// add local_max to current clique
		neighbor_weights -= local_max_weight;
		if (neighbor_weights <= target_weight) return true;

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
	auto config = br_alg->config;
    if (config.disable_critical_set) return false;

	auto& status = br_alg->status;
	size_t n = br_alg->status.n;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

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

	// std::cout << "cs redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.critical_set_reduced_nodes += oldn-status.remaining_nodes; */
	return oldn != status.remaining_nodes;
}


bool fold1_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_fold1) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);
		if (status.node_status[v] == IS_status::not_set && br_alg->deg(v) == 0) 
		{
            br_alg->set(v, IS_status::included);
		} else if (status.node_status[v] == IS_status::not_set && br_alg->deg(v) == 1) {
            NodeID neighbor = status.graph[v][0];
            if (status.weights[neighbor] <= status.weights[v]) {
                br_alg->set(v, IS_status::included, true);
            } else {
                fold(br_alg, {v, neighbor});
            }
		}
	}

	// std::cout << "fold1 redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.fold1_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

void fold1_reduction::fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {

    auto& status = br_alg->status;

	restore_vec.push_back({nodes, status.weights[nodes.deg1_node], status.graph[nodes.fold_node]});
    br_alg->set(nodes.deg1_node, IS_status::folded);

    status.reduction_offset += status.weights[nodes.deg1_node];
    status.weights[nodes.fold_node] -= status.weights[nodes.deg1_node];

    status.folded_queue.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.fold_node);
    br_alg->add_next_level_neighborhood(nodes.fold_node);
}

void fold1_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

    br_alg->unset(data.nodes.deg1_node);

	status.weights[data.nodes.fold_node] += data.deg1_weight;
	status.reduction_offset -= data.deg1_weight; 

	restore_vec.pop_back();
}

void fold1_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
	auto main_status = status.node_status[nodes.fold_node];
	restore(br_alg);

	if (main_status == IS_status::included) {
		status.node_status[nodes.fold_node] = IS_status::included;
		status.node_status[nodes.deg1_node] = IS_status::excluded;

		status.is_weight += status.weights[nodes.fold_node];

	} else if (main_status == IS_status::excluded) {
		status.node_status[nodes.fold_node] = IS_status::excluded;
		status.node_status[nodes.deg1_node] = IS_status::included;

		status.is_weight += status.weights[nodes.deg1_node];
	}
}

bool clique_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_clique) return false;

	auto& status = br_alg->status;
	auto& set_1 = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	auto& isolated = br_alg->buffers[1];
	std::vector<NodeID> non_isolated;

	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

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


   // remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -->we can always pick "max_isolated" over them
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

	// std::cout << "clique redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.clique_reduced_nodes += (oldn - status.remaining_nodes); */
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

bool triangle_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_triangle) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set && br_alg->deg(v) == 2) {
			bool triangle = false;

            NodeID bigger  = status.graph[v][0];
            NodeID smaller = status.graph[v][1];

            if (status.weights[bigger] < status.weights[smaller]) {
                smaller = status.graph[v][0];
                bigger  = status.graph[v][1];
            }

            for (NodeID neighbor : status.graph[smaller]) {
                if (neighbor == bigger) 
                    triangle = true;
            }
            if (!triangle) continue;
            
            if (status.weights[bigger] <= status.weights[v]) {
            // triangle reduction case 1
                  br_alg->set(v, IS_status::included);
            } else if (status.weights[bigger] > status.weights[v] && status.weights[smaller] <= status.weights[v]) {
            // triangle reduction case 2
                fold_mid_weight(br_alg, {v, bigger, smaller});
            } else {
            // triangle reduction case 3
                fold_min_weight(br_alg, {v, bigger, smaller});
            }
		}
	}

	// std::cout << "triangle redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.triangle_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

void triangle_reduction::fold_mid_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;

    br_alg->set(nodes.deg2_node, IS_status::folded);
	br_alg->set(nodes.smaller, IS_status::excluded);

    restore_vec.push_back({nodes, status.weights[nodes.deg2_node], 1});
    
	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.bigger] -= status.weights[nodes.deg2_node];

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.bigger);
	br_alg->add_next_level_neighborhood(nodes.bigger);

}

void triangle_reduction::fold_min_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;

    br_alg->set(nodes.deg2_node, IS_status::folded);

    restore_vec.push_back({nodes, status.weights[nodes.deg2_node], 0});
    
	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.bigger] -= status.weights[nodes.deg2_node];
	status.weights[nodes.smaller] -= status.weights[nodes.deg2_node];

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.bigger);
	br_alg->add_next_level_neighborhood(nodes.bigger);

}

void triangle_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();
    br_alg->unset(data.nodes.deg2_node);
	status.reduction_offset -= data.deg2_weight;

    if (data.fold_case == 1) {  // mid
	    status.weights[data.nodes.bigger] += data.deg2_weight;
    }
    else {  // min
	    status.weights[data.nodes.bigger] += data.deg2_weight;
	    status.weights[data.nodes.smaller] += data.deg2_weight;
    }

	restore_vec.pop_back();
}

void triangle_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
	auto bigger_status = status.node_status[nodes.bigger];
	auto smaller_status = status.node_status[nodes.smaller];
	restore(br_alg);
    if (bigger_status == IS_status::included) {
    	status.node_status[nodes.bigger] = IS_status::included;
    	status.node_status[nodes.smaller] = IS_status::excluded;
    	status.node_status[nodes.deg2_node] = IS_status::excluded;
    	status.is_weight += status.weights[nodes.bigger];

    } else if (smaller_status == IS_status::included){
    	status.node_status[nodes.bigger] = IS_status::excluded;
    	status.node_status[nodes.smaller] = IS_status::included;
    	status.node_status[nodes.deg2_node] = IS_status::excluded;
    	status.is_weight += status.weights[nodes.smaller];

    } else {
    	status.node_status[nodes.bigger] = IS_status::excluded;
    	status.node_status[nodes.smaller] = IS_status::excluded;
    	status.node_status[nodes.deg2_node] = IS_status::included;
    	status.is_weight += status.weights[nodes.deg2_node];
    }
}

bool v_shape_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_v_shape_mid && config.disable_v_shape_max) return false;

    auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

    for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set && br_alg->deg(v) == 2) {

			bool triangle = false;
            NodeID bigger  = status.graph[v][0];
            NodeID smaller = status.graph[v][1];

            if (status.weights[bigger] < status.weights[smaller]) {
                smaller = status.graph[v][0];
                bigger  = status.graph[v][1];
            }

            for (NodeID neighbor : status.graph[smaller]) {
                if (neighbor == bigger) 
				{
					triangle = true;
					break;
				} 
            }

            if (triangle) continue;

            if (status.weights[v] >= status.weights[bigger]) {
                if (config.disable_v_shape_max) continue;
                if (status.weights[v] >= status.weights[bigger] + status.weights[smaller]) {
                    // v is always best choice of the three vertices
                    // auto current_n = status.remaining_nodes;
                    br_alg->set(v, IS_status::included);
                    continue;
                }
                fold_max_weight(br_alg, {v, {bigger, smaller}});
            } else {
                if (status.weights[v] >= status.weights[smaller]) {
                    if (config.disable_v_shape_mid) continue;
                    fold_mid_weight(br_alg, {v, {bigger, smaller}});
                }
            }
        }
    }

    // std::cout << "v_shape redu -> " << (oldn - status.remaining_nodes) << std::endl; 

	return oldn != status.remaining_nodes;
}

bool v_shape_min_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_v_shape_min) return false;

    auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

    for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set && br_alg->deg(v) == 2) {

			bool triangle = false;
            NodeID bigger  = status.graph[v][0];
            NodeID smaller = status.graph[v][1];

            if (status.weights[bigger] < status.weights[smaller]) {
                smaller = status.graph[v][0];
                bigger  = status.graph[v][1];
            }

            for (NodeID neighbor : status.graph[smaller]) {
                if (neighbor == bigger) {
                    triangle = true;
					break;
				}
            }

            if (triangle) continue;

            if (status.weights[v] < status.weights[smaller]) 
                fold(br_alg, {v, {bigger, smaller}});
        }
    }

    // std::cout << "v_shape_min redu weight -> " << (status.reduction_offset - oldw) << std::endl; 
/* 	status.v_shape_min_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldw != status.reduction_offset;
}

void v_shape_min_reduction::fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;

	status.v_shape_min_applied++; 
    if (status.v_shape_min_applied > status.n*status.m) {
       status.m++;
       status.modified_queue.resize((status.m+1)*status.n);
       status.folded_queue.resize((status.m+1)*status.n);
       status.branching_queue.resize((status.m+1)*status.n);
    }

    // push deg2_node to modified_queue
    status.changed_status[nodes.deg2_node] += 1;
    status.modified_queue.push_back(nodes.deg2_node);

	restore_vec.push_back({ nodes, status.weights[nodes.deg2_node], {}});

	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.neighbors[0]] = status.weights[nodes.neighbors[0]] - status.weights[nodes.deg2_node];
	status.weights[nodes.neighbors[1]] = status.weights[nodes.neighbors[1]] - status.weights[nodes.deg2_node];

    status.graph.hide_edge_undirected(nodes.deg2_node, nodes.neighbors[0]);
    status.graph.hide_edge_undirected(nodes.deg2_node, nodes.neighbors[1]);

	neighbors.clear();
    neighbors.add(nodes.deg2_node);

	for (size_t i = 0; i < 2; i++) {
		for (auto neighbor : status.graph[nodes.neighbors[i]]) {
			if (neighbors.add(neighbor)) {
               status.graph.add_edge_undirected(neighbor, nodes.deg2_node, true);
               restore_vec.back().node_vecs[i].push_back(neighbor);
			}
		}
	}

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.deg2_node);
	br_alg->add_next_level_node(nodes.neighbors[0]);
	br_alg->add_next_level_node(nodes.neighbors[1]);
	br_alg->add_next_level_neighborhood(nodes.neighbors[0]);
	br_alg->add_next_level_neighborhood(nodes.neighbors[1]);

}


void v_shape_reduction::fold_mid_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
    auto oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;
    NodeID bigger  = nodes.neighbors[0];
    NodeID smaller = nodes.neighbors[1];
    NodeWeight deg2_weight = status.weights[nodes.deg2_node];

	restore_vec.push_back({nodes, deg2_weight, 2, status.graph[smaller], {}});
	br_alg->set(nodes.deg2_node, IS_status::folded);

	status.reduction_offset += deg2_weight;
	status.weights[bigger]  -= deg2_weight;
	neighbors.clear();
	neighbors.add(nodes.deg2_node);
	neighbors.add(smaller);
	neighbors.add(bigger);

    std::vector<NodeID> new_neighbors;

    for (auto neighbor : status.graph[smaller]) {
        neighbors.add(neighbor);
    }

    for (auto neighbor : status.graph[bigger]) {
        if (neighbors.add(neighbor)) {
            status.graph.add_edge_undirected(smaller, neighbor, true);
            restore_vec.back().node_vecs[1].push_back(neighbor);
        }
    }

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(smaller);
	br_alg->add_next_level_node(bigger);
	br_alg->add_next_level_neighborhood(smaller);

/* 	status.v_shape_mid_reduced_nodes += (oldn - status.remaining_nodes); */
}


void v_shape_reduction::fold_max_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
    auto oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	restore_vec.push_back({ nodes, status.weights[nodes.deg2_node], 1, status.graph[nodes.deg2_node],{}});

	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.deg2_node] = status.weights[nodes.neighbors[0]] + status.weights[nodes.neighbors[1]] - status.weights[nodes.deg2_node];

	neighbors.clear();
	neighbors.add(nodes.deg2_node);

	for (size_t i = 0; i < 2; i++) {
		for (auto neighbor : status.graph[nodes.neighbors[i]]) {
			if (neighbors.add(neighbor)) {
                status.graph.add_edge_undirected(nodes.deg2_node, neighbor, true);
			    restore_vec.back().node_vecs[i].push_back(neighbor);
			}
		}
	}
	br_alg->set(nodes.neighbors[1], IS_status::folded,true);
	br_alg->set(nodes.neighbors[0], IS_status::folded,false);

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.neighbors[0]);
/* 	status.v_shape_max_reduced_nodes += (oldn - status.remaining_nodes); */
}

void v_shape_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();
    if (data.fold_case == 1) {  // max
    
		br_alg->unset(data.nodes.neighbors[0]);
		br_alg->unset(data.nodes.neighbors[1]);

        for (size_t i = 0; i < 2; i++) {
            for (NodeID neighbor : data.node_vecs[i]) {
                status.graph.remove_edge_undirected(data.nodes.deg2_node, neighbor);
            }
        }
 
 	   status.weights[data.nodes.deg2_node] = data.deg2_weight;
 	   status.reduction_offset -= data.deg2_weight;
 	   restore_vec.pop_back();
    }

    if (data.fold_case == 2) { // mid

        br_alg->unset(data.nodes.deg2_node);
		if (status.node_status[data.nodes.neighbors[0]] != IS_status::not_set)
        	br_alg->unset(data.nodes.neighbors[0], false);
		if (status.node_status[data.nodes.neighbors[1]] != IS_status::not_set)
	        br_alg->unset(data.nodes.neighbors[1], false);

        // restore neighbors from smaller node
        for (auto neighbor : data.node_vecs[1]) {
            status.graph.remove_edge_undirected(neighbor, data.nodes.neighbors[1]);
        }

        status.weights[data.nodes.neighbors[0]] += data.deg2_weight;
        status.reduction_offset -= data.deg2_weight;

	    restore_vec.pop_back();
    }
}

void v_shape_min_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();
    
    if (status.node_status[data.nodes.deg2_node] != IS_status::not_set)
	    br_alg->unset(data.nodes.deg2_node, false);

    for (size_t i = 0; i < 2; i++) {
        if (status.node_status[data.nodes.neighbors[i]] != IS_status::not_set)
                br_alg->unset(data.nodes.neighbors[i], false);
    }

    status.changed_status[data.nodes.deg2_node] -= 1;

	for (size_t i = 0; i < 2; i++) {
        for (auto second_neighbor : data.node_vecs[i]) {
                 status.graph.remove_edge_undirected(second_neighbor, data.nodes.deg2_node);
        }

        status.graph.restore_edge_and_replace(data.nodes.neighbors[i], data.nodes.deg2_node);
        status.graph.restore_edge_and_replace(data.nodes.deg2_node, data.nodes.neighbors[i]);
        status.weights[data.nodes.neighbors[i]] += data.deg2_weight;
	}

    status.reduction_offset -= data.deg2_weight;
    restore_vec.pop_back();
}

void v_shape_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
    int fold_case = restore_vec.back().fold_case;
    
	auto deg_2_status = status.node_status[nodes.deg2_node];
    auto nbh_0_status = status.node_status[nodes.neighbors[0]];
    auto nbh_1_status = status.node_status[nodes.neighbors[1]];
	restore(br_alg);

    if (fold_case == 1){ //max 
        if (deg_2_status == IS_status::included) {
            status.node_status[nodes.deg2_node]    = IS_status::excluded;
            status.node_status[nodes.neighbors[0]] = IS_status::included;
            status.node_status[nodes.neighbors[1]] = IS_status::included;

            status.is_weight += status.weights[nodes.neighbors[0]]; 
            status.is_weight += status.weights[nodes.neighbors[1]];

        } else {
        // (deg_2_status == IS_status::excluded)
            status.node_status[nodes.deg2_node]    = IS_status::included;
            status.node_status[nodes.neighbors[0]] = IS_status::excluded;
            status.node_status[nodes.neighbors[1]] = IS_status::excluded;

            status.is_weight += status.weights[nodes.deg2_node];
        }
    }


    if (fold_case >= 2){ //mid
        status.node_status[nodes.deg2_node]    = IS_status::excluded;
        status.node_status[nodes.neighbors[0]] = IS_status::excluded;
        status.node_status[nodes.neighbors[1]] = IS_status::excluded;

        if (nbh_1_status == IS_status::included) { // smaller included
            status.node_status[nodes.neighbors[0]] = IS_status::included;
            status.node_status[nodes.neighbors[1]] = IS_status::included;
            status.is_weight += status.weights[nodes.neighbors[0]];
            status.is_weight += status.weights[nodes.neighbors[1]];
        } else if (nbh_0_status == IS_status::included) {
            status.node_status[nodes.neighbors[0]] = IS_status::included;
            status.is_weight += status.weights[nodes.neighbors[0]];
        } else if (nbh_0_status == IS_status::excluded && nbh_1_status == IS_status::excluded) {
            status.node_status[nodes.deg2_node] = IS_status::included;
            status.is_weight += status.weights[nodes.deg2_node];
        } else if (nbh_0_status == IS_status::not_set || nbh_1_status == IS_status::not_set) {
            std::cerr << "ERROR restore mid nodes not set" << std::endl;
        }
    }
}

void v_shape_min_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
    
	auto deg_2_status = status.node_status[nodes.deg2_node];
    auto nbh_0_status = status.node_status[nodes.neighbors[0]];
    auto nbh_1_status = status.node_status[nodes.neighbors[1]];

	restore(br_alg);

	if (deg_2_status == IS_status::included) {
        status.node_status[nodes.deg2_node]    = IS_status::excluded;
        status.node_status[nodes.neighbors[0]] = IS_status::included;
        status.node_status[nodes.neighbors[1]] = IS_status::included;
        status.is_weight += status.weights[nodes.neighbors[0]];
        status.is_weight += status.weights[nodes.neighbors[1]];
	} else if (nbh_1_status == IS_status::included) {
        status.node_status[nodes.deg2_node]    = IS_status::excluded;
    	status.node_status[nodes.neighbors[1]] = IS_status::included;
   		status.is_weight += status.weights[nodes.neighbors[1]];
        if (nbh_0_status == IS_status::included) {
    	    status.node_status[nodes.neighbors[0]] = IS_status::included;
        	status.is_weight += status.weights[nodes.neighbors[0]];
    	} else if (nbh_0_status == IS_status::excluded){
    	    status.node_status[nodes.neighbors[0]] = IS_status::excluded;
    	}
    } else if (nbh_1_status == IS_status::excluded) {
        status.node_status[nodes.neighbors[1]] = IS_status::excluded;
        if (nbh_0_status == IS_status::included) {
           	status.node_status[nodes.deg2_node]    = IS_status::excluded;
    	    status.node_status[nodes.neighbors[0]] = IS_status::included;
            status.is_weight += status.weights[nodes.neighbors[0]];
        } else if (nbh_0_status == IS_status::excluded){
        	status.node_status[nodes.neighbors[0]] = IS_status::excluded;
	        status.node_status[nodes.deg2_node]    = IS_status::included;
	        status.is_weight += status.weights[nodes.deg2_node];
        }
    } else {
        std::cerr << "Error included v_shape_min_reduction::apply" <<  std::endl;

    }
}


bool single_edge_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_basic_se) return false;

	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;
    NodeWeight partial_neighbor_sum = 0;

	for (size_t node_idx = 0; node_idx < marker.current_size(); node_idx++) {
		NodeID node = marker.current_vertex(node_idx);
       

		if (status.node_status[node] == IS_status::not_set) {
            neighbors.clear();

            for (NodeID neighbor : status.graph[node]) {
                if (status.node_status[neighbor] == IS_status::not_set)
                    neighbors.add(neighbor);
            }

            for (NodeID neighbor : status.graph[node]) {
                if (status.weights[node] <= status.weights[neighbor]) { // otherwise not applicable to this edge
                                                                       //
                    // compute w(N(neighbor)\N(node))
                    partial_neighbor_sum = 0;
                    for (NodeID second_neighbor : status.graph[neighbor]) {
                        if (status.node_status[second_neighbor] == IS_status::not_set && neighbors.get(second_neighbor)) {
                            partial_neighbor_sum += status.weights[second_neighbor];
							if (partial_neighbor_sum >= status.weights[neighbor]) break;
                        }
                    }
                 
                    // note: weight of node is in partial_neighbor_sum included
                    if (partial_neighbor_sum <= status.weights[neighbor]) { 
                        br_alg->set(node, IS_status::excluded);
                        break;
                    }
                }
            }
        }
	}

	// std::cout << "single_edge redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.single_edge_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

bool extended_single_edge_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_extended_se) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;
    NodeWeight neighbor_sum = 0;
    NodeID max_neighbor = status.n;
    NodeWeight max_neighbor_weight = 0;
    bool skip = false;
	auto& neighbors = br_alg->set_1;
	auto& checked = br_alg->set_2;


	for (size_t node_idx = 0; node_idx < marker.current_size(); node_idx++) {
		NodeID node = marker.current_vertex(node_idx);

        if (status.node_status[node] == IS_status::not_set) {
            checked.clear();
            neighbors.clear();

            // compute w(N(v)) and max(N(v))
			neighbor_sum = br_alg->get_unset_neighbor_weight_sum_and_max(node, max_neighbor, status);
			max_neighbor_weight = status.weights[max_neighbor];

            if (status.weights[node] >= neighbor_sum) {
                br_alg->set(node, IS_status::included);
                continue;
            }

            for (NodeID neighbor : status.graph[node]) {
                if (status.node_status[neighbor] == IS_status::not_set) {
                    neighbors.add(neighbor);
                }
            }
            
            while (status.weights[node] >= neighbor_sum - max_neighbor_weight) {

                if (status.weights[node] >= neighbor_sum) {
                    br_alg->set(node, IS_status::included);
                    break;
                }
            
                for (NodeID neighbor : status.graph[max_neighbor]) {
                    if (neighbor == node)  continue; 
                    if (status.node_status[neighbor] == IS_status::not_set) {
                        // exclude neighborhood intersection and update neighborhood
                        if (neighbors.get(neighbor)) { 
                            br_alg->set(neighbor, IS_status::excluded);
                            neighbors.remove(neighbor);
                            neighbor_sum -= status.weights[neighbor];
                        }
                    }
                }


                // check if different edge satisfies reduction
                checked.add(max_neighbor);
                max_neighbor_weight = 0;
                for (NodeID neighbor : status.graph[node]) {
                    if (checked.get(neighbor)) continue;
                    if (max_neighbor_weight < status.weights[neighbor]) { 
                        max_neighbor_weight = status.weights[neighbor];
                        max_neighbor = neighbor;
                    }
                }
            }

        }
	}

	// std::cout << "extended single_edge redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.extended_single_edge_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}



bool twin_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_twin) return false;

	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& twin_candidates_set = br_alg->set_1;
	auto& tmp_set = br_alg->set_2;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

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

	// std::cout << "twin redu -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.twin_reduced_nodes += (oldn - status.remaining_nodes); */
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

bool generalized_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_generalized_neighborhood) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	graph_access neighborhood_graph;


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


	// std::cout << "generalized_neighborhood_reduction -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.generalized_neighborhood_reduced_nodes += (oldn - status.remaining_nodes); */
	return oldn != status.remaining_nodes;
}

bool heavy_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_heavy_set == 0) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	auto& build_graph_neighbors = br_alg->buffers[0];
	auto& reverse_mapping = br_alg->buffers[1];

	//sets to build different subgraphs
	auto& build_graph_neighbors_set = br_alg->set_1;
	auto& heavy_vertex_neighbors_set = br_alg->set_2;
    fast_set second_heavy_vertex_neighbors_set(status.n);

	graph_access neighborhood_graph;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {

            NodeID heavy_vertex = 0;
            NodeID second_heavy_vertex = 0;
			::NodeWeight neighbors_weight = 0;
			::NodeWeight second_heavy_vertex_weight = 0;
			::NodeWeight heavy_vertex_weight = 0;

			neighbors_weight = br_alg->get_unset_neighbor_weight_sum_and_max(v, heavy_vertex, status);

			if (status.weights[v] >= neighbors_weight) {
                br_alg->set(v, IS_status::included);
				continue;
			}

			for (NodeID neighbor : status.graph[v]) {
                if (br_alg->deg(neighbor) > config.disable_heavy_set-1) continue; //subgraph too large
			}

            if (heavy_vertex_weight==0) continue; //no neighbor with small enough degree

			build_graph_neighbors.clear();
			build_graph_neighbors_set.clear();
			heavy_vertex_neighbors_set.clear();

			for (NodeID neighbor : status.graph[heavy_vertex]) {
		        if (status.node_status[neighbor] == IS_status::not_set)
				{
				    build_graph_neighbors.push_back(neighbor);
				    build_graph_neighbors_set.add(neighbor);
				    heavy_vertex_neighbors_set.add(neighbor);
                }
            }

			//find second heavy vertex (not adjacent to heavy vertex)
			second_heavy_vertex_weight = 0;
			for (NodeID neighbor : status.graph[v]) {
		        if (status.node_status[neighbor] == IS_status::not_set) {
                    if (neighbor == heavy_vertex) continue; // look for different nodes
                    if (heavy_vertex_neighbors_set.get(neighbor)) continue; // look for non adjacent nodes
                    if (br_alg->deg(neighbor) + br_alg->deg(heavy_vertex) > config.disable_heavy_set) continue; //subgraph too large

				    if (status.weights[neighbor] > second_heavy_vertex_weight) {
				    	second_heavy_vertex_weight = status.weights[neighbor];
                        second_heavy_vertex = neighbor;
                    }
                }
            }

            if (second_heavy_vertex_weight == 0) continue;
			second_heavy_vertex_neighbors_set.clear();

			for (NodeID neighbor : status.graph[second_heavy_vertex])
			{
			 	if (status.node_status[neighbor] != IS_status::not_set) continue;
			    second_heavy_vertex_neighbors_set.add(neighbor);
			    build_graph_neighbors.push_back(neighbor);
			    build_graph_neighbors_set.add(neighbor);
            }

            //build neighborhood graph 
			config.time_limit = 8.0 / 10.0;

            //check different cases on the neighborhood graph:
			// case 1)
            // compute MWIS in N(heavy_vertex) \cup N(second_heavy_vertex):
			br_alg->build_induced_subgraph(neighborhood_graph, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);
			branch_and_reduce_algorithm neighborhood_br_alg_case1(neighborhood_graph, config, true);
            
			if (!neighborhood_br_alg_case1.run_branch_reduce()) {
				std::cerr << "heavy_vertex br_call time out" << std::endl;
				continue;
			}
            bool only_check_single = false;

			::NodeWeight MWIS_weight_case1 = neighborhood_br_alg_case1.get_current_is_weight();
            if (second_heavy_vertex_weight >= MWIS_weight_case1) {
                br_alg->set(heavy_vertex, IS_status::included);
                br_alg->set(second_heavy_vertex, IS_status::included);
                continue;
            } else if (heavy_vertex_weight + second_heavy_vertex_weight < MWIS_weight_case1)
                only_check_single = true;


			// case 2)
            // compute MWIS in N(heavy_vertex): 
            // first set graph to N(heavy_vertex)
			build_graph_neighbors.clear();
			std::for_each(status.graph[heavy_vertex].begin(), status.graph[heavy_vertex].end(), [&](NodeID neighbor) {
				if (status.node_status[neighbor] == IS_status::not_set) {
					build_graph_neighbors.push_back(neighbor);
				}
			});
			// TODO need clear graph?
			br_alg->build_induced_subgraph(neighborhood_graph, build_graph_neighbors, heavy_vertex_neighbors_set, reverse_mapping);

            //solve graph
			branch_and_reduce_algorithm neighborhood_br_alg_case2(neighborhood_graph, config, true);
			if (!neighborhood_br_alg_case2.run_branch_reduce()) {
				std::cerr << "heavy_vertex br_call time out" << std::endl;
				continue;
			}

			::NodeWeight MWIS_weight_case2 = neighborhood_br_alg_case2.get_current_is_weight();
            if (heavy_vertex_weight >= MWIS_weight_case2) {
                br_alg->set(heavy_vertex, IS_status::included);
                only_check_single = true;
            }
 
            if (!only_check_single) {
			    // case 3)
                // compute MWIS in N(heavy_vertex)\N(second_heavy_vertex): 

				build_graph_neighbors.clear();
				build_graph_neighbors_set.clear();
				std::for_each(status.graph[heavy_vertex].begin(), status.graph[heavy_vertex].end(), [&](NodeID neighbor) {
					if (status.node_status[neighbor] == IS_status::not_set && !second_heavy_vertex_neighbors_set.get(neighbor)) {
						build_graph_neighbors.push_back(neighbor);
						build_graph_neighbors_set.add(neighbor);
					}
				});
				br_alg->build_induced_subgraph(neighborhood_graph, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);

                //solve graph
			    branch_and_reduce_algorithm neighborhood_br_alg_case3(neighborhood_graph, config, true);
			    if (!neighborhood_br_alg_case3.run_branch_reduce()) {
			    	std::cerr << "heavy_vertex br_call time out" << std::endl;
			    	continue;
			    }
			    ::NodeWeight MWIS_weight_case3 = neighborhood_br_alg_case3.get_current_is_weight();

                    
			    // case 4)
                // compute MWIS in N(second_heavy_vertex)\N(heavy_vertex): 
                // first set graph 

				build_graph_neighbors.clear();
				build_graph_neighbors_set.clear();
				std::for_each(status.graph[second_heavy_vertex].begin(), status.graph[second_heavy_vertex].end(), [&](NodeID neighbor) {
					if (status.node_status[neighbor] == IS_status::not_set && !heavy_vertex_neighbors_set.get(neighbor)) {
						build_graph_neighbors.push_back(neighbor);
						build_graph_neighbors_set.add(neighbor);
					}
				});
				br_alg->build_induced_subgraph(neighborhood_graph, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);


                //solve graph
			    branch_and_reduce_algorithm neighborhood_br_alg_case4(neighborhood_graph, config, true);
			    if (!neighborhood_br_alg_case4.run_branch_reduce()) {
			    	std::cerr << "heavy_vertex br_call time out" << std::endl;
			    	continue;
			    }
			    ::NodeWeight MWIS_weight_case4 = neighborhood_br_alg_case4.get_current_is_weight();

                if (heavy_vertex_weight + second_heavy_vertex_weight >= MWIS_weight_case1 &&
                    heavy_vertex_weight >= MWIS_weight_case3 &&
                    second_heavy_vertex_weight >= MWIS_weight_case4) {
                    br_alg->set(heavy_vertex, IS_status::included);
                    br_alg->set(second_heavy_vertex, IS_status::included);
                    continue;
                }

            }

			// case 5)
            // compute MWIS in N(second_heavy_vertex):
            // first set graph 
            if (only_check_single) {
				build_graph_neighbors.clear();
				std::for_each(status.graph[second_heavy_vertex].begin(), status.graph[second_heavy_vertex].end(), [&](NodeID neighbor) {
					if (status.node_status[neighbor] == IS_status::not_set) {
						build_graph_neighbors.push_back(neighbor);
					}
				});
				br_alg->build_induced_subgraph(neighborhood_graph, build_graph_neighbors, second_heavy_vertex_neighbors_set, reverse_mapping);

	            //solve graph
				branch_and_reduce_algorithm neighborhood_br_alg_case5(neighborhood_graph, config, true);
				if (!neighborhood_br_alg_case5.run_branch_reduce()) {
					std::cerr << "heavy_vertex br_call time out" << std::endl;
					continue;
				}
				::NodeWeight MWIS_weight_case5 = neighborhood_br_alg_case5.get_current_is_weight();

	            if (second_heavy_vertex_weight >= MWIS_weight_case5) {
	                br_alg->set(second_heavy_vertex, IS_status::included);
	            }
            }
        }
    }
    /* cout_handler::enable_cout(); */

	// std::cout << "heavy_set_reduction -> " << (oldn - status.remaining_nodes) << std::endl;
    /* std::cout << "Remaining nodes: " <<status.remaining_nodes << std::endl; */
	return oldn != status.remaining_nodes;
}
// bool heavy_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
// 	auto config = br_alg->config;
//     if (config.disable_heavy_set == 0) return false;

// 	auto& status = br_alg->status;
// 	auto& neighbors = br_alg->buffers[0];
// 	auto& neighbors_set = br_alg->set_1;
// 	auto& MWIS_set = br_alg->set_2;
// 	auto& reverse_mapping = br_alg->buffers[1];
// 	size_t oldn = status.remaining_nodes;
// 	size_t oldw = status.reduction_offset;

// 	graph_access neighborhood_graph;


// 	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
// 		NodeID v = marker.current_vertex(v_idx);

// 		if (status.node_status[v] == IS_status::not_set) {
// 			neighbors.clear();
// 			neighbors_set.clear();
//             NodeID max_neighbor = 0;
//             NodeID second_max_neighbor = 0;
// 			NodeWeight heavy_set_weight = 0;
// 			NodeWeight neighbors_weight = 0;
// 			NodeWeight max_neighbor_weight = 0;
//             bool second_neighbor_found = false;

// 			for (NodeID neighbor : status.graph[v]) {
// 		        if (status.node_status[neighbor] == IS_status::not_set) {
//                     neighbors_weight += status.weights[neighbor];
//                     if (br_alg->deg(neighbor) > config.heavy_set) continue; //subgraph too large
// 				    if (status.weights[neighbor] > max_neighbor_weight) {
// 				    	max_neighbor_weight = status.weights[neighbor];
//                         max_neighbor = neighbor;
//                     }
//                 } 
// 			}

// 			if (status.weights[v] >= neighbors_weight) {
//                 br_alg->set(v, IS_status::included);
// 				continue;
// 			}

//             if (max_neighbor_weight==0) continue; //no neighbor with small enough degree

//             heavy_set_weight += max_neighbor_weight;
// 			for (NodeID neighbor : status.graph[max_neighbor]) {
// 		        if (status.node_status[neighbor] == IS_status::not_set) {
// 				    neighbors.push_back(neighbor);
// 				    neighbors_set.add(neighbor);
//                 }
//             }

// 			max_neighbor_weight = 0;
// 			for (NodeID neighbor : status.graph[v]) {
// 		        if (status.node_status[neighbor] == IS_status::not_set) {
//                     if (neighbor == max_neighbor) continue; // look for different nodes
//                     if (neighbors_set.get(neighbor)) continue; // look for non adjacent nodes
//                     if (br_alg->deg(neighbor) + br_alg->deg(max_neighbor) > config.heavy_set) continue; //subgraph too large
// 				    if (status.weights[neighbor] > max_neighbor_weight) {
// 				    	max_neighbor_weight = status.weights[neighbor];
//                         second_max_neighbor = neighbor;
//                         second_neighbor_found=true;
//                     }
//                 }
//             }

//             if (!second_neighbor_found) continue;
//             heavy_set_weight += max_neighbor_weight;

// 			for (NodeID neighbor : status.graph[second_max_neighbor]) {
// 		        if (status.node_status[neighbor] == IS_status::not_set) {
//                     if (neighbors_set.get(neighbor)) continue;
// 				    neighbors.push_back(neighbor);
// 				    neighbors_set.add(neighbor);
//                 }
//             }

            

// 			// compute MWIS in N(v)
// 			config.time_limit = 8.0 / 10.0;

// 			br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);

// 			branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

// 			if (!neighborhood_br_alg.run_branch_reduce()) {
// 				std::cerr << "heavy_set_reduction br_call time out" << std::endl;
// 				continue;
// 			}

// 			NodeWeight MWIS_weight = neighborhood_br_alg.get_current_is_weight();
//             if (heavy_set_weight >= MWIS_weight) {
//                 br_alg->set(max_neighbor, IS_status::included);
//                 br_alg->set(second_max_neighbor, IS_status::included);
// 				continue;
// 			}
//         }
//     }

// 	// std::cout << "heavy_set_reduction -> " << (oldn - status.remaining_nodes) << std::endl;
// 	return oldn != status.remaining_nodes;
// }

bool generalized_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto config = br_alg->config;
    if (config.disable_generalized_fold) return false;

	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighbors_set = br_alg->set_1;
	auto& MWIS_set = br_alg->set_2;
	auto& reverse_mapping = br_alg->buffers[1];
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	graph_access neighborhood_graph;


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


	// std::cout << "generalized_fold_reduction -> " << (oldn - status.remaining_nodes) << std::endl;
/* 	status.generalized_fold_reduced_nodes += (oldn - status.remaining_nodes); */
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
			} //TODO EG check for restore also need else case in node_vec?
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
