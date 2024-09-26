/******************************************************************************
* reductions.cpp
*
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


#include "struction_reductions.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "flow_graph.h"
#include "push_relabel.h"
#include "key_functions.h"

#include <utility>

using namespace struction;

typedef branch_and_reduce_algorithm::IS_status IS_status;

template<typename F>
void struction_general_reduction::for_each_changed_vertex(branch_and_reduce_algorithm* br_alg, F f) {
    auto& status = br_alg->status;
    for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
        NodeID v = marker.current_vertex(v_idx);

        if (v < status.n && status.node_status[v] == IS_status::not_set) {
            f(v);
        }
    }
}
bool struction_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for_each_changed_vertex(br_alg, [&br_alg, &status](NodeID v) {
        NodeWeight neighbor_weights = 0;
        for (NodeID u : status.graph[v]) {
            neighbor_weights += status.weights[u];
        }

        if (status.weights[v] >= neighbor_weights) {
            br_alg->set(v, IS_status::included);
        }
    });

	//std::cout << "neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool struction_clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&br_alg, &status](NodeID v) {
        auto& neighbors = br_alg->buffers[0];
        auto& neighborhood = br_alg->set_1;
        NodeWeight neighbor_weights = 0;
        neighbors.clear();
        neighborhood.clear();

        for (NodeID neighbor : status.graph[v]) {
            neighbor_weights += status.weights[neighbor];
            neighbors.push_back(neighbor);
            neighborhood.add(neighbor);
        }

        if (status.weights[v] >= neighbor_weights) {
            br_alg->set(v, IS_status::included);
            return;
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

	});

	//std::cout << "clique neighbor fast redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool struction_clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	this->br_alg = br_alg;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;


    for_each_changed_vertex(br_alg, [&](NodeID v) {
        if (partition_into_cliques(v)) {
			br_alg->set(v, IS_status::included);
		}
	});

	//std::cout << "clique neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool struction_clique_neighborhood_reduction::partition_into_cliques(NodeID v) {
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

bool struction_clique_neighborhood_reduction::expand_clique(NodeID max_neighbor, sized_vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set) {
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

bool struction_critical_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->blowing_up)
        return false;

    auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    auto &mapping = br_alg->buffers[0];
    auto &inverse_mapping = br_alg->buffers[1];

    NodeID n = 0;
    for (NodeID node = 0; node < br_alg->status.n; node++) {
        if (status.node_status[node] == IS_status::not_set) {
            mapping[node] = n;
            inverse_mapping[n] = node;
            ++n;
        }
    }

    // build bipartite flow graph
    // node '+ n' shows that we refer to the node in the rest[1] partition
    flow_graph fg;
    fg.start_construction(2 * n + 2);

    const NodeID source = 2 * n;
    const NodeID sink = source + 1;

    for (NodeID id = 0; id < n; id++) {
        NodeID node = inverse_mapping[id];
        // add source and target edges
        fg.new_edge(source, id, status.weights[node]);
        fg.new_edge(id + n, sink, status.weights[node]);

        // add edges between node and its neighbors in the other partition
        for (NodeID neighbor : status.graph[node]) {
            // each outgoing edge has enough capacity to support full flow of the single incoming edge into 'node'
            fg.new_edge(id, mapping[neighbor] + n, status.weights[node]);
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
        NodeID id = fg.getEdgeTarget(source, edge);
        if (fg.getEdgeFlow(source, edge) < fg.getEdgeCapacity(source, edge)) {
            max_cs_set.add(inverse_mapping[id]);
        }
    } endfor

    // isolated nodes in the maximum critical set form the maximum independent critical set
    for (NodeID id = 0; id < n; id++) {
        NodeID node = inverse_mapping[id];
        if (max_cs_set.get(node)) {
            bool isolated = true;

            for (NodeID neighbor : status.graph[node]) {
                if (max_cs_set.get(neighbor)) {
                    isolated = false;
                    break;
                }
            }

            if (isolated) {
                // found isolated node
                br_alg->set(node, IS_status::included);
            }
        }
    }

    //std::cout << "cs redu -> " << (oldn - status.remaining_nodes) << std::endl;

    return oldn != status.remaining_nodes;
}

bool struction_fold2_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	std::array<NodeID, 2> neighbors;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
		if (br_alg->deg(v) == 2) {
            size_t neighbor_count = 0;

			for (NodeID u : status.graph[v]) {
				if (status.weights[u] > status.weights[v]) {
					return;
				}

				neighbors[neighbor_count++] = u;
			}

			if (status.weights[v] >= status.weights[neighbors[0]] + status.weights[neighbors[1]]) {
				// v is always best choice of the three vertices
				br_alg->set(v, IS_status::included);
				return;
			}

			for (NodeID u : status.graph[neighbors[0]]) if (u == neighbors[1]) {
				// clique of size 3
				if (status.weights[v] >= status.weights[u] && status.weights[v] >= status.weights[neighbors[0]]) {
					br_alg->set(v, IS_status::included);
				}
				return;
			}

			fold(br_alg, { v,{ neighbors[0], neighbors[1] } });
		}
	});

	//std::cout << "fold2 redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void struction_fold2_reduction::fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
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
				status.graph.add_edge_directed(neighbor, nodes.main);
				restore_vec.back().node_vecs[i].push_back(neighbor);
			}
		}
	}

	status.graph[nodes.main] = struction_dynamic_graph::neighbor_list(std::move(new_neighbors));
	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.main);
	br_alg->add_next_level_neighborhood(nodes.main);
}

void struction_fold2_reduction::restore(branch_and_reduce_algorithm* br_alg) {
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

void struction_fold2_reduction::apply(branch_and_reduce_algorithm* br_alg) {
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

bool struction_clique_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& set_1 = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	auto& isolated = br_alg->buffers[1];
	std::vector<NodeID> non_isolated;

	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID node) {
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

            if (count != neighbors.size()) return;
        }

        // one of "isolated" members has highest weight of clique: Add to IS
        // also handles completely isolated cliques
        if (max_isolated.weight >= max_non_isolated.weight) {
            br_alg->set(max_isolated.node, IS_status::included);
            return;
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
    });


	//std::cout << "clique redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void struction_clique_reduction::fold(branch_and_reduce_algorithm* br_alg, const weighted_node& isolated, std::vector<NodeID>&& non_isolated) {
	auto& status = br_alg->status;

	br_alg->set(isolated.node, IS_status::folded);
	status.reduction_offset += isolated.weight;

	for (auto node : non_isolated) {
		status.weights[node] -= isolated.weight;
		br_alg->add_next_level_neighborhood(node);
	}

	status.folded_stack.push_back(get_reduction_type());
	br_alg->add_next_level_neighborhood(non_isolated);

	restore_vec.emplace_back(isolated, std::move(non_isolated));
}

void struction_clique_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.isolated.node);
	status.reduction_offset -= data.isolated.weight;

	for (auto node : data.non_isolated) {
		status.weights[node] += data.isolated.weight;
	}

	restore_vec.pop_back();
}

void struction_clique_reduction::apply(branch_and_reduce_algorithm* br_alg) {
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

bool struction_twin_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& twin_candidates_set = br_alg->set_1;
	auto& tmp_set = br_alg->set_2;
	size_t oldn = status.remaining_nodes;

	NodeID twin;
	NodeWeight neighbors_weight;


    for_each_changed_vertex(br_alg, [&](NodeID v) {
        neighbors.clear();
        neighbors_weight = 0;

        for (NodeID neighbor : status.graph[v]) {
            neighbors.push_back(neighbor);
            neighbors_weight += status.weights[neighbor];
        }

        if (status.weights[v] >= neighbors_weight) {
            br_alg->set(v, IS_status::included);
            return;
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
            return;

        if (status.weights[v] + status.weights[twin] >= neighbors_weight) {
            br_alg->set(v, IS_status::included);
            br_alg->set(twin, IS_status::included);
        } else {
            fold(br_alg, v, twin);
        }
    });


	//std::cout << "twin redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void struction_twin_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin) {
	auto& status = br_alg->status;

	restore_vec.push_back({ main, twin });

	br_alg->set(twin, IS_status::folded, true);
	status.weights[main] += status.weights[twin];

	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(main);
	br_alg->add_next_level_neighborhood(main);
}

void struction_twin_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.twin);
	status.weights[data.main] -= status.weights[data.twin];

	restore_vec.pop_back();
}

void struction_twin_reduction::apply(branch_and_reduce_algorithm* br_alg) {
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

bool struction_domination_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;


    for_each_changed_vertex(br_alg, [&](NodeID v) {
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
            return;
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
    });


	//std::cout << "domination redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool struction_generalized_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;
	auto config = br_alg->config;

	cout_handler::disable_cout();


    for_each_changed_vertex(br_alg, [&](NodeID v) {
        NodeWeight neighbors_weight = 0;
        NodeWeight max_neighbor_weight = 0;

        for (NodeID neighbor : status.graph[v]) {
            neighbors_weight += status.weights[neighbor];

            if (status.weights[neighbor] > max_neighbor_weight)
                max_neighbor_weight = status.weights[neighbor];
        }

        if (status.weights[v] >= neighbors_weight) {
            br_alg->set(v, IS_status::included);
            return;
        }
        if (status.weights[v] < max_neighbor_weight) {
            //MWIS in N(v) >= max_neighbor_weight > w(v)
            return;
        }

        if (status.graph[v].size() > branch_and_reduce_algorithm::REDU_RECURSION_LIMIT) {
            return;
        }

        //std::cerr << "%gnr try neighborhood of size " << status.graph[v].size() << " start node: " << v << std::endl;

        // compute MWIS in N(v)
        config.time_limit = status.graph[v].size() / 10.0;
        br_alg->build_induced_neighborhood_subgraph(neighborhood_graph, v);
        branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

        if (!neighborhood_br_alg.run_branch_reduce()) {
            std::cerr << "%generalized_neighborhood_reduction br_call time out" << std::endl;
            return;
        }

        if (status.weights[v] >= neighborhood_br_alg.get_current_is_weight())
            br_alg->set(v, IS_status::included);
    });

	cout_handler::enable_cout();

	//std::cout << "generalized_neighborhood_reduction -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool struction_generalized_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighbors_set = br_alg->set_1;
	auto& MWIS_set = br_alg->set_2;
	auto& reverse_mapping = br_alg->buffers[1];
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;
	auto config = br_alg->config;

	cout_handler::disable_cout();


    for_each_changed_vertex(br_alg, [&](NodeID v) {
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
            return;
        }

        if (status.graph[v].size() > branch_and_reduce_algorithm::REDU_RECURSION_LIMIT) {
            return;
        }

        //std::cerr << "%gfr try neighborhood of size " << status.graph[v].size() << " start node: " << v << std::endl;

        // compute MWIS in N(v)
        config.time_limit = status.graph[v].size() / 10.0;

        br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);
        branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

        if (!neighborhood_br_alg.run_branch_reduce()) {
            std::cerr << "%generalized_fold_reduction br_call time out" << std::endl;
            return;
        }

        NodeWeight MWIS_weight = neighborhood_br_alg.get_current_is_weight();
        NodeWeight min_MWIS_neighbor_weight = std::numeric_limits<NodeWeight>::max();

        if (status.weights[v] >= MWIS_weight) {
            // same as in generalized_neighborhood_reduction
            br_alg->set(v, IS_status::included);
            return;
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
            return;
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
            return;
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
    });

	cout_handler::enable_cout();

	//if (oldn - status.remaining_nodes != 0)
		//std::cout << "%generalized_fold_reduction improvement: " << (oldn - status.remaining_nodes) << std::endl;

	//std::cout << "generalized_fold_reduction -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void struction_generalized_fold_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main_node, fast_set& MWIS_set, NodeWeight MWIS_weight) {
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
				status.graph.add_edge_directed(neighbor, nodes.main);
				node_vec.push_back(neighbor);
			}
		}

		data.MWIS_node_vecs.push_back(std::move(node_vec));
	}

	status.graph[nodes.main] = struction_dynamic_graph::neighbor_list(std::move(new_neighbors));
	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.main);
	br_alg->add_next_level_neighborhood(nodes.main);
}

void struction_generalized_fold_reduction::restore(branch_and_reduce_algorithm* br_alg) {
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

void struction_generalized_fold_reduction::apply(branch_and_reduce_algorithm* br_alg) {
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

template<typename struction_type, reduction_type type, int vertex_increase>
bool iterative_struction<struction_type, type, vertex_increase>::reduce(branch_and_reduce_algorithm* br_alg) {
    auto &status = br_alg->status;
    bool applied = false;
    for_each_changed_vertex(br_alg, [&](NodeID v) {
        if (br_alg->deg(v) > br_alg->config.struction_degree || !s.reduce(br_alg, v, s.removed_vertices(br_alg, v) + vertex_increase)) return;

        status.folded_stack.push_back(get_reduction_type());
        applied = true;
    });

    return applied;
}

template class iterative_struction<extended_struction<false>, reduction_type::struction_decrease, -1>;
template class iterative_struction<extended_struction<false>, reduction_type::struction_plateau, 0>;
template class iterative_struction<extended_struction<true>, reduction_type::struction_decrease, -1>;
template class iterative_struction<extended_struction<true>, reduction_type::struction_plateau, 0>;
template class iterative_struction<original_struction<false>, reduction_type::struction_decrease, -1>;
template class iterative_struction<original_struction<false>, reduction_type::struction_plateau, 0>;
template class iterative_struction<original_struction<true>, reduction_type::struction_decrease, -1>;
template class iterative_struction<original_struction<true>, reduction_type::struction_plateau, 0>;

template<typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::reduce(branch_and_reduce_algorithm *br_alg) {
    this->br_alg = br_alg;
    auto &status = br_alg->status;

    init_blow_up_phase();

    while (!is_done() && clean_up_queue()) {
        NodeID n = queue.maxElement();
        Gain key = denoise(queue.maxValue());
        size_t set_limit_by_key = f.set_limit(br_alg, n, key);
        size_t plain_set_limit = br_alg->config.set_limit;
        if (f.set_estimate(br_alg, n, key) > plain_set_limit)
            break;

        queue.deleteMax();
        if (s.reduce(br_alg, n, std::min(set_limit_by_key, plain_set_limit))) {
            //struction was successfully executed
            if (br_alg->config.backtrack_style != ::mmwis::MISConfig::Backtrack_Type::IMMEDIATE_EXCLUDE && update_set.add(n))
                update_list.emplace_back(n, key);
            status.folded_stack.push_back(get_reduction_type());
            ++blow_ups;
            //"blow up" actually can reduce remaining nodes if multiple blow ups per phase are enabled.
            br_alg->min_kernel = std::min(status.remaining_nodes, br_alg->min_kernel);
        } else if (plain_set_limit <= set_limit_by_key) {
            break;
        } else {
            //reinsert with tighter bound
            update_queue_by_key(n, f.key_by_set_estimate(br_alg, n, set_limit_by_key + 1));
        }
    }

    return blow_ups != 0;
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::update_queue(NodeID n) {
    if (br_alg->deg(n) > br_alg->config.struction_degree) return;
    update_queue_by_key(n, f.key(br_alg, n));
}


template<typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::clean_up_queue() {
    auto &status = br_alg->status;
    update_set.resize(status.n);

    //Update candidate set
    for (NodeID n : marker.next)
        update_queue(n);
    marker.clear_next();

    while (queue.size()) {
        NodeID n = queue.maxElement();
        if (status.node_status[n] == IS_status::not_set && br_alg->deg(n) <= br_alg->config.struction_degree) return true;
        if (update_set.add(n))
            update_list.emplace_back(n, denoise(queue.maxValue()));
        queue.deleteMax();
    }
    return false;
}

template<typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::is_done() {
    auto &config = br_alg->config;
    auto &status = br_alg->status;
    size_t cur_kernel = status.remaining_nodes;
    size_t max_kernel = config.phase_blow_up_factor * phase_start_kernel;
    return max_kernel < cur_kernel || blow_ups == config.phase_blow_ups;
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::restore(branch_and_reduce_algorithm *br_alg) {
    s.restore(br_alg);
    restored = true;
}


template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::reset(branch_and_reduce_algorithm* br_alg, size_t comp_size) {
    restored = false;
    queue.clear();
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::init_blow_up_phase() {
    auto &status = br_alg->status;

    update_set.resize(status.n);
    blow_ups = 0;
    if (restored) {
        for (NodeID n : added_list)
            if (queue.contains(n))
                queue.deleteNode(n);
        for (auto &e : update_list)
            update_queue_by_key(e.first, e.second);
    } else {
        for_each_changed_vertex(br_alg, [&](NodeID n) {
            update_queue(n);
        });
    }
    phase_start_kernel = status.remaining_nodes;
    update_list.clear();
    added_list.clear();
    update_set.clear();

    restored = false;
}

template class blow_up_struction<DegreeKey>;
template class blow_up_struction<IncreaseKey>;
template class blow_up_struction<ApproximateIncreaseKey>;
template class blow_up_struction<RandomKey>;


bool path_reduction::reduce(branch_and_reduce_algorithm *br_alg) {
    this->br_alg = br_alg;
    auto &status = br_alg->status;
    size_t old_n = status.remaining_nodes;

    br_alg->buffers[0].clear();
    br_alg->buffers[1].clear();

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        enqueue_node<false>(v);
    });

    while (reduce_degree_one_node() || reduce_path());
    return old_n != status.remaining_nodes;
}


void path_reduction::restore(branch_and_reduce_algorithm* br_alg) {
    auto &status = br_alg->status;
    auto &restore_data = restore_vec.back();
    auto &path = restore_data.path;
    NodeID v = path[restore_data.start], w = path[restore_data.end - 1];
    for (auto &n_w : restore_data.node_weights)
        status.weights[n_w.first] = n_w.second;
    for (size_t i = restore_data.start; i < restore_data.end; ++i) {
        status.node_status[path[i]] = IS_status::not_set;
        status.remaining_nodes++;
    }
    if (restore_data.relink) {
        status.graph.relink_directed(v, w, path[restore_data.start + 1]);
        status.graph.relink_directed(w, v, path[restore_data.end - 2]);
    } else {
        status.graph.restore_node(w);
        status.graph.restore_node(v);
    }
    status.reduction_offset -= restore_data.offset;
    restore_vec.pop_back();
}

void path_reduction::apply(branch_and_reduce_algorithm* br_alg) {
    auto &restore_data = restore_vec.back();
    auto path = restore_data.path;
    restore(br_alg);
    apply(path);
}

bool path_reduction::dequeue_node(sized_vector<NodeID> &queue, NodeID &n, size_t degree) {
    for (;;) {
        if (queue.empty()) return false;
        n = queue.back();
        queue.pop_back();
        if (br_alg->status.node_status[n] == IS_status::not_set && br_alg->deg(n) == degree) return true;
    }
}

bool path_reduction::reduce_degree_one_node() {
    NodeID n;
    if (!dequeue_node(br_alg->buffers[0], n, 1)) return false;

    auto &status = br_alg->status;

    NodeID u = status.graph[n][0];
    if (status.weights[n] >= status.weights[u]) {
        status.node_status[n] = IS_status::included;
        status.node_status[u] = IS_status::excluded;
        status.graph.hide_node(n);
        status.graph.hide_node(u);
        status.remaining_nodes -= 2;
        status.is_weight += status.weights[n];
        status.modified_stack.push_back(n);
        status.modified_stack.push_back(u);
        for (NodeID v : status.graph[u])
            enqueue_node(v);
    } else {
        restore_vec.emplace_back();
        auto &restore_data = restore_vec.back();

        br_alg->set(n, IS_status::folded);
        add_reduction_offset(status.weights[n], restore_data);
        reassign_weight(u, status.weights[u] - status.weights[n], restore_data);
        enqueue_node(u);
    }
    return true;

    /*sized_vector<NodeID> &path = br_alg->buffers[2];
    find_max_deg_1_path(n, path);
    NodeWeight w_e, w_i;
    find_MIS_on_deg_1_path(w_i, w_e, path);
    bool keep_v = br_alg->deg(path.back()) != 1 && w_i > w_e;
    for (size_t i = 0; i < path.size() - keep_v; ++i) {
        status.node_status[path[i]] = IS_status::folded;
        status.remaining_nodes--;
    }
    status.graph.hide_node(path[path.size() - 1 - keep_v]);
    if (keep_v) {
        status.weights[path.back()] = w_i - w_e;
        enqueue_node(path.back());
        status.reduction_offset += w_e;
    } else {
        status.reduction_offset += std::max(w_i, w_e);
        for (NodeID u : status.graph[path.back()])
            enqueue_node(u);
    }
    return true;*/
}

bool path_reduction::reduce_path() {
    NodeID n;
    if (!dequeue_node(br_alg->buffers[1], n, 2)) return false;

    auto &status = br_alg->status;
    auto &path = br_alg->buffers[2];
    find_max_path(n, path);
    find_MIS_on_path(path);

    NodeID v = path.front(), w = path.back();
    restore_vec.emplace_back();
    auto &restore_data = restore_vec.back();
    if (v == w) {
        NodeWeight w_e = w_e_e;
        NodeWeight w_i = w_i_i - status.weights[v];
        if (br_alg->deg(v) != 2 && w_i > w_e) {
            //fold
            add_reduction_offset(w_e, restore_data);
            fold_path(1, path.size() - 1, restore_data);
            reassign_weight(v, w_i - w_e, restore_data);
            enqueue_node(v);
        } else {
            //include/exclude v
            fold_path(0, path.size() - 1, restore_data);
            status.node_status[v] = w_i > w_e ? IS_status::folded : IS_status::excluded; //set node state to folded instead of included to be consistent with mis_weight
            add_reduction_offset(std::max(w_i, w_e), restore_data);
        }
    } else {
        bool connected = are_connected(v, w);
        if (connected || w_i_i <= w_i_e || w_i_i <= w_e_i) {
            bool keep_v = w_i_e > w_e_e;
            bool keep_w = w_e_i > w_e_e;

            add_reduction_offset(w_e_e, restore_data);
            fold_path(keep_v, path.size() - keep_w,restore_data);
            reassign_weight(v, w_i_e - w_e_e, restore_data);
            reassign_weight(w,  w_e_i - w_e_e, restore_data);
            if (!connected && keep_v && keep_w)
                reconnect(v, w, restore_data);
        } else {
            int gap = static_cast<int>(w_e_e - w_e_i - w_i_e + w_i_i);
            NodeID c = path[1];
            if (gap >= 0) {
                if (path.size() <= 3) return true;
                add_reduction_offset(w_e_e - gap,restore_data);
                fold_path(2, path.size() - 1,restore_data);
                reassign_weight(v, w_i_i - w_e_i, restore_data);
                reassign_weight(c, gap, restore_data);
                reassign_weight(w, w_i_i - w_i_e, restore_data);
                reconnect(c, w, restore_data);
            } else {
                if (path.size() <= 4) return true;
                NodeID d = path[path.size() - 2];
                add_reduction_offset(w_e_e + gap, restore_data);
                fold_path(2, path.size() - 2,restore_data);
                reassign_weight(v, w_i_e - w_e_e, restore_data);
                reassign_weight(c, -gap, restore_data);
                reassign_weight(d, -gap, restore_data);
                reassign_weight(w, w_e_i - w_e_e, restore_data);
                reconnect(c, d, restore_data);
            }
        }
        enqueue_node(v);
        enqueue_node(w);
    }

    restore_vec.back().path.resize(path.size());
    for (NodeID n : path)
        restore_vec.back().path.push_back(n);
    return true;
}

void path_reduction::reassign_weight(NodeID n, NodeWeight w, restore_data &data) {
    br_alg->status.weights[n] = w;
    data.node_weights.emplace_back(n, w);
}
void path_reduction::reconnect(NodeID v, NodeID w, restore_data &data) {
    auto &status = br_alg->status;
    status.graph.add_edge_undirected(v, w);
    data.relink = true;
}

void path_reduction::fold_path(size_t start, size_t end, restore_data &data) {
    auto &path = br_alg->buffers[2];
    auto &status = br_alg->status;
    for (size_t i = start; i < end; ++i) {
        status.node_status[path[i]] = IS_status::folded;
        status.remaining_nodes--;
    }
    status.graph.hide_node(path[start]);
    status.graph.hide_node(path[end - 1]);
    br_alg->add_next_level_neighborhood(path[start]);
    br_alg->add_next_level_neighborhood(path[end - 1]);
    status.modified_stack.push_back(path[end - 1]);
    data.start = start;
    data.end = end;
}

void path_reduction::add_reduction_offset(size_t offset, restore_data &data) {
    br_alg->status.reduction_offset += offset;
    data.offset = offset;
}
bool path_reduction::are_connected(NodeID v, NodeID w) const {
    auto &neighbors = br_alg->status.graph[v];
    return std::find(neighbors.begin(), neighbors.end(), w) != neighbors.end();
}

void path_reduction::find_MIS_on_path(sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    NodeWeight w_i = status.weights[path[1]], w_e = 0;
    find_MIS_on_path(w_i, w_e, path);
    w_e_i = w_i;
    w_e_e = w_e;
    w_e = status.weights[path[0]];
    w_i = 0;
    find_MIS_on_path(w_i, w_e, path);
    w_i_i = w_i;
    w_i_e = w_e;
}

template<bool track_choices>
void path_reduction::find_MIS_on_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    for (size_t i = 2; i < path.size(); ++i) {
        if (track_choices)
            br_alg->bool_buffer[i - 1] = w_i > w_e;
        NodeWeight next_e = std::max(w_e, w_i);
        w_i = w_e + status.weights[path[i]];
        w_e = next_e;
    }
}
template<bool track_choices>
void path_reduction::find_MIS_on_deg_1_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    w_i = status.weights[path[0]];
    w_e = 0;
    NodeWeight next_w_e;
    for (size_t i = 1; i < path.size(); ++i) {
        if (track_choices)
            br_alg->bool_buffer[i - 1] = w_i > w_e;
        next_w_e = std::max(w_e, w_i);
        w_i = w_e + status.weights[path[i]];
        w_e = next_w_e;
    }
}
void path_reduction::find_max_path(NodeID n, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    path.clear();
    path.push_back(n);

    NodeID last;
    NodeID next;
    for (NodeID current : status.graph[n]) {
        last = n;
        std::reverse(path.begin(), path.end());
        path.push_back(current);
        while (next_node_on_path(current, last, n, next)) {
            path.push_back(next);
            last = current;
            current = next;
        }
        if (path.front() == path.back())
            break;
    }
}
void path_reduction::find_max_deg_1_path(NodeID n, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;

    NodeID last = n;
    NodeID current = status.graph[n][0];
    NodeID next;

    path.clear();
    path.push_back(last);
    path.push_back(current);
    while (next_node_on_path(current, last, n, next)) {
        path.push_back(next);
        last = current;
        current = next;
    }
}
bool path_reduction::next_node_on_path(NodeID current, NodeID last, NodeID first, NodeID &next) {
    auto &status = br_alg->status;
    if (br_alg->deg(current) != 2 || current == first)
        return false;

    auto &neighbors = status.graph[current];
    next = neighbors[0] != last ? neighbors[0] : neighbors[1];
    return true;
}

template<bool add_global>
void path_reduction::enqueue_node(NodeID n) {
    if (br_alg->status.node_status[n] != IS_status::not_set)
        return;
    if (br_alg->deg(n) == 1)
        br_alg->buffers[0].push_back(n);
    /*else if (br_alg->deg(n) == 2)
        br_alg->buffers[1].push_back(n);*/
    else if (add_global)
        br_alg->add_next_level_node(n);
}
void path_reduction::apply(sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    NodeID v = path.front(), w = path.back();
    NodeWeight w_e = status.node_status[v] != IS_status::excluded ? status.weights[v] : 0;
    NodeWeight w_i = status.node_status[v] != IS_status::excluded ? 0 : status.weights[path[1]];
    find_MIS_on_path<true>(w_i, w_e, path);

    auto &choices = br_alg->bool_buffer;
    size_t i = path.size() - 2;
    bool include_v_i = status.node_status[w] == IS_status::included ? 0 : choices[i];
    for (; i >= 1; --i) {
        NodeID n = path[i];
        if (include_v_i) {
            status.node_status[n] = IS_status::included;
            status.is_weight += status.weights[n];
            include_v_i = false;
        } else {
            status.node_status[n] = IS_status::excluded;
            include_v_i = choices[i - 1];
        }
    }
}


template<reduction_type type, int new_nodes>
reduction_ptr struction::make_iterative_struction(const ::mmwis::MISConfig &config, size_t n) {
    const auto s = config.struction_type;
    if (s == Struction_Type::ORIGINAL) return reduction_ptr(new iterative_struction<original_struction<false>,type,new_nodes>(n));
    if (s == Struction_Type::MODIFIED) return reduction_ptr(new iterative_struction<original_struction<true>,type,new_nodes>(n));
    if (s == Struction_Type::EXTENDED_REDUCED) return reduction_ptr(new iterative_struction<extended_struction<true>,type,new_nodes>(n));
    return reduction_ptr(new iterative_struction<extended_struction<false>,type,new_nodes>(n));
};

reduction_ptr struction::make_decreasing_struction(const ::mmwis::MISConfig &config, size_t n) {
    return struction::make_iterative_struction<reduction_type::struction_decrease, -1>(config, n);
};

reduction_ptr struction::make_plateau_struction(const ::mmwis::MISConfig &config, size_t n) {
    return struction::make_iterative_struction<reduction_type::struction_plateau, 0>(config, n);
};

reduction_ptr struction::make_increasing_struction(const ::mmwis::MISConfig &config, size_t n) {
    if (config.key_type == Key_Type::RANDOM) return reduction_ptr(new blow_up_struction<RandomKey>(config, n));
    if (config.key_type == Key_Type::DEGREE) return reduction_ptr(new blow_up_struction<DegreeKey>(config, n));
    if (config.key_type == Key_Type::INCREASE) return reduction_ptr(new blow_up_struction<IncreaseKey>(config, n));
    return reduction_ptr(new blow_up_struction<ApproximateIncreaseKey>(config, n));
};
