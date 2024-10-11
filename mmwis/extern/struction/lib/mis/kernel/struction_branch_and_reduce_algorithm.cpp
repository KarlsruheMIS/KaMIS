/******************************************************************************
* Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
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

#include "struction_branch_and_reduce_algorithm.h"
#include "struction_reductions.h"
#include "ils/ils.h"
#include "hils/hils.h"
#include "ils/local_search.h"
#include "strongly_connected_components.h"
#include "graph_extractor.h"
#include "ReduceAndPeel.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <numeric>

using namespace struction;

constexpr NodeID branch_and_reduce_algorithm::BRANCHING_TOKEN;
std::streambuf* cout_handler::cout_rdbuf_backup = std::cout.rdbuf();
std::stringstream cout_handler::buffered_output;
int cout_handler::disable_count = 0;

branch_and_reduce_algorithm::~branch_and_reduce_algorithm() {
}

branch_and_reduce_algorithm::branch_and_reduce_algorithm(graph_access& G, const ::mmwis::MISConfig& config, bool called_from_fold)
        : config(config), global_status(G), set_1(global_status.n), set_2(global_status.n), double_set(global_status.n * 2),
        buffers(4, sized_vector<NodeID>(global_status.n)), bool_buffer(global_status.n), zero_vec(global_status.n, 0) {
    if (called_from_fold) {
        global_status.transformations = make_reduction_vector<
                struction_neighborhood_reduction
                , struction_fold2_reduction
                , struction_clique_reduction
                , struction_domination_reduction
                , struction_twin_reduction>(global_status.n);
        global_status.num_reductions = global_status.transformations.size();
    } else if (config.reduction_style == ::mmwis::MISConfig::StructionReduction_Style::DENSE) {
        global_status.transformations = make_reduction_vector<
                struction_neighborhood_reduction
                , struction_fold2_reduction
                , struction_clique_reduction
                , struction_domination_reduction
                , struction_twin_reduction
                , struction_generalized_fold_reduction>(global_status.n);
        global_status.num_reductions = global_status.transformations.size();
    } else {
        // ::mmwis::MISConfig::Reduction_Style::NORMAL
        if (!config.plain_struction) {
            global_status.transformations.emplace_back(new struction_neighborhood_reduction(global_status.n));
            global_status.transformations.emplace_back(new struction_fold2_reduction(global_status.n));
            if (!config.disable_clique)
                global_status.transformations.emplace_back(new struction_clique_reduction(global_status.n));
            global_status.transformations.emplace_back(new struction_domination_reduction(global_status.n));
            global_status.transformations.emplace_back(new struction_twin_reduction(global_status.n));
            if (!config.disable_clique_neighborhood)
                global_status.transformations.emplace_back(new struction_clique_neighborhood_reduction(global_status.n));
            if (!config.disable_generalized_fold)
                global_status.transformations.emplace_back(new struction_generalized_fold_reduction(global_status.n));
        }
        if (config.struction_type != ::mmwis::MISConfig::Struction_Type::NONE) {
            global_status.transformations.push_back(make_decreasing_struction(config, global_status.n));
            global_status.transformations.push_back(make_plateau_struction(config, global_status.n));
        }
        if (!config.plain_struction && !config.disable_critical_set) {
            global_status.transformations.emplace_back(new struction_critical_set_reduction(global_status.n));
        }


        global_status.num_reductions = global_status.transformations.size();
        if (!config.disable_blow_up) {
            global_status.transformations.push_back(make_increasing_struction(config, global_status.n));
        }
    }

    global_transformation_map.resize(STRUCTION_REDUCTION_NUM);
    for (size_t i = 0; i < global_status.transformations.size(); i++) {
        global_transformation_map[global_status.transformations[i]->get_reduction_type()] = i;
    }

    set_local_reductions = [this, called_from_fold, &config]() {
        if (this->config.reduction_style == ::mmwis::MISConfig::StructionReduction_Style::DENSE) {
            status.transformations = make_reduction_vector<
                    struction_neighborhood_reduction
                    , struction_fold2_reduction
                    , struction_clique_reduction
                    , struction_domination_reduction
                    , struction_twin_reduction
                    , struction_clique_neighborhood_reduction_fast>(status.n);
        } else {
            // ::mmwis::MISConfig::StructionReduction_Style::NORMAL
            if (called_from_fold) {
                status.transformations = make_reduction_vector<
                        struction_neighborhood_reduction
                        , struction_fold2_reduction
                        , struction_clique_reduction
                        , struction_domination_reduction
                        , struction_twin_reduction
                        , struction_clique_neighborhood_reduction
                        , struction_critical_set_reduction>(status.n);
            } else {
                status.transformations.clear();
                status.transformations.emplace_back(new struction_neighborhood_reduction(global_status.n));
                status.transformations.emplace_back(new struction_fold2_reduction(global_status.n));
                status.transformations.emplace_back(new struction_clique_reduction(global_status.n));
                status.transformations.emplace_back(new struction_domination_reduction(global_status.n));
                status.transformations.emplace_back(new struction_twin_reduction(global_status.n));
                status.transformations.push_back(make_decreasing_struction(config, global_status.n));
                status.transformations.emplace_back(new struction_critical_set_reduction(global_status.n));
            }
        }
        status.num_reductions = status.transformations.size();
        local_transformation_map.resize(STRUCTION_REDUCTION_NUM);

        for (size_t i = 0; i < status.transformations.size(); i++) {
            local_transformation_map[status.transformations[i]->get_reduction_type()] = i;
        }
    };
}

void branch_and_reduce_algorithm::push_nodes(size_t nodes) {
    resize(status.n + nodes);
    status.push_nodes(nodes);
}
void branch_and_reduce_algorithm::pop_nodes(size_t nodes) {
    resize(status.n - nodes);
    status.pop_nodes(nodes);
}
void branch_and_reduce_algorithm::resize(size_t size) {
    set_1.resize(size);
    set_2.resize(size);
    double_set.resize(2 * size);
    bool_buffer.resize(size + 1);
    zero_vec.resize(size, 0);
    for (auto &buffer : buffers) {
        buffer.resize(size);
    }
    for (auto& transformation : status.transformations) {
        transformation->marker.resize(size);
    }
}

size_t branch_and_reduce_algorithm::deg(NodeID node) const {
	return status.graph[node].size();
}

template<bool update_dependency_checking>
void branch_and_reduce_algorithm::set(NodeID node, IS_status mis_status, bool push_modified) {
	status.node_status[node] = mis_status;
	status.remaining_nodes--;
	status.graph.hide_node(node);

	if (push_modified)
		status.modified_stack.push_back(node);

	if (mis_status == IS_status::included) {
		status.is_weight += status.weights[node];

		for (auto neighbor : status.graph[node]) {
			status.node_status[neighbor] = IS_status::excluded;
			status.remaining_nodes--;
			status.graph.hide_node(neighbor);

			status.modified_stack.push_back(neighbor);
			//if (update_dependency_checking)
			    add_next_level_neighborhood(neighbor);
		}
	}
	else /*if (update_dependency_checking)*/ {
		add_next_level_neighborhood(node);
	}
}
template void branch_and_reduce_algorithm::set<true>(NodeID node, IS_status mis_status, bool push_modified);
template void branch_and_reduce_algorithm::set<false>(NodeID node, IS_status mis_status, bool push_modified);

void branch_and_reduce_algorithm::unset(NodeID node, bool restore) {
	if (status.node_status[node] == IS_status::included) {
		status.is_weight -= status.weights[node];
	}

	status.node_status[node] = IS_status::not_set;
	status.remaining_nodes++;

	if (restore)
		status.graph.restore_node(node);
}


void branch_and_reduce_algorithm::fill_global_greedy() {
	auto& nodes = buffers[0];
	nodes.clear();

	for (size_t node = 0; node < global_status.n; node++) {
		if (global_status.node_status[node] == IS_status::not_set)
			nodes.push_back(node);
	}

	// sort in descending order by node weights
	std::sort(nodes.begin(), nodes.end(), [this](NodeID lhs, NodeID rhs) {
		return global_status.weights[lhs] > global_status.weights[rhs];
	});

	for (NodeID node : nodes) {
		bool free_node = true;

		for (NodeID neighbor : global_status.graph[node]) {
			if (global_status.node_status[neighbor] == IS_status::included) {
				free_node = false;
				break;
			}
		}

		if (free_node) {
			global_status.node_status[node] = IS_status::included;
			global_status.is_weight += global_status.weights[node];
		}
		else {
			global_status.node_status[node] = IS_status::excluded;
		}
	}
}

void branch_and_reduce_algorithm::greedy_initial_is(graph_access& G, sized_vector<NodeID>& tmp_buffer) {
	auto nodes = tmp_buffer;
	nodes.set_size(G.number_of_nodes());

	for (size_t i = 0; i < nodes.size(); i++) {
		nodes[i] = i;
	}

	// sort in descending order by node weights
	std::sort(nodes.begin(), nodes.end(), [&](NodeID lhs, NodeID rhs) {
		return G.getNodeWeight(lhs) > G.getNodeWeight(rhs);
	});

	for (NodeID node : nodes) {
		bool free_node = true;

		forall_out_edges(G, edge, node) {
			NodeID neighbor = G.getEdgeTarget(edge);
			if (G.getPartitionIndex(neighbor) == 1) {
				free_node = false;
				break;
			}
		} endfor

		if (free_node) G.setPartitionIndex(node, 1);
	}
}


void branch_and_reduce_algorithm::compute_ils_pruning_bound() {
	is_ils_best_solution = true;

    std::cout << "started " << t.elapsed() << std::endl;
    cout_handler::disable_cout();

	auto config_cpy = config;
	config_cpy.time_limit = config_cpy.time_limit * status.n / total_ils_node_count / 100;

	best_weight = status.reduction_offset + status.is_weight + run_ils(config_cpy, *local_graph, buffers[0], 1000);
    best_is_time = t.elapsed();
	cout_handler::enable_cout();
	std::cout << "finished " << t.elapsed() << std::endl;
	std::cout << (get_current_is_weight() + best_weight) << "," << t.elapsed() << std::endl;
}

NodeWeight branch_and_reduce_algorithm::compute_cover_pruning_bound() {
	// Gather remaining nodes
	auto& nodes = buffers[0];
	nodes.clear();

	for (size_t node = 0; node < status.n; node++) {
		if (status.node_status[node] == IS_status::not_set)
			nodes.emplace_back(node);
	}

	// Sort by descending weight
	// Break ties by degree
	std::sort(nodes.begin(), nodes.end(), [&](NodeID lhs, NodeID rhs) {
		return (status.weights[lhs] > status.weights[rhs] || (status.weights[lhs] == status.weights[rhs] && deg(lhs) > deg(rhs)));
	});

	// Compute node mapping
	NodeID current_node = 0;
	auto& node_mapping = buffers[1];
	node_mapping.set_size(status.n);

	for (NodeID node : nodes)
		node_mapping[node] = current_node++;

	// Init cliques
	auto& clique = buffers[2];
	auto& clique_weight = buffers[3];
	auto& clique_sizes = buffers[4];
	auto& covered = buffers[5];

	clique.set_size(nodes.size());
	clique_weight.set_size(nodes.size());
	clique_sizes.set_size(nodes.size());
	covered.set_size(nodes.size());

	std::fill(clique_weight.begin(), clique_weight.end(), 0);
	std::fill(clique_sizes.begin(), clique_sizes.end(), 0);
	std::fill(covered.begin(), covered.end(), false);

	for (NodeID node : nodes)
		clique[node_mapping[node]] = node_mapping[node];

	auto& neigh_clique_sizes = buffers[6];
	neigh_clique_sizes.set_size(nodes.size());

	for (NodeID node : nodes) {
		NodeID v = node_mapping[node];
		// Find heaviest neighboring clique
		NodeID heaviest_clique = v;
		std::fill(neigh_clique_sizes.begin(), neigh_clique_sizes.end(), 0);

		for (NodeID neighbor : status.graph[node]) {
			if (status.node_status[neighbor] == IS_status::not_set) {
				NodeID w = node_mapping[neighbor];

				if (covered[w]) {
					NodeID c = clique[w];
					neigh_clique_sizes[c]++;

					if (neigh_clique_sizes[c] == clique_sizes[c] && clique_weight[c] > clique_weight[heaviest_clique])
						heaviest_clique = c;
				}
			}
		}

		// Update clique weights/sizes
		clique[v] = heaviest_clique;
		clique_weight[heaviest_clique] = std::max(clique_weight[heaviest_clique], status.weights[node]);
		clique_sizes[heaviest_clique]++;

		// Node is covered
		covered[v] = true;
	}

	// Check clique cover conditions
	// for (NodeID node : nodes) {
	//   NodeID v = node_mapping[node];
	//   NodeID c = clique[v];
	//   NodeID cs = clique_sizes[c];
	//   NodeID num_cn = 1;
	// 	for (NodeID neighbor : status.graph[node]) {
	// 		if (status.node_status[neighbor] == IS_status::not_set) {
	//       NodeID w = node_mapping[neighbor];
	//       if (clique[w] == c) num_cn++;
	//     }
	//   }
	// }

	// Keep unique cliques
	std::sort(clique.begin(), clique.end());
	auto end_iter = std::unique(clique.begin(), clique.end());

	// Add clique weights to get upper bound
	NodeWeight upper_bound = 0;
	for (auto iter = clique.begin(); iter != end_iter; iter++)
		upper_bound += clique_weight[*iter];

	return upper_bound;
}

size_t branch_and_reduce_algorithm::run_ils(const ::mmwis::MISConfig& config, graph_access& G, sized_vector<NodeID>& tmp_buffer, size_t max_swaps) {
    double time_limit = 0.01*(config.time_limit - t.elapsed());
    if (!config.perform_hils) {
        greedy_initial_is(G, tmp_buffer);
        mmwis::ils local_search(config);
        local_search.perform_ils(G, max_swaps, time_limit);
    } else {
        if (config.reduce_and_peel) {
            ReduceAndPeel reducer(G);
            reducer.reduceAndPeel();
        } else
            greedy_initial_is(G, tmp_buffer);
        hils local_search(config);
        local_search.perform_ils(G, max_swaps, time_limit);
    }

	size_t solution_weight = 0;

	forall_nodes(G, node) {
		if (G.getPartitionIndex(node) == 1) {
			solution_weight += G.getNodeWeight(node);
		}
	} endfor

	return solution_weight;
}

void branch_and_reduce_algorithm::init_transformation_step(reduction_ptr &reduction) {
	if (!reduction->has_run) {
        reduction->marker.fill_current_ascending(status.n);
        reduction->marker.clear_next();
        reduction->has_run = true;
	}
	else {
        reduction->marker.get_next();
	}
}

void branch_and_reduce_algorithm::add_next_level_node(NodeID node) {
	// mark node for next round of status.reductions
	for (auto& reduction : status.transformations) {
		if (reduction->has_run) {
			reduction->marker.add(node);
		}
	}
}

void branch_and_reduce_algorithm::add_next_level_neighborhood(NodeID node) {
	// node has been excluded in mis -> neighboring vertices are interesting for next round of reduction
	for (auto neighbor : status.graph[node]) {
		add_next_level_node(neighbor);
	}
}

void branch_and_reduce_algorithm::add_next_level_neighborhood(const std::vector<NodeID>& nodes) {
	for (auto node : nodes) {
		add_next_level_neighborhood(node);
	}
}

void branch_and_reduce_algorithm::initial_reduce() {
	std::swap(global_transformation_map, local_transformation_map);
	status = std::move(global_status);
	reduce_graph_internal();

    min_kernel = status.remaining_nodes;
	if (!config.disable_blow_up && status.transformations.size() != status.num_reductions)
	    cyclic_blow_up();

    status.modified_stack.push_back(BRANCHING_TOKEN);

	global_status = std::move(status);
	std::swap(global_transformation_map, local_transformation_map);
}

void branch_and_reduce_algorithm::reduce_graph() {
	initial_reduce();
	status = std::move(global_status);
}

void branch_and_reduce_algorithm::reduce_graph_internal() {
    size_t active_reduction_index = 0;
    while (active_reduction_index < status.num_reductions && t.elapsed() <= config.time_limit) {
        auto& reduction = status.transformations[active_reduction_index];

        init_transformation_step(reduction);
        bool progress = reduction->reduce(this);
        active_reduction_index = progress ? 0 : active_reduction_index + 1;
    }
    timeout |= t.elapsed() > config.time_limit;

	//std::cout << "\ncurrent weight: " << status.is_weight << "  weight offset: " << status.reduction_offset << "  remaining nodes: " << status.remaining_nodes << std::endl;
}


bool branch_and_reduce_algorithm::blow_up_graph_internal() {
    bool blown_up = false;
    for (size_t active_blow_up_index = status.num_reductions; active_blow_up_index < status.transformations.size(); ++active_blow_up_index) {
        auto &blow_up = status.transformations[active_blow_up_index];
        init_transformation_step(blow_up);
        blown_up |= blow_up->reduce(this);
        if (t.elapsed() > config.time_limit) break;
    }
    timeout |= t.elapsed() > config.time_limit;
    return blown_up;
}

void branch_and_reduce_algorithm::cyclic_blow_up() {
    blowing_up = true;
    size_t phase_count = 0;

#ifdef OUTPUT_GRAPH_CONVERGENCE
    std::cout << t.elapsed() << "," << min_kernel << std::endl;
#endif
    const double alpha = config.global_blow_up_factor;
    const int X = config.max_unimproving_phases;
    while (status.remaining_nodes < alpha * min_kernel  && phase_count < X) {
        if (timeout || !blow_up_graph_internal() || timeout) break;
        //REDUCE
        reduce_graph_internal();

        if (status.remaining_nodes < min_kernel) {
            min_kernel = status.remaining_nodes;
#ifdef OUTPUT_GRAPH_CONVERGENCE
            std::cout << t.elapsed() << "," << min_kernel << std::endl;
#endif
            phase_count = 0;
        } else if (config.backtrack_style == ::mmwis::MISConfig::IMMEDIATE_EXCLUDE ||
                   config.backtrack_style == ::mmwis::MISConfig::IMMEDIATE_TIE_BREAKING) {
            ++phase_count;
            undo_blow_up();
        }
    }
    if (config.backtrack_style != ::mmwis::MISConfig::NO_BACKTRACK)
        undo_blow_up();
#ifdef OUTPUT_GRAPH_CONVERGENCE
    std::cout << t.elapsed() << "," << min_kernel << std::endl;
#endif

    blowing_up = false;
    //std::cout << "\ncurrent weight: " << status.is_weight << "  weight offset: " << status.reduction_offset << "  remaining nodes: " << status.remaining_nodes << std::endl;
}

bool branch_and_reduce_algorithm::branch_reduce_recursive() {
    recursive_mapping.resize(status.n, 0);
    recursive_comp_map.resize(status.n, 0);
	build_graph_access(recursive_graph, recursive_mapping);
	size_t comp_count = strongly_connected_components().strong_components(recursive_graph, recursive_comp_map);

	if (comp_count == 1)
		return false;

	graph_extractor extractor;
	const auto time_limit = config.time_limit;

	status.modified_stack.pop_back();

	for (size_t node = 0; node < status.remaining_nodes; node++) {
		recursive_graph.setPartitionIndex(node, recursive_comp_map[node] + 2);
	}

	for (size_t i = 0; i < comp_count; i++) {
		if (t.elapsed() > time_limit) {
			timeout = true;
			break;
		}

		recursive_local_mapping.clear();
		graph_access G;
		extractor.extract_block(recursive_graph, G, i + 2, recursive_local_mapping);

		config.time_limit = time_limit - t.elapsed();

		branch_and_reduce_algorithm br_alg(G, config, true);
        cout_handler::disable_cout();
		timeout = !br_alg.run_branch_reduce();
        cout_handler::enable_cout();

		if (timeout)
			break;

		br_alg.apply_branch_reduce_solution(G);

		forall_nodes(G, node) {
			recursive_graph.setPartitionIndex(recursive_local_mapping[node], G.getPartitionIndex(node));
		} endfor
	}

	config.time_limit = time_limit;

	if (timeout) {
		return false;
	}

	forall_nodes(recursive_graph, node) {
		if (recursive_graph.getPartitionIndex(node) == 1)
			set<false>(recursive_mapping[node], IS_status::included);
	} endfor

	return true;
}

void branch_and_reduce_algorithm::branch_reduce_single_component() {
	best_weight = 0;

	if (status.n == 0) {
		return;
	}
	if (status.n == 1) {
		set<false>(0, IS_status::included);
		std::cout << (get_current_is_weight() + status.is_weight + status.reduction_offset) << "," << t.elapsed() << std::endl;
		return;
	}

    resize(status.n);
	sized_vector<NodeID> node_order(status.n);
    node_order.set_size(status.n);
	std::iota(node_order.begin(), node_order.end(), 0);
	std::sort(node_order.begin(), node_order.end(), [&](const NodeID lhs, const NodeID rhs) {
		return deg(lhs) > deg(rhs) || (deg(lhs) == deg(rhs) && status.weights[lhs] > status.weights[rhs]);
	});

	if (status.n > ILS_SIZE_LIMIT)
		compute_ils_pruning_bound();

	size_t i = 0;
	while (i < status.n) {
		if (t.elapsed() > config.time_limit && best_weight != 0) {
			timeout = true;
			break;
		}

		NodeID branch_node = node_order[i];
		if (i == status.n - 1) {
			if (status.node_status[branch_node] == IS_status::not_set)
				set<false>(branch_node, IS_status::included);

			update_best_solution();
			reverse_branching();
			i = status.branching_stack.back().pos;
			continue;
		}

		bool improvement_possible = status.is_weight + status.reduction_offset + compute_cover_pruning_bound() > best_weight;
		if (!improvement_possible) {
            //no improvement_possible --> backtrack
            if (!status.branching_stack.empty() && status.branching_stack.back().node == branch_node) {
                status.branching_stack.pop_back();
                unset(branch_node);

                if (i == 0)
                    break;
            }

            reverse_branching();
            i = status.branching_stack.back().pos;
            continue;
		}

        //improvement_possible
        if (status.node_status[branch_node] == IS_status::not_set) {
            // first time reaching this node in this branch --> include it to IS
            set(branch_node, IS_status::included, false);
            status.branching_stack.emplace_back(branch_node, i);
        } else if (!status.branching_stack.empty() && status.branching_stack.back().node == branch_node) {
            if (status.node_status[branch_node] == IS_status::included) {
                // second time reaching this node in this branch --> exclude it from IS; was included before
                status.node_status[branch_node] = IS_status::excluded;
                status.is_weight -= status.weights[branch_node];
                add_next_level_neighborhood(branch_node);
            } else {
                // third time reaching this node in this branch --> backtack
                status.branching_stack.pop_back();
                unset(branch_node);

                if (i == 0)
                    break;

                reverse_branching();
                i = status.branching_stack.back().pos;
                continue;
            }
        } else {
            i++;
            continue;
        }

		size_t old_n = status.n;
		reduce_graph_internal();

        status.modified_stack.push_back(BRANCHING_TOKEN);
		if (old_n < status.n) {
            node_order.resize(status.n);
            std::iota(node_order.begin() + old_n, node_order.begin() + status.n, old_n);
            std::sort(node_order.begin() + old_n, node_order.begin() + status.n, [&](const NodeID lhs, const NodeID rhs) {
                return deg(lhs) > deg(rhs) || (deg(lhs) == deg(rhs) && status.weights[lhs] > status.weights[rhs]);
            });

            resize(status.n);
        }

		if (status.remaining_nodes > SPLIT_CC_LIMIT && branch_reduce_recursive()) {
			status.modified_stack.push_back(BRANCHING_TOKEN);
			update_best_solution();
			reverse_branching();
			i = status.branching_stack.back().pos;
		}
		else {
			i++;
		}

	}

	restore_best_local_solution();
}

graph_access& branch_and_reduce_algorithm::kernelize() {
    initial_reduce();
    build_global_graph_access();
    return global_graph;
}

bool branch_and_reduce_algorithm::run_branch_reduce() {
    t.restart();
    initial_reduce();

    cout_handler::disable_cout();
    std::cout << "%reduction_nodes " << global_status.remaining_nodes << "\n";
    std::cout << "%reduction_offset " << global_status.is_weight + global_status.reduction_offset << "\n";
    std::cout << "%reduction_time " << t.elapsed() << "\n";
    cout_handler::enable_cout();

    kernelization_time = t.elapsed();
    best_is_time = t.elapsed();
    if (global_status.remaining_nodes == 0) {
        std::cout << get_current_is_weight() << "," << t.elapsed() << std::endl;
        max_min_kernel_comp = 0;
        restore_best_global_solution();
        return true;
    }

    build_global_graph_access();

    std::vector<int> comp_map(global_status.remaining_nodes, 0);
    size_t comp_count = strongly_connected_components().strong_components(global_graph, comp_map);

    //std::cout << "%components " << comp_count << "\n";

    for (size_t node = 0; node < global_status.remaining_nodes; node++) {
        global_graph.setPartitionIndex(node, comp_map[node]);
    }

    std::vector<size_t> comp_size(comp_count, 0);
    for (auto comp_id : comp_map) {
        comp_size[comp_id]++;
    }

    total_ils_node_count = 0;
    std::vector<size_t> comp_idx(comp_count, 0);
    for (size_t i = 0; i < comp_count; i++) {
        comp_idx[i] = i;
        if (comp_size[i] > ILS_SIZE_LIMIT) {
            total_ils_node_count += comp_size[i];
        }
    }

    std::sort(comp_idx.begin(), comp_idx.end(), [&comp_size](const size_t lhs, const size_t rhs) { return comp_size[lhs] < comp_size[rhs]; });

    //std::cout << "%max_component " << comp_size[comp_idx.back()] << "\n";
    max_min_kernel_comp = comp_size[comp_idx.back()];

    graph_extractor extractor;

    buffers.resize(7, sized_vector<NodeID>(global_status.n));

    for (size_t i : comp_idx) {
        if (t.elapsed() > config.time_limit) {
            timeout = true;
            break;
        }

        //std::cout << "%connected component " << i << ":  " << comp_size[i] << std::endl;

        local_mapping.clear();
        graph_access G;
        extractor.extract_block(global_graph, G, i, local_mapping);
        local_graph = &G;

        status = graph_status(*local_graph);

        set_local_reductions();

        branch_reduce_single_component();

        for (size_t node = 0; node < local_mapping.size(); node++) {
            // if (status.node_status[node] == IS_status::not_set || status.node_status[node] == IS_status::folded) {
            // 	std::cerr << "error: node is not_set / folded after restore" << std::endl;
            // 	exit(1);
            // }

            global_status.node_status[global_mapping[local_mapping[node]]] = status.node_status[node];
        }

        global_status.is_weight += best_weight;
        global_status.remaining_nodes -= status.n;
        /* local_graph.reset(); */
        local_graph = nullptr;
    }

    if (timeout) {
        fill_global_greedy();
    }

    restore_best_global_solution();
    return !timeout;
}


void branch_and_reduce_algorithm::update_best_solution() {
	NodeWeight current_weight = status.is_weight + status.reduction_offset;
	if (current_weight > best_weight) {
		best_solution_status = status;
		best_weight = current_weight;
        best_is_time = t.elapsed();
		is_ils_best_solution = false;

		std::cout << (get_current_is_weight() + current_weight) << "," << t.elapsed() << std::endl;
	}
}

void branch_and_reduce_algorithm::undo_blow_up() {
    while (status.remaining_nodes > min_kernel) {
        NodeID node = status.modified_stack.back();
        status.modified_stack.pop_back();

        if (status.node_status[node] == IS_status::folded) {
            auto type = status.folded_stack.back();
            status.folded_stack.pop_back();
            status.transformations[local_transformation_map[type]]->restore(this);
        } else {
            unset(node);
        }
    }
}

void branch_and_reduce_algorithm::reverse_branching() {
	// discard topmost branching token
	if (!status.modified_stack.empty()) {
		status.modified_stack.pop_back();
	}
	else {
		return;
	}

	while (!status.modified_stack.empty() && !is_token(status.modified_stack.back())) {
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();

		if (status.node_status[node] == IS_status::folded) {
			//std::cout << "reversed fold" << std::endl;
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();

			status.transformations[local_transformation_map[type]]->restore(this);
			continue;
		}

		unset(node);
	}
}

void branch_and_reduce_algorithm::restore_best_local_solution() {
	if (is_ils_best_solution) {
		for (size_t node = 0; node < status.n; node++) {
			if (local_graph->getPartitionIndex(node) == 1) {
				status.node_status[node] = IS_status::included;
			}
			else {
				status.node_status[node] = IS_status::excluded;
			}
		}
		return;
	}

	status = best_solution_status;
	status.modified_stack.pop_back();

	while (!status.modified_stack.empty()) {
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();

		if (node == BRANCHING_TOKEN) {
			status.graph.restore_node(status.branching_stack.back().node);
			status.branching_stack.pop_back();
			continue;
		}

		if (status.node_status[node] == IS_status::folded) {
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->apply(this);
		}
		else {
			status.graph.restore_node(node);
		}
	}
}

void branch_and_reduce_algorithm::restore_best_global_solution() {
    status = std::move(global_status);
	status.modified_stack.pop_back();

    while (!status.modified_stack.empty()) {
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();

        if (status.node_status[node] == IS_status::folded) {
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[global_transformation_map[type]]->apply(this);
		}
		else {
			status.graph.restore_node(node);
		}
    }
}

NodeWeight branch_and_reduce_algorithm::get_current_is_weight() const {
	return global_status.is_weight + global_status.reduction_offset;
}

void branch_and_reduce_algorithm::build_global_graph_access() {
	global_mapping.resize(global_status.remaining_nodes, 0);
	std::swap(status, global_status);
	build_graph_access(global_graph, global_mapping);
	std::swap(status, global_status);
}

void branch_and_reduce_algorithm::build_graph_access(graph_access & G, std::vector<NodeID> & reverse_mapping) const {
	std::vector<NodeID> mapping(status.graph.size(), UINT_MAX);
	size_t edge_count = 0;

	// Get number of edges and reorder nodes
	size_t node_counter = 0;
	for (NodeID node = 0; node < status.graph.size(); ++node) {
		if (status.node_status[node] == IS_status::not_set) {
            edge_count += status.graph[node].size();

			mapping[node] = node_counter;
			reverse_mapping[node_counter] = node;
			node_counter++;
		}
	}

	// Create the adjacency array
	std::vector<int> xadj(status.remaining_nodes + 1);
	std::vector<int> adjncy(edge_count + 1);
	size_t adjncy_counter = 0;

	for (size_t i = 0; i < status.remaining_nodes; ++i) {
		xadj[i] = adjncy_counter;

		for (auto neighbor : status.graph[reverse_mapping[i]]) {
			if (mapping[neighbor] == i) continue;
			if (mapping[neighbor] == UINT_MAX) continue;
			adjncy[adjncy_counter++] = mapping[neighbor];
		}

		std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
	}
	xadj[status.remaining_nodes] = adjncy_counter;

	// Build the graph
	G.build_from_metis(status.remaining_nodes, &xadj[0], &adjncy[0]);

	forall_nodes(G, node) {
		G.setNodeWeight(node, status.weights[reverse_mapping[node]]);
	} endfor
}

void branch_and_reduce_algorithm::build_induced_neighborhood_subgraph(graph_access& G, NodeID source_node) {
	buffers[0].clear();
	set_1.clear();

	for (NodeID neighbor : status.graph[source_node]) {
		buffers[0].push_back(neighbor);
		set_1.add(neighbor);
	}

	build_induced_subgraph(G, buffers[0], set_1, buffers[1]);
}

void branch_and_reduce_algorithm::build_induced_subgraph(graph_access& G, const sized_vector<NodeID>& nodes, const fast_set& nodes_set, sized_vector<NodeID>& reverse_mapping) {
	size_t edge_count = 0;

	for (size_t i = 0; i < nodes.size(); i++) {
		NodeID node = nodes[i];
		for (auto neighbor : status.graph[node]) {
			if (nodes_set.get(neighbor))
				edge_count++;
		}

        reverse_mapping[node] = i;
	}

	G.start_construction(nodes.size(), edge_count);

	for (NodeID node : nodes) {
		NodeID new_node = G.new_node();
		G.setNodeWeight(new_node, status.weights[node]);

		for (auto neighbor : status.graph[node]) {
			if (nodes_set.get(neighbor)) {
				G.new_edge(new_node, reverse_mapping[neighbor]);
			}
		}
	}

	G.finish_construction();
}

void branch_and_reduce_algorithm::reverse_reduction(graph_access & G, graph_access & reduced_G, std::vector<NodeID> & reverse_mapping) {
	//build_global_graph_access();

	forall_nodes(reduced_G, node) {
		if (reduced_G.getPartitionIndex(node) == 1) {
			status.node_status[reverse_mapping[node]] = IS_status::included;
		} else {
			status.node_status[reverse_mapping[node]] = IS_status::excluded;
		}
	} endfor

	global_status = std::move(status);
	restore_best_global_solution();
	apply_branch_reduce_solution(G);
}

void branch_and_reduce_algorithm::apply_branch_reduce_solution(graph_access & G) {
	forall_nodes(G, node) {
		if (status.node_status[node] == IS_status::included) {
			G.setPartitionIndex(node, 1);
		} else {
			G.setPartitionIndex(node, 0);
		}
	} endfor
}
