/**
 * local_search.cpp
 * Purpose: Apply the local search algorithm to a maximum independent set.
 *
 * The original code from Andrade et. al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "ils/local_search.h"

#include <algorithm>

#include "random_functions.h"
#include "bucket_array.h"

constexpr NodeID local_search::INVALID_NODE;

local_search::local_search(const MISConfig& config) {
	sort_freenodes = config.sort_freenodes;
}

local_search::~local_search() {

}

void local_search::preprocess_graph(graph_access & G) {
    perm.construct(G);
	onetight.clear();
	neighbors.clear();
	onetight.resize(G.getMaxDegree());
	neighbors.resize(G.getMaxDegree());
    build_candidates(G);
    onetight.clear();
    neighbors.clear();
    onetight.resize(G.getMaxDegree());
    neighbors.resize(G.getMaxDegree());
}

void local_search::preprocess_graph_with_candidates(graph_access & G, std::vector<NodeID> cand, unsigned int cand_size) {
    perm.construct(G);
	candidates.clear();
	candidate_swaps.clear();
    insert_candidates(G, cand, cand_size);
    onetight.clear();
    neighbors.clear();
    onetight.resize(G.getMaxDegree());
    neighbors.resize(G.getMaxDegree());
}

void local_search::force(graph_access & G, unsigned int k) {
    for (unsigned int i = 0; i < k; ++i) {
        unsigned int upper_limit = G.number_of_nodes() - perm.get_solution_size();
        unsigned int random = random_functions::nextInt(0, upper_limit - 1);
        NodeID x = perm.get_non_solution_node(random);
        force_node(G, x);
        make_maximal(G);
        direct_improvement(G, true, x);
    }
}

void local_search::force_node(graph_access & G, NodeID node) {
    forall_out_edges(G, edge, node) {
        NodeID w = G.getEdgeTarget(edge);
        if (perm.is_solution_node(w)) {
            perm.remove_from_solution(w, G);
        }
    } endfor

    perm.add_to_solution(node, G);
    if (!candidates.contains(node)) add_candidate(node, G);
}

void local_search::simple_improvement(graph_access & G, bool forced, NodeID forced_node) {
	// TODO: adapt to weighted case
    //
    //if (candidates.empty()) build_candidates(G);

    //while(1) {
    //    NodeID x;
    //    if (candidates.empty()) {
    //        if (forced) {
    //            x = forced_node;
    //            forced = false;
    //        }
    //        else break;
    //    } else {
    //        x = candidates.deleteMax();
    //        if (!perm.is_solution_node(x)) continue;
    //        if (forced && x == forced_node) continue;
    //    }

    //    build_onetight(x, G);

    //    if (onetight_size < 2) continue;
    //    int solutions_found = 0;
    //    NodeID improv_v = 0;
    //    NodeID improv_w = 0;
    //    perm.remove_from_solution(x, G);

    //    for (int i = onetight_size - 1; i >= 0; i--) {
    //        if (i == 0) {
    //            if (solutions_found == 0) continue;
    //            if (onetight_size == 2) continue;
    //        }
    //        NodeID v = onetight[i];
    //        perm.add_to_solution(v, G);
    //        unsigned int remaining_free = perm.get_free_size();
    //        if (remaining_free > 0) {
    //            solutions_found += remaining_free;
    //            int random = random_functions::nextInt(0, perm.get_free_size() - 1);
    //            improv_v = v;
    //            improv_w = perm.get_free_node(random);
    //        }
    //        perm.remove_from_solution(v, G);
    //    }

    //    if (solutions_found == 0) {
    //        perm.add_to_solution(x, G);
    //    } else {
    //        perm.add_to_solution(improv_v, G);
    //        perm.add_to_solution(improv_w, G);
    //        if (!candidates.contains(improv_v)) add_candidate(improv_v, G);
    //        if (!candidates.contains(improv_w)) add_candidate(improv_w, G);

    //        if (!perm.is_maximal()) make_maximal(G);  
    //        // Incremental
    //        update_candidates(x, G);
    //    }
    //    ASSERT_TRUE(perm.check_permutation());
    //} 
}

void local_search::direct_improvement(graph_access & G, bool forced, NodeID forced_node) {
	// TODO: figure out what to do with a fored node
	forced = false;
    
    if (candidates.empty()) build_candidates(G);

    while(1) {
        NodeID x;
		Swap_1_2 swap;
        if (candidates.empty()) {
            if (forced) {
                x = forced_node;
				swap = find_best_swap(x, G);
                forced = false;
            }
            else break;
        } else {
            x = candidates.deleteMax();
			swap = std::move(candidate_swaps[x]);
			candidate_swaps.erase(x);

            if (!perm.is_solution_node(x)) continue;
            if (forced && x == forced_node) continue;
        }

        perm.remove_from_solution(x, G);

		if (swap.first != INVALID_NODE) {
			perm.add_to_solution(swap.first, G);
			if (!candidates.contains(swap.first)) add_candidate(swap.first, G);
		}

		if (swap.second != INVALID_NODE) {
			perm.add_to_solution(swap.second, G);
			if (!candidates.contains(swap.second)) add_candidate(swap.second, G);
		}

        if (!perm.is_maximal()) make_maximal(G);  
            
        // Incremental
        update_candidates(x, G);
        
        ASSERT_TRUE(perm.check_permutation());
    }
}

void local_search::insert_candidates(graph_access & G, std::vector<NodeID> cand, unsigned int num_cand) {
    for (unsigned int i = 0; i < num_cand; ++i) {
        if (!candidates.contains(cand[i])) {
            add_candidate(cand[i], G);
        }
    }
}

void local_search::make_maximal(graph_access & G) {
	if (sort_freenodes) {
		std::vector<NodeID> free_ndoes(perm.get_free_size());
		for (size_t i = 0; i < perm.get_free_size(); i++) {
			free_ndoes[i] = perm.get_free_node(i);
		}

		sort_by_weight(G, free_ndoes.begin(), free_ndoes.end());

		for (auto n : free_ndoes) {
			if (perm.is_free_node(n)) {
				perm.add_to_solution(n, G);
				if (!candidates.contains(n)) add_candidate(n, G);
			}
		}

		return;
	}
	
    while(perm.get_free_size() > 0) {
        int random = random_functions::nextInt(0, perm.get_free_size() - 1);
        NodeID free_node = perm.get_free_node(random);
        perm.add_to_solution(free_node, G);
        if (!candidates.contains(free_node)) add_candidate(free_node, G);
    }
}

void local_search::build_onetight(NodeID node, graph_access & G) {
    onetight_size = 0; 
    forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        int target_tight = perm.get_tightness(target);
        if (target_tight == 1) {
            onetight[onetight_size++] = target;
        }
    } endfor
}

void local_search::update_swaps(NodeID node, graph_access & G) {
	// check if any neighboring nodes are contained in any (1,2)-swaps for other
	// candidate nodes and update them if this is the case

	// TODO: optimize runtime -> use hashmap / vector to map from onetight nodes to their candidate node
	forall_out_edges(G, edge, node) {
		NodeID old_onetight = G.getEdgeTarget(edge);
		if (perm.get_tightness(old_onetight) == 2) {
			// former onetight node might be contained in a (1,2)-swap
			forall_out_edges(G, edge2, old_onetight) {
				NodeID candidate = G.getEdgeTarget(edge2);
				if (candidate == node) continue;

				auto iter = candidate_swaps.find(candidate);
				if (iter != candidate_swaps.end()) {
					if (iter->second.first == old_onetight || iter->second.second == old_onetight) {
						// conflict found
						auto swap = find_best_swap(candidate, G);
						int gain = insert_swap(candidate, swap, G);

						if (gain > 0) {
							candidates.changeKey(candidate, gain);
						} else {
							candidates.deleteNode(candidate);
							candidate_swaps.erase(candidate);
						}
					}
					// only one other candidate node other than "node"
					break;
				}
			} endfor
		}
	} endfor
}

void local_search::print_onetight() {
    for (unsigned int i = 0; i < onetight_size; ++i) {
        printf("Node: %d\n", onetight[i]);
    }
}

void local_search::build_neighbors(NodeID node, graph_access & G) {
    neighbors_size = 0;
    forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        neighbors[neighbors_size++] = target;
    } endfor
}

void local_search::build_candidates(graph_access & G) {
    candidates.clear();
	candidate_swaps.clear();

    unsigned int solution_size = perm.get_solution_size();
    for (unsigned int i = 0; i < solution_size; ++i) {
		add_candidate(perm.get_solution_node(i), G);
    }
}

void local_search::add_candidate(NodeID node, graph_access & G) {
	Swap_1_2 swap = find_best_swap(node, G);
	int gain = insert_swap(node, swap, G);
	if (gain > 0) {
		candidates.insert(node, gain);
	}
	update_swaps(node, G);
}

int local_search::insert_swap(NodeID node, Swap_1_2 swap, graph_access & G) {
	int gain = -G.getNodeWeight(node);
	if (swap.first != INVALID_NODE) gain += G.getNodeWeight(swap.first);
	if (swap.second != INVALID_NODE) gain += G.getNodeWeight(swap.second);

	if (gain > 0) candidate_swaps[node] = std::move(swap);
	
	return gain;
}

void local_search::update_candidates(NodeID node, graph_access & G) {
    forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        // Skip if neighbor is not 1-tight
        if (perm.get_tightness(target) != 1) continue;
        forall_out_edges(G, target_edge, target) {
            NodeID candidate = G.getEdgeTarget(target_edge);
            if (perm.is_solution_node(candidate)) {
                if (!candidates.contains(candidate)) add_candidate(candidate, G);
				else update_candidate(candidate, target, G);
                // There can only be one valid candidate
                break;
            }
        } endfor
    } endfor
}

void local_search::update_candidate(NodeID node, NodeID one_tight, graph_access & G) {
	build_onetight(node, G);

	sort_by_weight(G, onetight.begin(), std::next(onetight.begin(), onetight_size));

	Swap_1_2 new_swap = { one_tight, INVALID_NODE };
	int new_gain = G.getNodeWeight(one_tight) - G.getNodeWeight(node);

	for (size_t i = 0; i < onetight_size; i++) {
		if (onetight[i] == one_tight) continue;

		// check if 1-tight neighbors are connected
		bool connected = false;
		forall_out_edges(G, edge, one_tight) {
			connected = G.getEdgeTarget(edge) == onetight[i];
			if (connected) break;
		} endfor

		if (connected) continue;

		new_gain += G.getNodeWeight(onetight[i]);
		new_swap.second = onetight[i];
		break;
	}
	

	auto iter = candidate_swaps.find(node);
	if (iter != candidate_swaps.end()) {
		if (candidates.getKey(node) >= new_gain) return; // old swap still good

		iter->second = std::move(new_swap);
		candidates.changeKey(node, new_gain);
	} else {
		// node was not in candidate list
		candidate_swaps[node] = std::move(new_swap);
		candidates.insert(node, new_gain);
	}
}

void local_search::print_permutation() {
    perm.print(0);
    perm.check_permutation();
}

local_search::Swap_1_2 local_search::find_best_swap(NodeID node, graph_access & G) {
	build_onetight(node, G);

	sort_by_weight(G, onetight.begin(), std::next(onetight.begin(), onetight_size));

	Swap_1_2 best_swap{ INVALID_NODE, INVALID_NODE };
	NodeWeight best_gain = 0;

	for (size_t i = 0; i < onetight_size && 2 * G.getNodeWeight(onetight[i]) > best_gain; i++) {
		NodeWeight first_gain = G.getNodeWeight(onetight[i]);

		if (best_swap.first == INVALID_NODE) {
			// possible (1,1)-swap
			best_swap.first = onetight[i];
			best_gain = first_gain;
		}

		for (size_t j = i + 1; j < onetight_size && G.getNodeWeight(onetight[j]) + first_gain > best_gain; j++) {
			// check if 1-tight neighbors are connected
			bool connected = false;
			forall_out_edges(G, edge, onetight[i]) {
				connected = G.getEdgeTarget(edge) == onetight[j];
				if (connected) break;
			} endfor

			if (connected) continue;

			if (first_gain + G.getNodeWeight(onetight[j]) > best_gain) {
				best_gain = first_gain + G.getNodeWeight(onetight[j]);
				best_swap = { onetight[i], onetight[j] };
				// loop stops
			}
		}
	}

	return best_swap;
}

void local_search::sort_by_weight(graph_access & G, std::vector<NodeID>::iterator begin, std::vector<NodeID>::iterator end) {
	const size_t size = end - begin;
	NodeWeight max_weight = 0;
	std::vector<NodeWeight> input;
	input.reserve(size);

	for (auto iter = begin; iter != end; iter++) {
		input.push_back(*iter);
		if (max_weight < G.getNodeWeight(*iter)) {
			max_weight = G.getNodeWeight(*iter);
		}
	}

	std::vector<size_t> count(max_weight + 1);

	// count occurences
	for (auto node : input) {
		count[G.getNodeWeight(node)]++;
	}

	size_t count_tmp;
	size_t total = 0;
	// prefix sum
	for (int weight = max_weight; weight >= 0; weight--) {
		count_tmp = count[weight];
		count[weight] = total;
		total += count_tmp;
	}

	// extract sorted nodes
	for (auto node : input) {
		*std::next(begin, count[G.getNodeWeight(node)]) = node;
		count[G.getNodeWeight(node)]++;
	}
}
