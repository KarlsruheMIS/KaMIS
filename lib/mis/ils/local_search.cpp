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

local_search::local_search() {

}

local_search::~local_search() {

}

void local_search::preprocess_graph(graph_access & G) {
    perm.construct(G); 
    build_candidates(G);
    onetight.clear();
    neighbors.clear();
    onetight.resize(G.getMaxDegree());
    neighbors.resize(G.getMaxDegree());
}

void local_search::preprocess_graph_with_candidates(graph_access & G, std::vector<NodeID> cand, unsigned int cand_size) {
    perm.construct(G); 
    candidates.init(G.number_of_nodes());
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
    if (!candidates.contains(node)) candidates.insert(node);
}

void local_search::simple_improvement(graph_access & G, bool forced, NodeID forced_node) {
    
    if (candidates.is_empty()) build_candidates(G);

    while(1) {
        NodeID x;
        if (candidates.is_empty()) {
            if (forced) {
                x = forced_node;
                forced = false;
            }
            else break;
        } else {
            x = candidates.remove_random();
            if (!perm.is_solution_node(x)) continue;
            if (forced && x == forced_node) continue;
        }

        build_onetight(x, G);

        if (onetight_size < 2) continue;
        int solutions_found = 0;
        NodeID improv_v = 0;
        NodeID improv_w = 0;
        perm.remove_from_solution(x, G);

        for (int i = onetight_size - 1; i >= 0; i--) {
            if (i == 0) {
                if (solutions_found == 0) continue;
                if (onetight_size == 2) continue;
            }
            NodeID v = onetight[i];
            perm.add_to_solution(v, G);
            unsigned int remaining_free = perm.get_free_size();
            if (remaining_free > 0) {
                solutions_found += remaining_free;
                int random = random_functions::nextInt(0, perm.get_free_size() - 1);
                improv_v = v;
                improv_w = perm.get_free_node(random);
            }
            perm.remove_from_solution(v, G);
        }

        if (solutions_found == 0) {
            perm.add_to_solution(x, G);
        } else {
            perm.add_to_solution(improv_v, G);
            perm.add_to_solution(improv_w, G);
            if (!candidates.contains(improv_v)) candidates.insert(improv_v);
            if (!candidates.contains(improv_w)) candidates.insert(improv_w);

            if (!perm.is_maximal()) make_maximal(G);  
            // Incremental
            update_candidates(x, G);
        }
        ASSERT_TRUE(perm.check_permutation());
    } 
}

void local_search::direct_improvement(graph_access & G, bool forced, NodeID forced_node) {
    
    if (candidates.is_empty()) build_candidates(G);

    while(1) {

        NodeID x;
        if (candidates.is_empty()) {
            if (forced) {
                x = forced_node;
                forced = false;
            }
            else break;
        } else {
            x = candidates.remove_random();
            if (!perm.is_solution_node(x)) continue;
            if (forced && x == forced_node) continue;
        }

        // Build L(x)
        build_onetight(x, G);

        if (onetight_size < 2) continue;
        int solutions_found = 0;
        NodeID improv_v = 0;
        NodeID improv_w = 0;

        for (int i = onetight_size - 1; i >= 0; i--) {
            NodeID v = onetight[i];
            // No need to last neighbor if no swap was found so far
            if (i == 0) {
                // No swap can be found
                if (solutions_found == 0) break;
                // Swap would have been found already
                if (onetight_size == 2) break;
            }
            // Build A(v)
            build_neighbors(v, G);
            unsigned int v_pos = 0;
            unsigned int x_pos = 0;

            while(1) {
                // No more neighbors to check 
                if (x_pos >= onetight_size) break; 

                NodeID x_neighbor = onetight[x_pos];
                // Skip if the neighbor is v itself
                if (x_neighbor == v) {
                    x_pos++;
                    continue;
                }

                // Reached the end of A(v) so all remaining nodes are possible improvements
                if (v_pos >= neighbors_size) {
                    solutions_found++;
                    if (solutions_found == 1 || random_functions::nextInt(1, solutions_found) == 1) {
                        improv_v = v;
                        improv_w = x_neighbor;
                    }
                    x_pos++;
                    continue;
                }

                NodeID v_neighbor = neighbors[v_pos];
                // Skip if the neighbor is x itself
                if (v_neighbor == x) {
                    v_pos++;
                    continue;
                }

                // Skip if non candidate neighbor
                if (v_neighbor < x_neighbor) {
                    v_pos++;
                    continue;
                }

                // Skip if both neighbors are the same
                if (v_neighbor == x_neighbor) {
                    v_pos++;
                    x_pos++;
                    continue;
                }

                // Candidate found but still continue looking
                if (v_neighbor > x_neighbor) {
                    solutions_found++;
                    if (random_functions::nextInt(1, solutions_found) == 1) {
                        improv_v = v;
                        improv_w = x_neighbor;
                    }
                    x_pos++;
                }
            }
            
        }

        if (solutions_found > 0) {
            perm.remove_from_solution(x, G);
            perm.add_to_solution(improv_v, G);
            perm.add_to_solution(improv_w, G);
            if (!candidates.contains(improv_v)) candidates.insert(improv_v);
            if (!candidates.contains(improv_w)) candidates.insert(improv_w);

            if (!perm.is_maximal()) make_maximal(G);  
            
            // Incremental
            update_candidates(x, G);
        }
        ASSERT_TRUE(perm.check_permutation());
    }
}

void local_search::insert_candidates(graph_access & G, std::vector<NodeID> cand, unsigned int num_cand) {
    for (unsigned int i = 0; i < num_cand; ++i) {
        if (!candidates.contains(cand[i])) {
            candidates.insert(cand[i]);
        }
    }
}

void local_search::make_maximal(graph_access & G) {
    while(perm.get_free_size() > 0) {
        int random = random_functions::nextInt(0, perm.get_free_size() - 1);
        NodeID free_node = perm.get_free_node(random);
        perm.add_to_solution(free_node, G);
        if (!candidates.contains(free_node)) candidates.insert(free_node);
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
    candidates.init(G.number_of_nodes());
    unsigned int solution_size = perm.get_solution_size();
    for (unsigned int i = 0; i < solution_size; ++i) {
        candidates.insert(perm.get_solution_node(i));
    }
}

void local_search::update_candidates(NodeID node, graph_access & G) {
    forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        // Skip if neighbor is not 1-tight
        if (perm.get_tightness(target) != 1) continue;
        forall_out_edges(G, target_edge, target) {
            NodeID candidate = G.getEdgeTarget(target_edge);
            if (perm.get_position(candidate) < perm.get_solution_size()) {
                if (!candidates.contains(candidate)) candidates.insert(candidate);
                // There can only be one valid candidate
                break;
            }
        } endfor
    } endfor
}

void local_search::print_permutation() {
    perm.print(0);
    perm.check_permutation();
}

