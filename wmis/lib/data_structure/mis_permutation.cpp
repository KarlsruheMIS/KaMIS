/**
 * mis_permutation.cpp
 * Purpose: Data structure used for the local search algorithm.
 *          The structure itself is a permutation of all nodes divided into three blocks.
 *          The first block contains the solution nodes.
 *          The second block contains free nodes.
 *          The third block contains all remaining (non-free, non-solution) nodes.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "mis_permutation.h"

#include <stdio.h>
#include <algorithm>

#include "macros_assertions.h"
#include "data_structure/operation_log.h"

mis_permutation::mis_permutation() {

}

mis_permutation::~mis_permutation() {
}

void mis_permutation::construct(graph_access & G) {
	G_ptr = &G;
    inconsistencies = 0;
    solution_size = 0;
    free_size = 0;
    total_size = G.number_of_nodes();
    nodes.clear();
    tightness.clear();
    position.clear();
    nodes.resize(total_size);
    tightness.resize(total_size);
    position.resize(total_size);
    onetight_all.init(G.number_of_nodes());

    // Insert solution nodes
    forall_nodes(G, n) {
        nodes[n] = n;
        position[n] = n;
        unsigned int index = G.getPartitionIndex(n);
        // Maybe implement tightness calculations here
        if (index == 1) {
            move_to_solution(n, G);
        } else {
            int tight = calculate_tightness(n, G);
            tightness[n] = tight;
            if (tight == 0) move_to_free(n, G);
            else move_to_non_free(n, G);
            if (tight == 1) onetight_all.insert(n);
        }
    } endfor

}

int mis_permutation::calculate_tightness(NodeID node, graph_access & G) {
    int tightness = 0;
    forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);
        bool target_index = G.getPartitionIndex(target);
        if (target_index == 1) {
            tightness++;
        }
    } endfor

    return tightness;
}

NodeID mis_permutation::get_solution_node(unsigned int i) {
    ASSERT_LT(i, solution_size);
    return nodes[i];
}

NodeID mis_permutation::get_free_node(unsigned int i) {
    ASSERT_LT(i, free_size);
    return nodes[solution_size + i]; 
}

NodeID mis_permutation::get_non_free_node(unsigned int i) {
    ASSERT_LT(i, total_size - (free_size + solution_size));
    return nodes[solution_size + free_size + i]; 
}

NodeID mis_permutation::get_non_solution_node(unsigned int i) {
    ASSERT_LT(i, total_size - solution_size);
    return nodes[solution_size + i]; 
}

unsigned int mis_permutation::get_solution_size() {
    return solution_size;
}

NodeWeight mis_permutation::get_solution_weight() {
	NodeWeight solution_weight = 0;
	for (size_t i = 0; i < solution_size; i++) {
		solution_weight += G_ptr->getNodeWeight(nodes[i]);
	}
	return solution_weight;
}

unsigned int mis_permutation::get_free_size() {
    return free_size;
}

bool mis_permutation::is_maximal() {
    return (free_size == 0);
}

void mis_permutation::add_to_solution(NodeID node, graph_access & G) {

    // Node is already in the solution
    if (is_solution_node(node)) return;

    // Add node to the solution
    G.setPartitionIndex(node, 1);
    move_to_solution(node, G);

    // Update neighbors
    forall_out_edges (G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);

        if (is_solution_node(target)) {
            inconsistencies++;
            printf("Inconsistency");
            exit(-1);
        }
        else {
            if (tightness[target] == 0) move_to_non_free(target, G);
        }

        if (tightness[target] == 1) onetight_all.remove(target);
        tightness[target]++;
        if (tightness[target] == 1) onetight_all.insert(target);

        // if (is_solution_node(target)) G.setPartitionIndex(target, 1);
        // if (!is_solution_node(target)) G.setPartitionIndex(target, 0);
    } endfor

    ASSERT_LT(position[node], solution_size);
    operation_log::instance()->report_insert(node);
}

void mis_permutation::remove_from_solution(NodeID node, graph_access & G) {

    // Node isn't in the solution
    if (!is_solution_node(node)) return;

    // Remove node from the solution
    G.setPartitionIndex(node, 0);
    bool is_free = (tightness[node] == 0);
    if (is_free) move_to_free(node, G);
    else move_to_non_free(node, G);

    // Update neighbors
    forall_out_edges(G, edge, node) {
        NodeID target = G.getEdgeTarget(edge);

        if (tightness[target] == 1) onetight_all.remove(target);
        tightness[target]--;
        if (tightness[target] == 1) onetight_all.insert(target);

        if (is_solution_node(node)) {
            inconsistencies--;
            printf("Inconsistency");
            exit(-1);
        }
        else {
            if (tightness[target] == 0) move_to_free(target, G);
        }

        // if (is_solution_node(target)) G.setPartitionIndex(target, 1);
        // if (!is_solution_node(target)) G.setPartitionIndex(target, 0);
    } endfor

    ASSERT_GEQ(position[node], solution_size);
    operation_log::instance()->report_remove(node);
}

void mis_permutation::move_to_solution(NodeID node, graph_access & G) {
    unsigned int node_position = position[node];

    // Node already in the solution then skip
    if (node_position < solution_size) return;
    // Is the node free?
    if (node_position < solution_size + free_size) {
        swap_nodes(solution_size, node_position);
        free_size--;
    } 
    // Non free node 
    else {
        swap_nodes(solution_size + free_size, position[node]);
        swap_nodes(solution_size, solution_size + free_size);
    }
    solution_size++;

	added_vertices++;
}

void mis_permutation::move_to_free(NodeID node, graph_access & G) {
    unsigned int node_position = position[node];

    // printf("Move %d to free\n", node);
    // Is the node in the solution?
    if (node_position < solution_size) {
        swap_nodes(solution_size - 1, node_position);
        solution_size--;
    }
    // Node already free then skip
    else if (node_position < solution_size + free_size) return;
    // Non free node
    else {
        swap_nodes(solution_size + free_size, position[node]);
    }
    free_size++;
}

void mis_permutation::move_to_non_free(NodeID node, graph_access & G) {
    unsigned int node_position = position[node];

    // printf("Move %d to non free\n", node);
    // Is the node already non free
    if (node_position >= solution_size + free_size) return;
    // Is the node free?
    if (node_position >= solution_size) {
        swap_nodes(solution_size + free_size - 1, node_position);
        free_size--;
    }
    // Solution node
    else {
        swap_nodes(solution_size - 1, position[node]);
        swap_nodes(solution_size + free_size - 1, solution_size - 1);
        solution_size--;
    }
}

void mis_permutation::swap_nodes(unsigned int first_pos, unsigned int second_pos) {
    if (first_pos == second_pos) return;

    NodeID first_node = nodes[first_pos]; 
    NodeID second_node = nodes[second_pos];
    nodes[first_pos] = second_node;
    nodes[second_pos] = first_node;
    position[first_node] = second_pos;
    position[second_node] = first_pos;
}

bool mis_permutation::is_solution_node(NodeID node) {
    return (position[node] < solution_size);
}

bool mis_permutation::is_free_node(NodeID node) {
    return (position[node] >= solution_size && position[node] < free_size + solution_size);
}

bool mis_permutation::is_non_solution_node(NodeID node) {
    return (position[node] >= solution_size + free_size && position[node] < total_size);
}

void mis_permutation::print(bool details) {
    printf("\n");
    printf("********************\n");
    printf("Solution size: %d\n", solution_size);
    printf("Free size: %d\n", free_size);
    printf("Maximal: %d\n", is_maximal());
    if (details) {
        printf("\n");
        printf("Solution nodes: \n");
        for(unsigned int i = 0; i < solution_size; ++i) {
            printf(" %d ", nodes[i]);
        }
        printf("\n");

        printf("Free nodes: \n");
        for(unsigned int i = solution_size; i < solution_size + free_size; ++i) {
            printf(" %d ", nodes[i]);
        }
        printf("\n");

        printf("Remaining nodes: \n");
        for(unsigned int i = solution_size + free_size; i < total_size; ++i) {
            printf(" %d ", nodes[i]);
        }

        printf("\n");
    }

    printf("********************\n");
    printf("\n");
}

void mis_permutation::print_position() {
    printf("\n");

    for(unsigned int i = 0; i < solution_size; ++i) {
        printf("Position: %d\n", position[i]);
    }

    for(unsigned int i = solution_size; i < solution_size + free_size; ++i) {
        printf("Position: %d\n", position[i]);
    }

    for(unsigned int i = solution_size + free_size; i < total_size; ++i) {
        printf("Position: %d\n", position[i]);
    }

    printf("\n");
}

void mis_permutation::print_tightness() {
    printf("\n");

    for(unsigned int i = 0; i < solution_size; ++i) {
        printf("Tightness: %d\n", tightness[i]);
    }

    for(unsigned int i = solution_size; i < solution_size + free_size; ++i) {
        printf("Tightness: %d\n", tightness[i]);
    }

    for(unsigned int i = solution_size + free_size; i < total_size; ++i) {
        printf("Tightness: %d\n", tightness[i]);
    }

    printf("\n");
}

bool mis_permutation::check_permutation() {
    bool solution_correct = true;

    for(unsigned int i = 0; i < solution_size; ++i) {
        NodeID node = nodes[i]; 
        if (tightness[node] != 0) {
            printf("Tightness error in solution\n");
            solution_correct = false;
        }
    }

    for(unsigned int i = solution_size; i < solution_size + free_size; ++i) {
        NodeID node = nodes[i]; 
        if (tightness[node] != 0) {
            printf("Tightness error in free\n");
            solution_correct = false;
        }
    }

    for(unsigned int i = solution_size + free_size; i < total_size; ++i) {
        NodeID node = nodes[i]; 
        if (tightness[node] < 1) {
            printf("Tightness error in non free\n");
            solution_correct = false;
        }
    }

    return solution_correct;
}

bool mis_permutation::check_consistency(graph_access & G) {
    if (inconsistencies > 0) return false;
    forall_nodes(G, node) {
        if (is_solution_node(node) && G.getPartitionIndex(node) != 1) return false;
        if (!is_solution_node(node) && G.getPartitionIndex(node) != 0) return false;
    } endfor
    return true;
}

