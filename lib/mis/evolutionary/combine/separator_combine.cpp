/**
 * separator_combine.h
 * Purpose: Combine two individuals based on node separators.
 *
 *****************************************************************************/

#include "separator_combine.h"

#include "data_structure/priority_queues/bucket_array.h"

separator_combine::separator_combine() {
    
}

separator_combine::~separator_combine() {
    
}

void separator_combine::combine(MISConfig & config, graph_access & G, separator_pool *pool, individuum_mis & first, individuum_mis & second, individuum_mis & out_first, individuum_mis & out_second) {
    // Generate node seperator V = S u A u B
    // Retrieve random separator from the pool and apply it to G
    population_mis pop;
    G.resizeSecondPartitionIndex(G.number_of_nodes());

    apply_separator_kahip(config, G, pool);

    // Build intersections
    // Use both partition indices to build two offsprings at once
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) {
            // A and I1 on first index
            G.setPartitionIndex(node, first.solution[node]);
            // A and I2 on second index
            G.setSecondPartitionIndex(node, second.solution[node]);
        }
        else if (G.getPartitionIndex(node) == 1) {
            // B and I2 on first index
            G.setPartitionIndex(node, second.solution[node]);
            // B and I1 on second index
            G.setSecondPartitionIndex(node, first.solution[node]);
        }
        else if (G.getPartitionIndex(node) == 2) {
            G.setPartitionIndex(node, 2);
            G.setSecondPartitionIndex(node, 2);
        }
        else printf("Error in partitioning separator.");
    } endfor

    std::vector<NodeID> node_candidates;
    build_separator_candidates(config, G, node_candidates);

    if (!config.optimize_candidates) local.preprocess_graph(G);
    else local.preprocess_graph_with_candidates(G, node_candidates, node_candidates.size());
    local.make_maximal(G);
    local.direct_improvement(G);

    NodeID *solution_first = new NodeID[G.number_of_nodes()];
    unsigned int solution_first_size = pop.create_solution(G, solution_first);
    out_first.solution = solution_first;
    out_first.solution_size = solution_first_size;
    ASSERT_TRUE(pop.is_mis(config, G, out_first));
    
    // Swap the partition index and repeat this process to obtain 
    // a second offspring
    forall_nodes(G, node) {
        G.setPartitionIndex(node, G.getSecondPartitionIndex(node));
    } endfor

    node_candidates.clear();
    build_separator_candidates(config, G, node_candidates);

    if (!config.optimize_candidates) local.preprocess_graph(G);
    else local.preprocess_graph_with_candidates(G, node_candidates, node_candidates.size());
    local.make_maximal(G);
    local.direct_improvement(G);

    NodeID *solution_second = new NodeID[G.number_of_nodes()];
    unsigned int solution_second_size = pop.create_solution(G, solution_second);
    out_second.solution = solution_second;
    out_second.solution_size = solution_second_size;
    ASSERT_TRUE(pop.is_mis(config, G, out_second));
}

void separator_combine::build_separator_candidates(MISConfig & config, graph_access & G, std::vector<NodeID> & candidates) {
    // Initialize a bucket queue for the separator
    bucket_array *buckets = new bucket_array(G.number_of_nodes());
    forall_nodes(G, node) {
        // Limit the bucket queue to S
        if (G.getPartitionIndex(node) != 2) buckets->remove(node);
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            // If their is a neighbor in the solution, remove the node
            if (G.getPartitionIndex(target) == 1) {
                G.setPartitionIndex(node, 0);
                buckets->remove(node);
            }
        } endfor
    } endfor

    // Count degrees for the remaining nodes
    forall_nodes(G, node) {
        if (buckets->contains(node)) {
            forall_out_edges(G, edge, node) {
                NodeID target = G.getEdgeTarget(edge);
                if (buckets->contains(target)) {
                    buckets->increment(node);
                }
            } endfor
        }
    } endfor
    
    // The bucket queue now only contains the nodes of S that have 
    // no neighbor in the solution.
    while(1) {
    // Pick the one with the least neighbors within S
        int start_node = buckets->pickSmallest();
        if (start_node == -1) break;
        buckets->remove(start_node);
        // Insert it in the solution 
        G.setPartitionIndex(start_node, 1);
        // Add to candidates
        if (config.optimize_candidates) candidates.push_back(start_node);
        forall_out_edges(G, edge, start_node) {
            NodeID target = G.getEdgeTarget(edge);
            // Remove any remaining neighbors
            if (buckets->contains(target)) {
                buckets->remove(target);
                if (config.optimize_candidates && G.getPartitionIndex(target) == 1) candidates.push_back(target);
                G.setPartitionIndex(target, 0);
                // Decrease the next neighbors 
                forall_out_edges(G, target_edge, target) {
                    NodeID target_neighbor = G.getEdgeTarget(target_edge);
                    if (buckets->contains(target_neighbor)) buckets->decrement(target_neighbor);
                } endfor
            }
        } endfor
    }
    delete buckets;
}

void separator_combine::apply_separator_kahip(MISConfig & config, graph_access & G, separator_pool *pool) {
    separator sep;
    pool->get_random_separator(sep);
    pool->apply_separator(config, G, sep);
}

