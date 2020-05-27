/**
 * multiway_combine.h
 * Purpose: Combine a list of individuals using k-separators/partitions.
 *
 *****************************************************************************/

#include "multiway_combine.h"

#include "data_structure/priority_queues/bucket_array.h"

multiway_combine::multiway_combine() {
    
}

multiway_combine::~multiway_combine() {
    
}

void multiway_combine::combine(MISConfig & config, graph_access & G, separator_pool *pool, population_mis & pop, individuum_mis & out) {
    std::vector<individuum_mis> candidates;
    if (config.use_multiway_vc) {
        partition part;
        apply_k_partition_kahip(config, G, pool, part);
        candidates.resize(G.get_partition_count());
        pool->get_best_candidates(config, G, pop, candidates, part);
    }
    else {
        separator sep;
        apply_k_separator_kahip(config, G, pool, sep);
        candidates.resize(G.get_partition_count());
        pool->get_best_candidates(config, G, pop, candidates, sep);
    }

    std::vector<NodeID> node_candidates;
    if (config.use_multiway_vc) build_combined_solution_vc(config, G, candidates, node_candidates);
    else build_combined_solution_ns(config, G, candidates, node_candidates);

    // Improve the solution
    if (!config.optimize_candidates) local.preprocess_graph(G);
    else local.preprocess_graph_with_candidates(G, node_candidates, node_candidates.size());
    local.make_maximal(G);
    local.direct_improvement(G);

    // Generate solution
    NodeID *solution = new NodeID[G.number_of_nodes()];
    unsigned int solution_size = pop.create_solution(G, solution);
    out.solution = solution;
    out.solution_size = solution_size;
    ASSERT_TRUE(pop.is_mis(config, G, out));
}

void multiway_combine::build_combined_solution_vc(MISConfig & config, graph_access & G, std::vector<individuum_mis> & best_parents, std::vector<NodeID> & candidates) {
    // Build combined solution
    forall_nodes(G, node) {
        unsigned int partition = best_parents[G.getPartitionIndex(node)].solution[node];

        if (partition == 0) G.setPartitionIndex(node, 0);
        else if (partition == 1) G.setPartitionIndex(node, 1);
        else printf("Error in partitioning multiway (vc).");
    } endfor

    // Build complement
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) G.setPartitionIndex(node, 1); 
        else G.setPartitionIndex(node, 0);
    } endfor

    // Fix the complement by adding uncovered nodes
    build_vertex_cover_candidates(config, G, candidates);

    // Build complement
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) G.setPartitionIndex(node, 1); 
        else G.setPartitionIndex(node, 0);
    } endfor

    // Add boundary of previous candidates
    if (config.optimize_candidates) {
        std::vector<NodeID> new_candidates;
        for (NodeID node : candidates) {
            forall_out_edges(G, edge, node) {
                NodeID target = G.getEdgeTarget(edge);
                if (G.getPartitionIndex(target) == 1) new_candidates.push_back(target);
            } endfor
        }
        candidates = new_candidates;
    }
}

void multiway_combine::build_combined_solution_ns(MISConfig & config, graph_access & G, std::vector<individuum_mis> & best_parents, std::vector<NodeID> & candidates) {
    // Build combined solution
    forall_nodes(G, node) {
        unsigned int partition = 0;
        if (G.getPartitionIndex(node) == G.get_partition_count()) partition = 2;
        else partition = best_parents[G.getPartitionIndex(node)].solution[node];

        if (partition == 0) G.setPartitionIndex(node, 0);
        else if (partition == 1) G.setPartitionIndex(node, 1);
        else if (partition == 2) G.setPartitionIndex(node, 2);
        else printf("Error in partitioning multiway (ns).");
    } endfor

    // Add separator nodes greedily
    build_separator_candidates(config, G, candidates);
}

void multiway_combine::create_complement(graph_access & G, individuum_mis & in, individuum_mis & out) {
    NodeID *complement = new NodeID[G.number_of_nodes()]; 
    int size = 0;

    forall_nodes(G, node) {
        G.setPartitionIndex(node, 1);
        if (in.solution[node] == 1) {
            G.setPartitionIndex(node, 0);
            complement[node] = 0;
        }
        if (in.solution[node] == 0) {
            G.setPartitionIndex(node, 1);
            complement[node] = 1;
            size++;
        }
    } endfor

    out.solution = complement;
    out.solution_size = size;
}

void multiway_combine::build_separator_candidates(MISConfig & config, graph_access & G, std::vector<NodeID> & candidates) {
    // Initialize a bucket queue for the separator
    bucket_array *buckets = new bucket_array(G.number_of_nodes());
    forall_nodes(G, node) {
        // Limit the bucket queue to S
        if (G.getPartitionIndex(node) != 2) buckets->remove(node);
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            // Remove nodes that have a neighboring solution node
            if (G.getPartitionIndex(target) == 1) {
                buckets->remove(node);
                G.setPartitionIndex(node, 0);
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

void multiway_combine::build_vertex_cover_candidates(MISConfig & config, graph_access & G, std::vector<NodeID> & candidates) {
    bucket_array *buckets = new bucket_array(G.number_of_nodes());
    unsigned int uncovered_nodes = G.number_of_nodes();
    // Only keep nodes that still have an uncovered edge left 
    // and are not in the solution
    forall_nodes(G, node) {
        bool covered = true;
        buckets->increment(node, G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            // Decrement the nodes value for each uncovered edge
            // A nodes value therefore is #total_edges - #uncovered_edges
            if (G.getPartitionIndex(node) == 0 && G.getPartitionIndex(target) == 0) {
                covered = false;
                buckets->decrement(node);
            }
        } endfor
        // If it is covered completely, remove it
        if (covered) {
            buckets->remove(node);
            uncovered_nodes--;
        }
    } endfor
    
    while(1) {
        // Pick the node with the most uncovered edges
        int start_node = buckets->pickSmallest();
        if (start_node == -1) break;
        buckets->remove(start_node);
        // Insert it into the solution
        G.setPartitionIndex(start_node, 1);
        if (config.optimize_candidates) candidates.push_back(start_node);
        
        forall_out_edges(G, edge, start_node) {
            NodeID target = G.getEdgeTarget(edge);
            // Remove any neighbors that are now completely covered
            if (buckets->contains(target)) {
                buckets->increment(target);
                bool covered = true;
                forall_out_edges(G, target_edge, target) {
                    NodeID target_neighbor = G.getEdgeTarget(target_edge);
                    if (buckets->contains(target_neighbor) &&
                            G.getPartitionIndex(target) == 0 &&
                            G.getPartitionIndex(target_neighbor) == 0) {
                        covered = false;
                    }
                } endfor
                if (covered) {
                    if (config.optimize_candidates && G.getPartitionIndex(target) == 1) candidates.push_back(target);
                    buckets->remove(target);
                }
            }
        } endfor
    }
    delete buckets;
}

void multiway_combine::apply_k_partition_kahip(MISConfig & config, graph_access & G, separator_pool *pool, partition & part) {
    pool->get_random_k_partition(part);
    pool->apply_partition(config, G, part);
}

void multiway_combine::apply_k_separator_kahip(MISConfig & config, graph_access & G, separator_pool *pool, separator & sep) {
    pool->get_random_k_separator(sep);
    pool->apply_separator(config, G, sep);
}

