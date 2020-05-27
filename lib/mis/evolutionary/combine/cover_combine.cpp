/**
 * cover_combine.cpp
 * Purpose: Combine two individuals with a vertex cover based mechanism.
 *
 *****************************************************************************/

#include "cover_combine.h"

#include "graph_extractor.h"
#include "random_functions.h"
#include "hopcroft/bipartite_vertex_cover.h"
#include "data_structure/priority_queues/bucket_array.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

cover_combine::cover_combine() {
    
}

cover_combine::~cover_combine() {
    
}

void cover_combine::combine(MISConfig & config, graph_access & G, separator_pool *pool, individuum_mis & first, individuum_mis & second, individuum_mis & out_first, individuum_mis & out_second) {
    // Create the input covers
    population_mis pop;
    individuum_mis first_cover, second_cover;
    create_complement(G, first, first_cover);
    create_complement(G, second, second_cover);

    ASSERT_TRUE(pop.is_vertex_cover(config, G, first_cover));
    ASSERT_TRUE(pop.is_vertex_cover(config, G, second_cover));

    apply_partition_kahip(config, G, pool);

    // Create a second graph to store the vertex cover
    graph_access cover_graph;
    G.copy(cover_graph);
    cover_graph.resizeSecondPartitionIndex(cover_graph.number_of_nodes());

    // Build crossovers
    // Use the cover graph to store the vertex cover
    // First partition index contains A + B
    // Second partition index contains B + A
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 0) {
            cover_graph.setPartitionIndex(node, first_cover.solution[node]);
            cover_graph.setSecondPartitionIndex(node, second_cover.solution[node]);
        }
        else {
            cover_graph.setPartitionIndex(node, second_cover.solution[node]);
            cover_graph.setSecondPartitionIndex(node, first_cover.solution[node]);
        }
    } endfor

    // Fix the partition boundary
    if (config.use_hopcroft) {
        build_hopcroft_cover(config, G, cover_graph);
    } else {
        build_vertex_cover_candidates(config, cover_graph, 0);
        build_vertex_cover_candidates(config, cover_graph, 1);
    }

    // Create individuals for the final vertex covers
    NodeID *first_vertex_cover = new NodeID[cover_graph.number_of_nodes()];
    NodeID *second_vertex_cover = new NodeID[cover_graph.number_of_nodes()];
    unsigned int first_vertex_cover_size = 0;
    unsigned int second_vertex_cover_size = 0;

    forall_nodes(cover_graph, node) {
        if (cover_graph.getPartitionIndex(node) == 1) {
            first_vertex_cover[node] = 1;
            first_vertex_cover_size++;
        }
        else first_vertex_cover[node] = 0;
        if (cover_graph.getSecondPartitionIndex(node) == 1) {
            second_vertex_cover[node] = 1;
            second_vertex_cover_size++;
        }
        else second_vertex_cover[node] = 0;
    } endfor

    individuum_mis first_out_complement;
    first_out_complement.solution = first_vertex_cover;
    first_out_complement.solution_size = first_vertex_cover_size;
    
    ASSERT_TRUE(pop.is_vertex_cover(config, cover_graph, first_out_complement));
    
    individuum_mis second_out_complement;
    second_out_complement.solution = second_vertex_cover;
    second_out_complement.solution_size = second_vertex_cover_size;
    
    ASSERT_TRUE(pop.is_vertex_cover(config, cover_graph, second_out_complement));

    // Create mis 
    vertex_cover_to_mis(config, cover_graph, first_out_complement, out_first);

    // Repeat the procedure for the second offspring
    forall_nodes(cover_graph, node) {
        cover_graph.setPartitionIndex(node, cover_graph.getSecondPartitionIndex(node));
    } endfor
    vertex_cover_to_mis(config, cover_graph, second_out_complement, out_second);

    delete [] first_cover.solution;
    first_cover.solution = NULL;
    delete [] second_cover.solution;
    second_cover.solution = NULL;
    delete [] first_vertex_cover;
    first_vertex_cover = NULL;
    delete [] second_vertex_cover;
    second_vertex_cover = NULL;
}

void cover_combine::create_complement(graph_access & G, individuum_mis & in, individuum_mis & out) {
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

void cover_combine::build_hopcroft_cover(MISConfig & config, graph_access & G, graph_access & cover) {
    // Improve the vertex cover by finding a minimal vertex cover for the boundary
    // Create and extract the boundary
    std::vector<std::pair<std::vector<NodeID>, std::vector<NodeID>>> boundaries;
    extract_boundaries(G, boundaries);
    // Calculate minimal vertex cover for the boundary
    bipartite_vertex_cover bipartite;
    for (unsigned int i = 0; i < boundaries.size(); ++i) {
        std::vector<NodeID> vertex_cover;
        std::vector<NodeID> lhs = boundaries[i].first;
        std::vector<NodeID> rhs = boundaries[i].second;
        bipartite.hopcroft_cover(G, lhs, rhs, vertex_cover);
        // Add nodes from the bipartite cover
        for (unsigned int i = 0; i < vertex_cover.size(); ++i) {
            NodeID node = vertex_cover[i];
            cover.setPartitionIndex(node, 1);
            cover.setSecondPartitionIndex(node, 1);
        }
    }
}

void cover_combine::build_vertex_cover_candidates(MISConfig & config, graph_access & G, unsigned int partition_index) {
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
            if (partition_index == 0) {
                if (G.getPartitionIndex(node) == 0 && G.getPartitionIndex(target) == 0) {
                    covered = false;
                    buckets->decrement(node);
                }
            } else {
                if (G.getSecondPartitionIndex(node) == 0 && G.getSecondPartitionIndex(target) == 0) {
                    covered = false;
                    buckets->decrement(node);
                }
            }
        } endfor
        // If it is covered completely, remove it
        if (covered) {
            buckets->remove(node);
            uncovered_nodes--;
        }
    } endfor
    
    unsigned int nodes_added = 0;
    while(1) {
        // Pick the node with the most uncovered edges
        int start_node = buckets->pickSmallest();
        if (start_node == -1) break;
        buckets->remove(start_node);
        // Insert it into the solution
        if (partition_index == 0) G.setPartitionIndex(start_node, 1);
        else G.setSecondPartitionIndex(start_node, 1);
        nodes_added++;
        
        forall_out_edges(G, edge, start_node) {
            NodeID target = G.getEdgeTarget(edge);
            // Remove any neighbors that are now completely covered
            if (buckets->contains(target)) {
                buckets->increment(target);
                bool covered = true;
                forall_out_edges(G, target_edge, target) {
                    NodeID target_neighbor = G.getEdgeTarget(target_edge);
                    if (partition_index == 0) {
                        if (buckets->contains(target_neighbor) &&
                                G.getPartitionIndex(target) == 0 &&
                                G.getPartitionIndex(target_neighbor) == 0) {
                            covered = false;
                        }
                    }
                    else {
                        if (buckets->contains(target_neighbor) &&
                                G.getSecondPartitionIndex(target) == 0 &&
                                G.getSecondPartitionIndex(target_neighbor) == 0) {
                            covered = false;
                        }
                    }
                } endfor
                if (covered) buckets->remove(target);
            }
        } endfor
    }
    // printf("Nodes added: %d, Uncovered: %d\n", nodes_added, uncovered_nodes);
    delete buckets;
}

void cover_combine::extract_boundaries(graph_access & G, std::vector<std::pair<std::vector<NodeID>, std::vector<NodeID>>> & boundaries) {
    // Extract the complete boundary
    std::vector<bool> boundary(G.number_of_nodes(), 0);
    std::vector<NodeID> candidates;
    unsigned int total_nodes = 0;
    forall_nodes(G, node) {
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            if (G.getPartitionIndex(node) != G.getPartitionIndex(target)) {
                if (!boundary[node]) {
                    boundary[node] = true;
                    candidates.push_back(node);
                    total_nodes++;
                }
                if (!boundary[target]) {
                    boundary[target] = true;
                    candidates.push_back(target);
                    total_nodes++;
                }
            }
        } endfor
    } endfor

    // Split the boundaries
    std::vector<bool> processed(G.number_of_nodes(), false);
    random_functions::permutate_vector_good(candidates, false);
    // Continue while there are unprocessed nodes left
    while(total_nodes > 0) {
        std::vector<NodeID> lhs(G.number_of_nodes(), 0);
        std::vector<NodeID> rhs(G.number_of_nodes(), 0);
        // Choose a random starting node
        // TODO: More efficient selection method
        NodeID start_node;
        do {
            start_node = candidates.back();   
            candidates.pop_back();   
        } while (processed[start_node]);

        // Create a queue with the starting node
        std::queue<NodeID> bound_queue;
        bound_queue.push(start_node);
        processed[start_node] = true;
        total_nodes--;
        if (G.getPartitionIndex(start_node) == 0) {
            lhs[start_node] = true;
            rhs[start_node] = false;
        }
        else {
            lhs[start_node] = false;
            rhs[start_node] = true;
        }
        // Grow the boundary
        while (!bound_queue.empty()) {
            NodeID current_node = bound_queue.front(); 
            bound_queue.pop();
            forall_out_edges(G, edge, current_node) {
                NodeID target = G.getEdgeTarget(edge);
                if (boundary[target] && !processed[target]) {
                    bound_queue.push(target);
                    processed[target] = true;
                    total_nodes--;
                    if (G.getPartitionIndex(current_node) == G.getPartitionIndex(target)) {
                        lhs[target] = lhs[current_node];
                        rhs[target] = rhs[current_node];
                    }
                    else {
                        lhs[target] = !lhs[current_node];
                        rhs[target] = !rhs[current_node];
                    }
                }
            } endfor
        }
        std::pair<std::vector<NodeID>, std::vector<NodeID>> bound = std::make_pair(lhs, rhs);
        unsigned int boundary_size = 0;
        forall_nodes(G, node){
            if (lhs[node] || rhs[node]) boundary_size++;
        } endfor
        boundaries.push_back(bound);
    }
}

unsigned int cover_combine::get_smaller_individuum(graph_access & G, std::vector<NodeID> mapping, individuum_mis & first, individuum_mis & second) {
    unsigned int first_count = 0;
    unsigned int second_count = 0;
    forall_nodes(G, node) {
        NodeID original_node = mapping[node];
        if (first.solution[original_node] == 1) first_count++;
        if (second.solution[original_node] == 1) second_count++;
    } endfor

    return (first_count < second_count)? 1: 2;
}

void cover_combine::vertex_cover_to_mis(MISConfig & config, graph_access & G, individuum_mis & cover, individuum_mis & ind) {
    create_complement(G, cover, ind);

    // Use local search
    local.preprocess_graph(G);
    local.make_maximal(G);
    local.direct_improvement(G);

    ind.solution_size = 0;
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 1) {
            ind.solution[node] = 1;
            ind.solution_size++;
        }
        else ind.solution[node] = 0;
    } endfor

}

void cover_combine::apply_partition_kahip(MISConfig & config, graph_access & G, separator_pool *pool) {
    partition part;
    pool->get_random_partition(part);
    pool->apply_partition(config, G, part);
}

