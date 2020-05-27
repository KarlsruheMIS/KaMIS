/**
 * reduction_evolution.cpp
 * Purpose: Recursive evolutionary framework including reduction techniques.
 *
 *****************************************************************************/

#include "reduction_evolution.h"

#include <limits.h>

#include "mis_log.h"
#include "diversifier.h"
#include "cover_combine.h"
#include "multiway_combine.h"
#include "separator_combine.h"
#include "graph_io.h"

template <typename reducer>
reduction_evolution<reducer>::reduction_evolution() {
    pool = new separator_pool();
    pool_counter = 0;
    reduction_counter = 0;
}

template <typename reducer>
reduction_evolution<reducer>::~reduction_evolution() {

}

template <typename reducer>
void reduction_evolution<reducer>::init(MISConfig & mis_config, graph_access & G) {
    // Set the RNG
    srand(mis_config.seed);
    random_functions::setSeed(mis_config.seed);
    // Start the timer
    mis_log::instance()->restart_evo_timer();
}

template <typename reducer>
void reduction_evolution<reducer>::extract_nodes(MISConfig & mis_config, 
                              graph_access & G, 
                              std::vector<NodeID> & is_nodes,
                              std::vector<NodeID> & other_nodes,
                              std::unique_ptr<reducer> &full_reducer) {
    std::vector<NodeID> nodes;
    island.get_best_individual_nodes(mis_config, G, nodes);
    full_reducer->force_into_independent_set(nodes);
}

template <typename reducer>
void reduction_evolution<reducer>::get_extract_nodes(MISConfig & mis_config, 
                              graph_access & G, 
                              std::vector<NodeID> & nodes) {
    island.get_best_individual_nodes(mis_config, G, nodes);
}

template <typename reducer>
void reduction_evolution<reducer>::init_pool(MISConfig & mis_config, graph_access & G) {
    pool->init(mis_config, G);
    pool->renew_pool(mis_config, G, true, false, island);
}

template <typename reducer>
unsigned int reduction_evolution<reducer>::reduce(MISConfig & mis_config, 
                               graph_access & G, 
                               graph_access & reduced, 
                               std::vector<NodeID> & is_nodes,
                               std::vector<NodeID> & other_nodes,
                               std::vector<NodeID> & reverse_mapping,
                               std::vector<NodeID> & best_nodes,
                               bool recursive,
                               std::unique_ptr<reducer> &full_reducer) {

    // extract mis nodes
    if (best_nodes.size() > 0)
        full_reducer->force_into_independent_set(best_nodes);

    // perform reduction
    full_reducer->reduce_graph();

    // retrieve reduced graph
    full_reducer->convert_adj_lists(reduced, reverse_mapping);
    return full_reducer->number_of_nodes_remaining();
}

template <typename reducer>
void reduction_evolution<reducer>::set_local_iterations(MISConfig & mis_config, graph_access & G) {
    mis_config.ils_iterations = std::min(G.number_of_nodes(), mis_config.ils_iterations);
}

template <typename reducer>
unsigned int reduction_evolution<reducer>::perform_mis_search(MISConfig & mis_config, 
                                   graph_access & G, 
                                   std::vector<bool> & independent_set,
                                   std::vector<NodeID> & best_nodes,
                                   bool recursive, 
                                   unsigned int is_size) {
    is_base = is_size;
    // Initialize for non-recursive calls
    if (!recursive) init(mis_config, G);
   
    // Reduce the graph 
    std::vector<NodeID> is_nodes;
    std::vector<NodeID> other_nodes;
    graph_access reduced;
    std::vector<NodeID> reverse_mapping(G.number_of_nodes(), 0);
    std::unique_ptr<reducer> full_reducer;

    // Perform reductions
    // initialize full reducer
    std::vector<std::vector<int>> adj(G.number_of_nodes());

    // Build adjacency vectors
    forall_nodes(G, node) {
        adj[node].reserve(G.getNodeDegree(node));
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
            adj[node].push_back(neighbor);
        } endfor
              } endfor

                    full_reducer = std::unique_ptr<reducer>(new reducer(adj, adj.size()));
    unsigned int remaining_size = reduce(mis_config, G, reduced, is_nodes, other_nodes, reverse_mapping, best_nodes, recursive, full_reducer);
    is_base += full_reducer->get_current_is_size_with_folds();

    mis_log::instance()->print_reduction(mis_config, is_base, reduced.number_of_nodes());

    // Reconfigure local search iterations
    if (recursive) set_local_iterations(mis_config, reduced);

    // Stop if the graph is completely reduced
    if (remaining_size == 0 && recursive) {
        add_reductions(mis_config, independent_set, full_reducer);
        return is_base;
    } else if (remaining_size == 0 && !recursive) {
        add_reductions(mis_config, independent_set, full_reducer);
        build_final_solution(mis_config, G, independent_set);
        return is_base;
    };

    std::vector<bool> coarser_is(reduced.number_of_nodes(), false);
    coarser_max = 0;
    
    perform_evolutionary(mis_config, reduced, coarser_is);

    perform_reverse_mapping(mis_config, reduced, independent_set, coarser_is, reverse_mapping);
    add_reductions(mis_config, independent_set, full_reducer);

    if (!recursive) build_final_solution(mis_config, G, independent_set);

    return std::max(finer_max, coarser_max);
}

template <typename reducer>
void reduction_evolution<reducer>::perform_evolutionary(MISConfig & mis_config,
                                               graph_access & G, 
                                               std::vector<bool> & independent_set) { 
    // Build a separator pool 
    init_pool(mis_config, G);

    // Initialize the population
    island.reset(mis_config, G);

    // Fill the population
    fill_population(mis_config, G);

    do {
        // Increment and print the current round
        mis_log::instance()->inc_rounds();
        mis_log::instance()->print_round(mis_config);

        if (pool_counter >= mis_config.pool_threshold) {
            pool->renew_pool(mis_config, G, false, true, island);
            pool_counter = 0;
        }

        // Perform the evolutionary algorithm
        mis_log::instance()->print_evolution_title();
        perform_local_mis(mis_config, G);

        // If converged, recurse
        if (reduction_counter >= mis_config.reduction_threshold) {
            // Stop if recursive call is done
            std::vector<NodeID> best_nodes; 
            if (mis_config.extract_best_nodes) get_extract_nodes(mis_config, G, best_nodes);
            reduction_evolution re;
            coarser_max = re.perform_mis_search(mis_config, G, independent_set, best_nodes, true, is_base);
            break;
        }
    } while (mis_log::instance()->get_evo_timer() <= mis_config.time_limit);

    delete pool;
    pool = NULL;
}

template <typename reducer>
void reduction_evolution<reducer>::add_reductions(MISConfig & mis_config, 
                               std::vector<bool> & independent_set, 
                               std::unique_ptr<reducer> & full_reducer) {
    unsigned int size = 0;
    for (bool is : independent_set) {
        if (is) size++;
    }
    full_reducer->extend_finer_is(independent_set);
    size = 0;
    for (bool is : independent_set) {
        if (is) size++;
    }
}

template <typename reducer>
void reduction_evolution<reducer>::perform_reverse_mapping(MISConfig & mis_config, 
                                                  graph_access & coarser,
                                                  std::vector<bool> & independent_set,
                                                  std::vector<bool> & coarser_is,
                                                  std::vector<NodeID> & reverse_mapping) {
    individuum_mis best;
    finer_max = collect_best_mis(mis_config, coarser, best) + is_base;

    // If finer solution is global max
    // just perform the reverse mapping
    if (coarser_max > finer_max) {
        forall_nodes(coarser, node) {
            independent_set[reverse_mapping[node]] = coarser_is[node];
        } endfor
    } 
    // If yours is the global max
    // use it as a starting point
    else {
        individuum_mis out;
        island.get_best_individuum(out);
        forall_nodes(coarser, node) {
            independent_set[reverse_mapping[node]] = out.solution[node];
        } endfor
    }
}

template <typename reducer>
void reduction_evolution<reducer>::build_final_solution(MISConfig & mis_config, 
                                               graph_access & G, 
                                               std::vector<bool> & independent_set) {
    // Create individuum for final independent set
    individuum_mis final_mis;
    NodeID *solution = new NodeID[G.number_of_nodes()];
    forall_nodes(G, node) {
        G.setPartitionIndex(node, independent_set[node]);
    } endfor

    // Apply ILS
    ils iterate;
    iterate.perform_ils(mis_config, G, mis_config.ils_iterations);

    final_mis.solution_size = island.create_solution(G, solution);
    final_mis.solution = solution;
    island.set_mis_for_individuum(mis_config, G, final_mis);
    forall_nodes(G, node) {
        independent_set[node] = final_mis.solution[node];
    } endfor
    mis_log::instance()->set_operator("Combine reduction");
    mis_log::instance()->reset_best_size();
    mis_log::instance()->set_best_size(mis_config, final_mis.solution_size);

    delete[] solution;
    solution = NULL;
} 

template <typename reducer>
unsigned int reduction_evolution<reducer>::collect_best_mis(MISConfig & mis_config, graph_access & G, individuum_mis & out) {
    island.get_best_individuum(out);
    return out.solution_size;
}

template <typename reducer>
void reduction_evolution<reducer>::fill_population(MISConfig & mis_config, graph_access & G) {
    // If the population is not filled simply create new individuals
    mis_log::instance()->print_init_title();
    while (!island.is_full()) {
        // Diversify?
        if (mis_config.diversify) {
            diversifier div;
            div.diversify(mis_config);
        }

        // Create new individuum
        individuum_mis ind;
        island.create_individuum(mis_config, G, ind);

        mis_log::instance()->set_operator("Initial");
        mis_log::instance()->set_result_operator(ind.solution_size);

        individuum_mis best;
        island.insert(mis_config, G, ind);

        // Set average and best solution and log information
        unsigned int best_after = collect_best_mis(mis_config, G, best);
        mis_log::instance()->set_best_size(mis_config, best_after + is_base);
    } 

    if (mis_config.use_multiway_vc) pool->calculate_partition_scores(mis_config, G, island);
    else pool->calculate_separator_scores(mis_config, G, island);
}

template <typename reducer>
void reduction_evolution<reducer>::perform_local_mis(MISConfig & mis_config, graph_access & G) {
    unsigned int repetitions = mis_config.repetitions;

    timer round_timer;
    round_timer.restart();
    // Main algorithm
    for (unsigned int i = 0; i < repetitions; ++i) {
        // Diversify?
        if (mis_config.diversify) {
            diversifier div;
            div.diversify(mis_config);
        }
        mis_log::instance()->inc_repetitions();

        individuum_mis first;
        individuum_mis second;
        individuum_mis out; 
        individuum_mis out_second; 
        mis_log::instance()->restart_operator_timer();
        int combine = random_functions::nextInt(0, 2);

        // Reproduction
        if (combine < 2) {
            if (mis_config.enable_tournament_selection) island.get_two_individuals_tournament(mis_config, first, second);
            else island.get_two_random_individuals(first, second);
        }

        // Crossover
        // Node Separator
        if (combine == 0) {
            mis_log::instance()->set_operator("Node separator");
            separator_combine combinator;
            combinator.combine(mis_config, G, pool, first, second, out, out_second); 
        }
        // Vertex cover
        else if (combine == 1) {
            mis_log::instance()->set_operator("Vertex cover");
            cover_combine combinator;
            combinator.combine(mis_config, G, pool, first, second, out, out_second);
        } 
        // Multiway
        else if (combine == 2) {
            mis_log::instance()->set_operator("Multiway");
            multiway_combine combinator;
            combinator.combine(mis_config, G, pool, island, out);
        }
        
        // Mutation
        int decision = random_functions::nextInt(0, 9);
        if (decision < mis_config.flip_coin) {
            island.mutate(mis_config, G, out);
            if (combine < 2) island.mutate(mis_config, G, out_second);
        }

        individuum_mis best;
        // Only select better offspring
        individuum_mis better_offspring = out;
        if (combine < 2) {
            if (out.solution_size < out_second.solution_size) {
                better_offspring = out_second;
                delete [] out.solution;
                out.solution = NULL;
            } else {
                delete [] out_second.solution;
                out_second.solution = NULL;
            }
        }

        // Insertion
        unsigned int best_before = collect_best_mis(mis_config, G, best);
        bool success = island.insert(mis_config, G, better_offspring);
        // Update the separator cache
        if (success) pool->update_scores_for_individuum(mis_config, G, better_offspring);
        // Only log the bigger offspring for the node separator or vertex cover
        mis_log::instance()->set_result_operator(better_offspring.solution_size);

        // New best?
        unsigned int best_after = collect_best_mis(mis_config, G, best);

        // Set average and best solution and log information
        mis_log::instance()->set_best_size(mis_config, best_after + is_base);
        if (best_after > best_before) {
            pool_counter = 0;
            reduction_counter = 0;
        }
        else {
            pool_counter++;
            reduction_counter++;
        }

        // Stop if time limit was reached
        if (mis_log::instance()->get_evo_timer() > mis_config.time_limit) break;
    }
}

template class reduction_evolution<branch_and_reduce_algorithm>;
template class reduction_evolution<full_reductions>;
