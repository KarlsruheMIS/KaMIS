/*******
 * reduction_evolution.cpp
 * Purpose: Recursive evolutionary framework including reduction techniques.
 *
 *****************************************************************************/

#include "reduction_evolution.h"

#include <cstddef>
#include <limits.h>

#include "graph_access.h"
#include "mmwis_log.h"
#include "diversifier.h"
#include "cover_combine.h"
#include "multiway_combine.h"
#include "reductions.h"
#include "separator_combine.h"
#include "graph_io.h"
#include "struction_branch_and_reduce_algorithm.h"
#include "struction/app/configuration_struction.h"

using namespace mmwis;

template <typename reducer>
reduction_evolution<reducer>::reduction_evolution() {
    pool = new separator_pool();
    pool_counter = 0;
    reduction_counter = 0;
}

template <typename reducer>
reduction_evolution<reducer>::~reduction_evolution() {
    delete pool;
    pool = nullptr; 
}

template <typename reducer>
void reduction_evolution<reducer>::init(MISConfig & mis_config, graph_access & G) {
    // Set the RNG
    srand(mis_config.seed);
    random_functions::setSeed(mis_config.seed);
    // Start the timer
    mmwis_log::instance()->restart_evo_timer();
}


template <typename reducer>
void reduction_evolution<reducer>::extract_nodes(MISConfig & mis_config, 
                              graph_access & G, 
                              std::vector<NodeID> & is_nodes,
                              std::vector<NodeID> & other_nodes,
                              std::unique_ptr<reducer> &full_reducer) {
    std::vector<NodeID> nodes;
    std::vector<NodeID> nodes_excluded;
    island.get_best_individual_nodes(mis_config, G, nodes, nodes_excluded);
    full_reducer->force_into_independent_set(nodes);
    if (nodes_excluded.size()> 0) full_reducer->exclude_nodes(nodes_excluded);
}

template <typename reducer>
void reduction_evolution<reducer>::get_extract_nodes(MISConfig & mis_config, 
                              graph_access & G, 
                              std::vector<NodeID> & nodes,
                              std::vector<NodeID> & nodes_excluded) {
    island.get_best_individual_nodes(mis_config, G, nodes, nodes_excluded);
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
                               std::vector<NodeID> & worse_nodes,
                               bool recursive,
                               std::unique_ptr<reducer> &full_reducer) {
    // extract mis nodes
    if (best_nodes.size() > 0) {
        full_reducer->force_into_independent_set(best_nodes);
    }
    if (worse_nodes.size()> 0) {
        full_reducer->exclude_nodes(worse_nodes);
    }

    // perform reduction
    full_reducer->reduce_graph();

    // retrieve reduced graph
	full_reducer->build_graph_access(reduced, reverse_mapping);

	NodeWeight is_weight = full_reducer->get_current_is_weight();
    int remaining_size = reduced.number_of_nodes();

    return remaining_size;
}

template <typename reducer>
void reduction_evolution<reducer>::set_local_iterations(MISConfig & mis_config, graph_access & G) {
    mis_config.ils_iterations = std::min(G.number_of_nodes(), mis_config.ils_iterations);
}


template <typename reducer>
NodeWeight reduction_evolution<reducer>::perform_mis_search(MISConfig & mis_config, 
                                   graph_access & G, 
                                   std::vector<bool> & independent_set,
                                   std::vector<NodeID> & best_nodes,
                                   std::vector<NodeID> & worse_nodes,
                                   bool& solved_exactly, 
                                   bool recursive = false, 
                                   NodeWeight weight_offset) {

    is_base = weight_offset;

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
    full_reducer = std::unique_ptr<reducer>(new reducer(G, mis_config));

    int remaining_size = reduce(mis_config, G, reduced, is_nodes, other_nodes, reverse_mapping, best_nodes, worse_nodes, recursive, full_reducer);
    ASSERT_TRUE(check_IS_partition(reduced));
    reduced_exact_ind_calculated=false;

    is_base += full_reducer->get_current_is_weight();
    mmwis_log::instance()->print_reduction(mis_config, is_base, reduced.number_of_nodes());
    if(reduced.number_of_nodes() == 0) {
        mmwis_log::instance()->set_operator("Fully Reduced");
        mmwis_log::instance()->set_result_operator(is_base);
        mmwis_log::instance()->set_best_weight(mis_config, is_base);
     // solve exactly if |V|<mis_config.V_solve_exact:
    } else if(reduced.number_of_nodes() < mis_config.V_solve_exact) {
        //set config for cyclicFast:
        MISConfig struction_misconfig;
        configuration_struction struction_config;
        struction_config.cyclicFast(struction_misconfig);
        
        double struction_time = mis_config.time_limit - mmwis_log::instance()->get_total_timer();

        // if struction time large, only use time_solve_exact and repeat
        // reduce V_exact threshold 
        if (mis_config.time_solve_exact <= struction_time ) {
            std::cout << "solve reduced graph (remaining number of nodes: " << reduced.number_of_nodes() <<" using: " <<mis_config.time_solve_exact<< " seconds from remaining time: " << struction_time  << ")" << std::endl;
            struction_time = mis_config.time_solve_exact;
            mis_config.V_solve_exact *=0.5;
        } else {
            std::cout << "solve reduced graph (remaining number of nodes: " << reduced.number_of_nodes() <<" using all remaining time: " << struction_time  << ")" << std::endl;
        }
        
        struction_misconfig.time_limit = struction_time;
        
        struction::cout_handler::disable_cout();
        struction::branch_and_reduce_algorithm exact_solver(reduced, struction_misconfig);
        bool timeout = !exact_solver.run_branch_reduce();
	    exact_solver.apply_branch_reduce_solution(reduced);
        NodeWeight solution_weight_check = exact_solver.get_current_is_weight() + is_base;
        struction::cout_handler::enable_cout();

        reduced_exact_ind_calculated=true;
        
        NodeID *solution = new NodeID[reduced.number_of_nodes()];
        NodeWeight solution_weight = 0;
        
        forall_nodes(reduced, node) {
            if (reduced.getPartitionIndex(node) == 1) {
                solution[node] = 1;
                solution_weight += reduced.getNodeWeight(node);
            } 
            else solution[node] = 0;
        } endfor

        red_exact_ind.solution = solution;
        red_exact_ind.solution_weight = solution_weight;

        mmwis_log::instance()->set_operator("Exact Solver On Reduced Graph");
        mmwis_log::instance()->set_result_operator(red_exact_ind.solution_weight);
        mmwis_log::instance()->set_best_weight(mis_config, red_exact_ind.solution_weight+is_base);


        if (timeout) {
            std::cout << "solved reduced graph to: " << solution_weight << std::endl;
        } else {
            std::cout << "solved reduced graph exactly: " << solution_weight << std::endl;
            solved_exactly = true;
        } 

        if (solved_exactly) { 
            add_exact_solution(mis_config, reduced, independent_set, red_exact_ind, reverse_mapping);
            add_reductions(mis_config, G, reduced, reverse_mapping, independent_set, full_reducer);
            build_final_solution(mis_config, G, independent_set, weight_offset);
            return is_base;
         }
    }

    // Reconfigure local search iterations
    if (recursive) set_local_iterations(mis_config, reduced);

    // Stop if time limit or the graph is completely reduced
    if (mis_config.time_limit <= mmwis_log::instance()->get_total_timer() || (remaining_size == 0 && recursive)) {
        add_reductions(mis_config, G, reduced, reverse_mapping, independent_set, full_reducer);
        build_final_solution(mis_config, G, independent_set, weight_offset);
        return is_base;

    } else if (remaining_size == 0 && !recursive) {
        add_reductions(mis_config, G, reduced, reverse_mapping, independent_set, full_reducer);
        build_final_solution(mis_config, G, independent_set, weight_offset);

        return is_base;
    }
    
 
    std::vector<bool> coarser_is(reduced.number_of_nodes(), false);
    coarser_max = 0;
 
    perform_evolutionary(mis_config, reduced, coarser_is);
    perform_reverse_mapping(mis_config, reduced, independent_set, coarser_is, reverse_mapping);
    add_reductions(mis_config, G, reduced, reverse_mapping, independent_set, full_reducer);

    if (!recursive) {
        build_final_solution(mis_config, G, independent_set);
    }
 
    return std::max(finer_max, coarser_max);
}

template <typename reducer>
void reduction_evolution<reducer>::perform_evolutionary(MISConfig & mis_config,
                                               graph_access & G, 
                                               std::vector<bool> & independent_set) {


    // Initialize the population
    island.reset(mis_config, G);
    fill_population(mis_config, G);

    if(mis_config.time_limit <= mmwis_log::instance()->get_total_timer())
        return;

    // Build a separator pool 
    init_pool(mis_config, G);
    calculate_population_scores(mis_config, G);

    mmwis_log::instance()->restart_evo_timer();

    bool solved_exactly = false;
    do {
        if (mis_config.time_limit <= mmwis_log::instance()->get_total_timer()) break;
        // Increment and print the current round
        mmwis_log::instance()->inc_rounds();
        mmwis_log::instance()->print_round(mis_config);

        if (pool_counter >= mis_config.pool_threshold) {
            pool->renew_pool(mis_config, G, false, true, island);
            pool_counter = 0;
        }

        if (mis_config.time_limit <= mmwis_log::instance()->get_total_timer()) break;
        // Perform the evolutionary algorithm
        perform_local_mis(mis_config, G);

        // If converged or time limit is reached, recurse
        if (reduction_counter >= mis_config.reduction_threshold || mmwis_log::instance()->get_evo_timer() > mis_config.evo_time_limit) {

                // Stop if recursive call is done
                std::vector<NodeID> best_nodes;
                std::vector<NodeID> worse_nodes;
                if (mis_config.extract_best_nodes) get_extract_nodes(mis_config, G, best_nodes, worse_nodes);
                reduction_evolution re;
                mmwis_log::instance()->restart_evo_timer();
                coarser_max = re.perform_mis_search(mis_config, G, independent_set, best_nodes, worse_nodes, solved_exactly, true, is_base);
                break;
        }

    } while (mmwis_log::instance()->get_total_timer() <= mis_config.time_limit && !solved_exactly);

    delete pool;
    pool = nullptr;
}

template <typename reducer>
void reduction_evolution<reducer>::add_reductions(MISConfig & mis_config, 
                               graph_access & G,
                               graph_access & reduced,
                               std::vector<NodeID> & reverse_mapping,
                               std::vector<bool> & independent_set, 
                               std::unique_ptr<reducer> & full_reducer) {

    full_reducer->set_node_status(independent_set, G, reduced, reverse_mapping);
    full_reducer->reverse_reduction(G, reduced, reverse_mapping);
    full_reducer->update_independent_set(independent_set);

}

template <typename reducer>
void reduction_evolution<reducer>::add_exact_solution(MISConfig & mis_config, 
                                                  graph_access & reduced,
                                                  std::vector<bool> & independent_set,
                                                  individuum_mis& exact_individuum,
                                                  std::vector<NodeID> & reverse_mapping) {
        forall_nodes(reduced, node) {
            independent_set[reverse_mapping[node]] = exact_individuum.solution[node];
        } endfor
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
                                               std::vector<bool> & independent_set,
                                               NodeWeight weight_offset) {
 
    // Create individuum for final independent set
    ASSERT_TRUE(check_IS_vector(independent_set, G));
    forall_nodes(G, node) {
        G.setPartitionIndex(node, independent_set[node]);
    } endfor


    // Apply HILS
    /* hils iterate(mis_config); */
    /* iterate.perform_ils(G, mis_config.ils_iterations, weight_offset); */

    individuum_mis final_mis;
    NodeID *solution = new NodeID[G.number_of_nodes()];
    final_mis.solution_weight = island.create_solution(G, solution);
    final_mis.solution = solution;
    island.set_mis_for_individuum(mis_config, G, final_mis);

    forall_nodes(G, node) {
        independent_set[node] = final_mis.solution[node];
    } endfor

    mmwis_log::instance()->set_operator("Combine reduction");
    /* mmwis_log::instance()->reset_best_size(); */
    mmwis_log::instance()->set_best_weight(mis_config, final_mis.solution_weight);

    delete[] solution;
    solution = NULL;
} 

template <typename reducer>
NodeWeight reduction_evolution<reducer>::collect_best_mis(MISConfig & mis_config, graph_access & G, individuum_mis & out) {
    island.get_best_individuum(out);
    return out.solution_weight;
}

template <typename reducer>
void reduction_evolution<reducer>::fill_population(MISConfig & mis_config, graph_access & G) {
    // If the population is not filled simply create new individuals
    double remaining_time = mis_config.time_limit - mmwis_log::instance()->get_total_timer();
    if (reduced_exact_ind_calculated) {
        island.insert(mis_config, G, red_exact_ind, remaining_time);
    }

    mmwis_log::instance()->print_init_title();
    bool found_optimal_individuum = false;
    while (!island.is_full() && mmwis_log::instance()->get_total_timer() <= mis_config.time_limit) {
        // Diversify?
        if (mis_config.diversify) {
            diversifier div;
            div.diversify(mis_config);
        }

        individuum_mis ind;
        double remaining_time = mis_config.time_limit - mmwis_log::instance()->get_total_timer();
        found_optimal_individuum = island.create_individuum(mis_config, G, ind, remaining_time); 

        mmwis_log::instance()->set_operator("Initial");
        mmwis_log::instance()->set_result_operator(ind.solution_weight);

        individuum_mis best;
        double time_limit = mis_config.time_limit - mmwis_log::instance()->get_total_timer();
        island.insert(mis_config, G, ind, time_limit);

        // Set average and best solution and log information
        unsigned int best_after = collect_best_mis(mis_config, G, best);
        mmwis_log::instance()->set_best_weight(mis_config, best_after + is_base);
        if (found_optimal_individuum) {
            mis_config.time_limit = 0;
            break;
        }
    } 
}

template <typename reducer>
void reduction_evolution<reducer>::calculate_population_scores(MISConfig & mis_config, graph_access & G) {
    if (mmwis_log::instance()->get_total_timer() > mis_config.time_limit) return;

    if (mis_config.use_multiway_vc) pool->calculate_partition_scores(mis_config, G, island);
    else pool->calculate_separator_scores(mis_config, G, island);
}

template <typename reducer>
void reduction_evolution<reducer>::perform_local_mis(MISConfig & mis_config, graph_access & G) {
    unsigned int repetitions = mis_config.repetitions;
    // Main evolutionary algorithm part
    for (unsigned int i = 0; i < repetitions; ++i) {
        // Diversify?
        if (mis_config.diversify) {
            diversifier div;
            div.diversify(mis_config);
        }
        mmwis_log::instance()->inc_repetitions();

        individuum_mis first;
        individuum_mis second;
        individuum_mis out; 
        individuum_mis out_second; 
        mmwis_log::instance()->restart_operator_timer();

        int combine = random_functions::nextInt(0, 2);
        // Reproduction
        if (combine < 2) {
            if (mis_config.enable_tournament_selection) island.get_two_individuals_tournament(mis_config, first, second);
            else island.get_two_random_individuals(first, second);
        }
        // Crossover
        // Node Separator
        if (combine == 0) {
            mmwis_log::instance()->set_operator("Node separator");
            separator_combine combinator(mis_config);
            combinator.combine(mis_config, G, pool, first, second, out, out_second); 
        }
        // Vertex cover
        else if (combine == 1) {
            mmwis_log::instance()->set_operator("Vertex cover");
            cover_combine combinator(mis_config);
            combinator.combine(mis_config, G, pool, first, second, out, out_second);
        } 
        // Multiway
        else if (combine == 2) {
            mmwis_log::instance()->set_operator("Multiway");
            multiway_combine combinator(mis_config);
            combinator.combine(mis_config, G, pool, island, out);
        }
        
        // Mutation
        int decision = random_functions::nextInt(0, 9);
        if (decision < mis_config.flip_coin) {
            island.mutate(mis_config, G, out, mis_config.ils_time_limit);
            if (combine < 2) island.mutate(mis_config, G, out_second, mis_config.ils_time_limit);
        }

        individuum_mis best;
        // Only select better offspring
        individuum_mis better_offspring = out;
        if (combine < 2) {
            if (out.solution_weight < out_second.solution_weight) {
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
        double remaining_time = mis_config.time_limit - mmwis_log::instance()->get_total_timer();
        bool success = island.insert(mis_config, G, better_offspring, remaining_time);
        // Update the separator cache
        if (success) pool->update_scores_for_individuum(mis_config, G, better_offspring);
        // Only log the bigger offspring for the node separator or vertex cover
        mmwis_log::instance()->set_result_operator(better_offspring.solution_weight);

        // New best?
        unsigned int best_after = collect_best_mis(mis_config, G, best);

        // Set average and best solution and log information
        mmwis_log::instance()->set_best_weight(mis_config, best_after + is_base);
        if (best_after > best_before) {
            pool_counter = 0;
            reduction_counter = 0;
        }
        else {
            pool_counter++;
            reduction_counter++;
        }

        // Stop if time limit was reached
        if (mmwis_log::instance()->get_evo_timer() > mis_config.evo_time_limit) break;
        if (mmwis_log::instance()->get_total_timer() > mis_config.time_limit) break;
    }
    ASSERT_TRUE(check_IS_partition(G));
}

template <typename reducer>
bool reduction_evolution<reducer>::check_IS_vector(std::vector<bool> & independent_set, graph_access& G) {
    for (size_t node =0; node < independent_set.size(); node++) {
        if (independent_set[node]) {
            forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if (independent_set[target]) {
                    std::cout << "ERROR: Vector not an independent set!" << std::endl;
                    return false;
                }
            }endfor
        }
    }
    return true;
}

template <typename reducer>
bool reduction_evolution<reducer>::check_IS_partition(graph_access & G) {
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node)) {
            forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if (G.getPartitionIndex(target)) {
                    std::cout << "ERROR: Partition not an independent set!" << std::endl;
                    return false;
                }
            } endfor
        }
    } endfor
    return true;
}


template class reduction_evolution<branch_and_reduce_algorithm>;
