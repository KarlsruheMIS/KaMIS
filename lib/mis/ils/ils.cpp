/**
 * ils.cpp
 * Purpose: Perform the iterated local search (ILS) as described by Andrade et al.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "ils.h"

#include "greedy_mis.h"
#include "diversifier.h"
#include "data_structure/operation_log.h"

ils::ils() {
    // Set config parameters
    plateau_down = 1;
    plateau_up = 0;
    plateau_best = 0;
    plateau_best_again = 0;
    pden_floor = 1;
    delta_penalty = 0;
    limit_plateau = false;
    swap_on_failure = true;

    force_list = NULL;
    best_solution.solution = NULL;
}

ils::~ils() {
    reset();
}

void ils::perform_ils(MISConfig & config, graph_access & G, unsigned int iteration_limit) {
    reset();
    // Init operation log
    operation_log::instance()->init(G.number_of_nodes());
    operation_log::instance()->reset();
    operation_log::instance()->deactivate();
            
    // Init last forced
    last_forced.clear();
    last_forced.resize(G.number_of_nodes(), 0);

    // Init RNG
    srand(config.seed);
    random_functions::setSeed(config.seed);

    t.restart();

    local.preprocess_graph(G);
    local.direct_improvement(G);

    // Add the solution to the "population"
    pop.init(config, G);
    NodeID *solution = new NodeID[G.number_of_nodes()];
    best_solution.solution_size = pop.create_solution(G, solution);
    best_solution.solution = solution;

    // Initialize "friends".
    cand = &local.candidates;
    perm = &local.perm;
    one = &local.perm.onetight_all;

    unsigned int denominator = 4; 
    unsigned int plateau = plateau_best * perm->get_solution_size();
    force_list = new candidate_list();
    force_list->init(G.number_of_nodes());

    unsigned int iterations = 0;
    while (iterations < iteration_limit) {
        iterations++;
        // Stop if the time limit was passed
        if (t.elapsed() > config.time_limit) break;
        // if (iterations % 10000 == 0) std::cout << perm->get_solution_size() << " [" << t.elapsed() << "]" << std::endl;

        plateau--;
        if (plateau < 0) plateau = 0;

        // Activate the operation log
        operation_log::instance()->reset();
        operation_log::instance()->activate();

        int non_solution = G.number_of_nodes() - perm->get_solution_size();
        if (non_solution <= 0) break;

        denominator = 2*perm->get_solution_size();
        unsigned int solution_before = perm->get_solution_size();
        unsigned int forced = 1;

        // Check if more than one vertex should be forced and if so
        // determine how many
        unsigned int outer_denominator = denominator;
        if (random_functions::nextInt(1, outer_denominator) == 1) {
            unsigned int inner_denominator = 2;
            forced++;
            while (random_functions::nextInt(1, inner_denominator) == 1) forced++;
        }

        // Look at a certain number of candidates for the forceful insertion
        // and pick the one that hasn't been in the solution the longest time.
        int v;
        int best_v = -1;
        non_solution = G.number_of_nodes() - perm->get_solution_size();
        for (unsigned int i = 0; i < config.force_cand; ++i) {
            unsigned int position = random_functions::nextInt(0, non_solution - 1);
            v = perm->get_non_solution_node(position);
            if (best_v == -1) best_v = v;
            else if (last_forced[v] < last_forced[best_v]) best_v = v;
        }
        v = best_v;

        // Actual insertion
        if (forced == 1) {
            force(config, G, v, NULL);
        }
        else {
            force_list->reset();
            for (unsigned int i = 0; i < forced; ++i) {
                // No vertex left then stop
                if (v == -1) {
                    forced = i;
                    break;
                }
                // Force the node
                force(config, G, v, force_list);
                v = -1;
                // Pick a node thats close to the ones we already removed
                unsigned int num_cand = force_list->get_size();
                if (num_cand == 0) break;
                force_list->random_permute();
                unsigned int valid_count = 0;
                for (unsigned int i = 0; i < num_cand; ++i) {
                    NodeID w = force_list->pick(i);
                    forall_out_edges(G, edge, w) {
                        NodeID u = G.getEdgeTarget(edge);
                        if (perm->is_solution_node(u)) continue;
                        if (force_list->contains(u)) continue;

                        valid_count++;
                        if (random_functions::nextInt(1, valid_count) == 1) v = u;
                    } endfor
                }
            }
        }

        if (forced == 0) printf("Should have at least inserted one vertex\n");
        if (forced != 1) v = -1;

        // See if the new solution can be further improved
        local.make_maximal(G);
        if (v >= 0) local.direct_improvement(G, true, v);
        else local.direct_improvement(G);

        // individuum_mis temp;
        // NodeID *temp_sol = new NodeID[G.number_of_nodes()];
        // temp.solution_size = pop.create_solution(G, temp_sol);
        // temp.solution = temp_sol;
        // ASSERT_TRUE(pop.is_mis(config, G, temp));

        // Check the difference
        unsigned int solution_after = perm->get_solution_size();
        int delta = solution_before - solution_after;
        int delta_best = best_solution.solution_size - solution_after;

        // New best?
        // printf("Insertion\n");
        // if (!pop.is_mis(config, G, best_solution)) printf("No MIS in %d\n", i);
        // else printf("Is MIS in %d\n", i);
        if (solution_after > best_solution.solution_size) {
            forall_nodes(G, node) {
                best_solution.solution[node] = G.getPartitionIndex(node);
            } endfor
            best_solution.solution_size = perm->get_solution_size();
            plateau = plateau_best * perm->get_solution_size();
        } 
        else if (solution_after > solution_before) {
            if (solution_after == best_solution.solution_size) {
                plateau = plateau_best_again * perm->get_solution_size();
            }
            else {
                plateau = plateau_up * perm->get_solution_size();
            }
        }
        // ASSERT_TRUE(pop.is_mis(config, G, best_solution));
    
        // printf("After insertion\n");
        // if (!pop.is_mis(config, G, best_solution)) printf("No MIS in %d\n", i);
        // else printf("Is MIS in %d\n", i);

        // If the new solution solution isn't better 
        // it's possible to stay put
        if (delta > 0 || (delta == 0 && limit_plateau)) {
            int c = delta_penalty;
            int pden = pden_floor + (delta_best + c) * (delta + c);

            // Keep the solution with a small probability
            if (plateau > 0 || random_functions::nextInt(1, pden) != 1) {
                // Undo operations that led to the new solution
                unwind(G);
                if (swap_on_failure) {
                    // Perform a random 1-swap
                    int x = -1;
                    if (!one->is_empty()) {
                        for (unsigned int i = 0; i < config.force_cand; ++i) {
                            int y = one->pick_random();
                            if (x == -1 || last_forced[y] < last_forced[x]) x = y;
                        }
                        force(config, G, x, NULL);
                        // Try to improve the solution
                        local.direct_improvement(G, true, x);
                        // Improvement found?
                        if (perm->get_solution_size() > best_solution.solution_size) {
                            forall_nodes(G, node) {
                                best_solution.solution[node] = G.getPartitionIndex(node);
                            } endfor
                            best_solution.solution_size = perm->get_solution_size();
                            plateau = plateau_best * perm->get_solution_size();
                        }
                    }
                }
            } else {
                plateau = plateau_down * perm->get_solution_size();
            }
        }

        // Update the "Forced list"
        for (unsigned int pos = operation_log::instance()->get_size(); pos > 0; pos--) {
            int v = operation_log::instance()->peek(pos - 1);
            if (v < 0) v = -v;
            if (!perm->is_solution_node(v)) last_forced[v] = iterations;
        }

        // Candidates should be empty anyway
        if (cand && !cand->is_empty()) {
            printf("Candidates should be empty\n");
            cand->reset();
        }

        ASSERT_TRUE(pop.is_mis(config, G, best_solution));
        ASSERT_TRUE(perm->check_consistency(G));
    }

    forall_nodes(G, node) {
        G.setPartitionIndex(node, best_solution.solution[node]);
    } endfor
}

void ils::force(MISConfig & config, graph_access & G, NodeID v, candidate_list *force_list) {
    // printf("Force %d\n", v);
    if (perm->is_solution_node(v)) printf("Vertex already in the solution.\n");
    if (force_list) {
        if (!force_list->contains(v)) force_list->insert(v);
        else printf("Vertex already in the forcelist.\n");
    }

    forall_out_edges(G, edge, v) {
        NodeID w = G.getEdgeTarget(edge);
        if (force_list) {
            if (!force_list->contains(w)) force_list->insert(w);
        }
        if (perm->is_solution_node(w)) {
            // printf("Remove %d\n", w);
            perm->remove_from_solution(w, G);
            if (cand) local.update_candidates(w, G);
            // operation_log::instance()->report_remove(w); 
        }
    } endfor

    perm->add_to_solution(v, G);
    // operation_log::instance()->report_insert(v); 
    if (cand && !cand->contains(v)) cand->insert(v);
}

void ils::unwind(graph_access & G) {
    // printf("Unwind\n");
    operation_log::instance()->deactivate(); 
    // operation_log::instance()->print();
    while (!operation_log::instance()->is_empty()) {
        int v = operation_log::instance()->unwind();
        if (v < 0) {
            // printf("Add %d\n", w);
            perm->add_to_solution(-v, G);
        }
        else {
            // printf("Remove %d\n", v);
            perm->remove_from_solution(v, G);
        }
    }
    operation_log::instance()->activate();
}

void ils::reset() {
    if (force_list != NULL) { 
        delete force_list;
        force_list = NULL;
    }
    if (best_solution.solution != NULL) {
        delete [] best_solution.solution;
        best_solution.solution = NULL;
    }
}

