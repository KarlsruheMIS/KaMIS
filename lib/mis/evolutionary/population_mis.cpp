/**
 * population_mis.cpp
 * Purpose: Represents the population used in the evolutionary framework.
 *          This includes combine, mutate and selection operators.
 *
 *****************************************************************************/

#include "population_mis.h"

#include <limits.h>

#include "ils/ils.h"
#include "greedy_mis.h"
#include "random_mis.h"
#include "greedy_vertex.h"
#include "random_functions.h"
#include "diversifier.h"
#include "macros_assertions.h"
#include "data_structure/mis_permutation.h"
#include "data_structure/priority_queues/bucket_array.h"

population_mis::population_mis() {

}

population_mis::~population_mis() {
    extinction();
}

void population_mis::init(MISConfig & config, graph_access & G) {
    internal_population.clear();
    remove_fraction = config.remove_fraction;
    insert_no_change = 0;
    population_size = config.population_size;
}

void population_mis::reset(MISConfig & config, graph_access & G) {
    extinction();
    init(config, G);
}

void population_mis::create_individuum(MISConfig & config, graph_access & G, individuum_mis & ind) {

    // Build solution
    int initial = random_functions::nextInt(0, 1);
    if (initial == 0) {
        random_mis init_solution;
        init_solution.initial_partition(config.seed, G);
    }
    else if (initial == 1) {
        greedy_mis init_solution;
        init_solution.initial_partition(config.seed, G);
    }
    // else if (initial == 2) {
    //     greedy_vertex init_solution;
    //     init_solution.initial_partition(config.seed, G);
    // }

    // Create solution for the individuum
    NodeID *solution = new NodeID[G.number_of_nodes()];
    unsigned int solution_size = create_solution(G, solution);
    ind.solution = solution;
    ind.solution_size = solution_size;
}

void population_mis::get_individuum(unsigned int id, individuum_mis & ind) {
    ind = internal_population[id]; 
}

void population_mis::get_random_individuum(individuum_mis & ind) {
    int random = random_functions::nextInt(0, population_size - 1);     
    ind = internal_population[random];
}

void population_mis::get_one_individuum_tournament(MISConfig & config, individuum_mis & ind) {

    // Create random candidates
    std::vector<unsigned int> cand(population_size);
    random_functions::permutate_vector_good(cand, true);

    unsigned int winner_id = 0;
    unsigned int winner_size = 0;
    for (unsigned int i = 0; i < config.tournament_size; ++i) {
        unsigned int ind_id = cand[i];
        if (internal_population[ind_id].solution_size > winner_size) {
            winner_id = ind_id;
            winner_size = internal_population[ind_id].solution_size;
        }
    }
    ind = internal_population[winner_id];
}

void population_mis::get_random_individuals(unsigned int k, std::vector<individuum_mis> & parents) {
    std::vector<unsigned int> pop_permutation(population_size); 
    random_functions::permutate_vector_good(pop_permutation, true);
    parents.clear();
    parents.resize(k);
    for (unsigned int i = 0; i < k; ++i) {
        parents[i] = internal_population[pop_permutation[i]];
    }
}

void population_mis::get_two_random_individuals(individuum_mis & first, individuum_mis & second) {
    int first_id = random_functions::nextInt(0, internal_population.size() - 1);
    first = internal_population[first_id];
    int second_id = 0;
    do {
        second_id = random_functions::nextInt(0, internal_population.size() - 1);
    } while (first_id == second_id);
    second = internal_population[second_id];
}

void population_mis::get_two_individuals_tournament(MISConfig & config, individuum_mis & first, individuum_mis & second) {
    get_one_individuum_tournament(config, first);
    do {
        get_one_individuum_tournament(config, second);
    } while (first.id == second.id);
}

void population_mis::get_best_individuum(individuum_mis & ind) {
    unsigned int max_solution = std::numeric_limits<unsigned int>::min();     
    unsigned int max_solution_id = 0;

    for (unsigned int i = 0; i < internal_population.size(); ++i) {
        if (internal_population[i].solution_size > max_solution) {
            max_solution = internal_population[i].solution_size;
            max_solution_id = i;
        }
    }

    ind = internal_population[max_solution_id];
}

void population_mis::set_mis_for_individuum(MISConfig & config, graph_access & G, individuum_mis & ind, bool secondary) {
    G.resizeSecondPartitionIndex(G.number_of_nodes());
    forall_nodes(G, node) {
        if (!secondary) G.setPartitionIndex(node, ind.solution[node]);
        else G.setSecondPartitionIndex(node, ind.solution[node]);
    } endfor
}

unsigned int population_mis::create_solution(graph_access & G, NodeID *solution) {
    unsigned int solution_size = 0;    
    forall_nodes(G, node) {
        if (G.getPartitionIndex(node) == 1) {
            solution[node] = 1;
            solution_size++;
        } 
        else solution[node] = 0;
    } endfor
    return solution_size;
}

void population_mis::set_population_size(unsigned int size) {
    population_size = size; 
}

void population_mis::mutate_random(MISConfig & config, graph_access & G, individuum_mis & ind) {
    get_random_individuum(ind); 
    mutate(config, G, ind);
}

void population_mis::mutate(MISConfig & config, graph_access & G, individuum_mis & ind) {
    set_mis_for_individuum(config, G, ind);
    delete [] ind.solution;
    ind.solution = NULL;

    ils iterate;
    iterate.perform_ils(config, G, config.ils_iterations);

    // Create solution for the individuum
    NodeID *solution = new NodeID[G.number_of_nodes()];
    unsigned int solution_size = create_solution(G, solution);
    ind.solution = solution;
    ind.solution_size = solution_size;
}

bool population_mis::insert(MISConfig & config, graph_access & G, individuum_mis & ind) {
    ASSERT_TRUE(is_mis(config, G, ind));
    bool successful_insertion = false;

    if (internal_population.size() < population_size) {
        ind.id = internal_population.size();
        internal_population.push_back(ind);
    } else {

        unsigned int worst_solution = std::numeric_limits<unsigned int>::max();
        unsigned int worst_id = 0;
        // Get worst solution
        for (unsigned int i = 0; i < internal_population.size(); ++i) {
            if (internal_population[i].solution_size < worst_solution) {
                worst_solution = internal_population[i].solution_size;            
                worst_id = i;
            }
        }
        individuum_mis worst = internal_population[worst_id];

        // Should the solution be inserted by force?
        if (insert_no_change > config.insert_threshold) {
            // std::cout << "Force" << std::endl;
            set_mis_for_individuum(config, G, ind);
            ils iterate;
            iterate.perform_ils(config, G, config.ils_iterations);

            ind.solution_size = 0;
            forall_nodes(G, node) {
                if (G.getPartitionIndex(node) == 1) ind.solution_size++; 
                ind.solution[node] = G.getPartitionIndex(node);
            } endfor
            replace(worst, ind);
            insert_no_change = 0;
            return true;
        }

        // Is the solution bigger than any so far?
        if (ind.solution_size <= worst_solution) {
            // std::cout << "Worse" << std::endl;
            insert_no_change++;
            delete [] ind.solution;
            ind.solution = NULL;
            return false;
        }

        // Else perform ILS and search most similar
        set_mis_for_individuum(config, G, ind);
        ils iterate;
        iterate.perform_ils(config, G, config.ils_iterations);

        ind.solution_size = 0;
        forall_nodes(G, node) {
            if (G.getPartitionIndex(node) == 1) ind.solution_size++;
            ind.solution[node] = G.getPartitionIndex(node);
        } endfor

        individuum_mis remove;
        bool valid_replacement = get_most_similar_replacement(config, G, ind, remove); 

        if (!valid_replacement) {
            // std::cout << "Not Valid" << ind.solution_size << std::endl;
            insert_no_change++;
            delete [] ind.solution;
            ind.solution = NULL;
            return false;
        } 
        else {
            // std::cout << "Valid: " << ind.solution_size << std::endl;
            replace(remove, ind);
            successful_insertion = true;
            insert_no_change = 0;
        }

    }
    ASSERT_TRUE(is_mis(config, G, ind));
    return successful_insertion;
}

bool population_mis::get_most_similar_replacement(MISConfig & config, graph_access & G, individuum_mis & ind, individuum_mis & replacement) {
    double min_ratio = 0.0;
    unsigned int max_similarity_id = 0;
    individuum_mis best;
    get_best_individuum(best);
    bool replacement_found = false;
    for (unsigned int i = 0; i < internal_population.size(); ++i) {

        if (i == best.id) continue;
        if (ind.solution_size < internal_population[i].solution_size) continue; 
        // measure similarity based on the intersection size of two solutions
        unsigned int similarities = 0;
        unsigned int nodes = 0;

        // Build intersection
        forall_nodes(G, node) {
            if (ind.solution[node] && internal_population[i].solution[node]) similarities++;
            if (ind.solution[node] || internal_population[i].solution[node]) nodes++;
        } endfor

        double similarity_ratio = (double)similarities / (double)nodes;

        if (similarity_ratio > min_ratio) {
            min_ratio = similarity_ratio;
            max_similarity_id = i;
            replacement_found = true;
        }
    }
    replacement = internal_population[max_similarity_id];
    return replacement_found;
}

void population_mis::replace(individuum_mis & remove, individuum_mis & insert) {
    insert.id = remove.id;
    delete [] internal_population[remove.id].solution;
    internal_population[remove.id].solution = NULL;
    internal_population[remove.id] = insert;
}

void population_mis::extinction() {
    for (unsigned int i = 0; i < internal_population.size(); ++i) {
        delete [] internal_population[i].solution;
        internal_population[i].solution = NULL;
    }
    internal_population.clear();
}

void population_mis::get_best_individual_nodes(MISConfig & config, graph_access & G, std::vector<NodeID> & nodes) {
    // Build PQ with nodes from best individual
    // Priority is the degree
    // Initial value is the population size
    bucket_array pq(G.number_of_nodes());
    individuum_mis best;
    get_best_individuum(best);

    NodePermutationMap permutation;
    permutation.resize(G.number_of_nodes());
    random_functions::permutate_vector_good(permutation, true);

    unsigned int is_nodes = 0;
    forall_nodes(G, n) {
        NodeID node = permutation[n];
        pq.increment(node, G.getNodeDegree(node));
        // Remove non independent set nodes
        if (best.solution[node] == 0) pq.remove(node);
        else is_nodes++;
    } endfor

    // Extract the most common nodes
    nodes.clear();
    unsigned int fraction = ceil((double) is_nodes * remove_fraction);
    unsigned int added_nodes = 0;
    while (added_nodes != fraction) {
        int cand = pq.pickSmallest();
        if (cand == -1) break;
        nodes.push_back(cand);
        pq.remove(cand);
        added_nodes++;
    }
}

void population_mis::print(MISConfig & config) {
    printf("Population:\n");
    printf("|-Size:\t\t\t%zd\n", internal_population.size());
    for (unsigned int i = 0; i < internal_population.size(); ++i) {
        printf("|-Individuum:\t\t%d\n", internal_population[i].id);
        printf("|--Solution size:\t%d\n", internal_population[i].solution_size);
    }
    printf("\n");
}

double population_mis::get_avg_solution_size() {
    double total_size = 0.0;
    for (unsigned int i = 0; i < internal_population.size(); ++i) {
        total_size += (double) internal_population[i].solution_size;
    }
    return total_size / (double) internal_population.size();
}

bool population_mis::is_full() {
    return internal_population.size() == population_size;
}

bool population_mis::is_mis(MISConfig & config, graph_access & G, individuum_mis & ind) {
    forall_nodes(G, node) {
        if (ind.solution[node] > 1) return false;
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            if (ind.solution[node] == 1 && ind.solution[target] == 1) {
                return false;
            }
        } endfor
    } endfor

    return true;
}

bool population_mis::is_vertex_cover(MISConfig & config, graph_access & G, individuum_mis & ind) {
    forall_nodes(G, node) {
        forall_out_edges(G, edge, node) {
            NodeID target = G.getEdgeTarget(edge);
            if (ind.solution[node] == 0 && ind.solution[target] == 0) {
                return false;
            }
        } endfor
    } endfor

    return true;
}

void population_mis::set_extract_fraction(double fraction) {
    remove_fraction = fraction;
}

