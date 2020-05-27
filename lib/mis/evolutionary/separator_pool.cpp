/**
 * separator_pool.cpp
 * Purpose: Create and maintain a number of node separators and according partitions.
 *
 *****************************************************************************/

#include "separator_pool.h"

#include "mis_log.h"
#include "diversifier.h"

separator_pool::separator_pool() {
    
}

separator_pool::~separator_pool() {
    clear_partitions();
    clear_separators();
    clear_k_partitions();
    clear_k_separators();
    internal_separators.clear();
    internal_partitions.clear();
    internal_k_separators.clear();
    internal_k_partitions.clear();
    delete [] xadj;
    xadj = NULL;
    delete [] adjncy;
    adjncy = NULL;
}

void separator_pool::clear_partitions() {
    for (unsigned int i = 0; i < partitions_size; ++i) {
        delete [] internal_partitions[i].partition_map;
        internal_partitions[i].partition_map = NULL;
    }
}

void separator_pool::clear_graph() {
    delete [] xadj;
    xadj = NULL;
    delete [] adjncy;
    adjncy = NULL;
}

void separator_pool::clear_k_partitions() {
    for (unsigned int i = 0; i < k_partitions_size; ++i) {
        delete [] internal_k_partitions[i].partition_map;
        internal_k_partitions[i].partition_map = NULL;
    }
}

void separator_pool::clear_separators() {
    for (unsigned int i = 0; i < separators_size; ++i) {
        delete [] internal_separators[i].separator_map;
        internal_separators[i].separator_map = NULL;
    }
}

void separator_pool::clear_k_separators() {
    for (unsigned int i = 0; i < k_separators_size; ++i) {
        delete [] internal_k_separators[i].separator_map;
        internal_k_separators[i].separator_map = NULL;
    }
}

void separator_pool::calculate_scores(MISConfig & config, graph_access & G, population_mis & pop) {
    clear_scores(config);
    if (config.use_multiway_vc) calculate_partition_scores(config, G, pop);
    else calculate_separator_scores(config, G, pop);
}

void separator_pool::renew_pool(MISConfig & config, graph_access & G, bool init, bool scores, population_mis & pop) {
    mis_log::instance()->print_pool_title();
    mis_log::instance()->restart_building_pool_timer();
    if (!init) {
        clear_partitions();
        clear_separators();
        clear_k_partitions();
        clear_k_separators();
    }
    
    generate_partitions(config, G);
    generate_separators(config, G);
    if (config.use_multiway_vc) generate_k_partitions(config, G);
    else generate_k_separators(config, G);

    if (scores) calculate_scores(config, G, pop);
    mis_log::instance()->print_separator();
}

void separator_pool::init(MISConfig & config, graph_access & G) {
    xadj = G.UNSAFE_metis_style_xadj_array();
    adjncy = G.UNSAFE_metis_style_adjncy_array();
    separators_size = config.number_of_separators; 
    partitions_size = config.number_of_partitions; 
    k_separators_size = config.use_multiway_vc ? 0 : config.number_of_k_separators;
    k_partitions_size = config.use_multiway_vc ? config.number_of_k_partitions : 0;
    if (config.use_multiway_vc) scores.resize(config.multiway_blocks*config.population_size*config.number_of_k_partitions);
    else scores.resize(config.multiway_blocks*config.population_size*config.number_of_k_separators);
    clear_scores(config);
}

void separator_pool::generate_separators(MISConfig & config, graph_access & G) {
    internal_separators.clear();
    internal_separators.resize(config.number_of_separators);
    for (unsigned int i = 0; i < separators_size; ++i) {
        diversifier div;
        div.diversify(config);
        // Add randomization to the imbalance
        config.imbalance = random_functions::nextDouble(0.05, 0.75);
        separator sep; 
        forall_nodes(G, node) {
            G.setPartitionIndex(node, 0);
        } endfor
        create_separator(config, G, sep);
        sep.id = i;
        internal_separators[i] = sep;
    }
}

void separator_pool::generate_partitions(MISConfig & config, graph_access & G) {
    internal_partitions.clear();
    internal_partitions.resize(config.number_of_partitions);
    for (unsigned int i = 0; i < partitions_size; ++i) {
        diversifier div;
        div.diversify(config);
        // Add randomization to the imbalance
        config.imbalance = random_functions::nextDouble(0.05, 0.75);
        partition part; 
        forall_nodes(G, node) {
            G.setPartitionIndex(node, 0);
        } endfor
        create_partition(config, G, part);
        part.id = i;
        internal_partitions[i] = part;
    }
}

void separator_pool::generate_k_partitions(MISConfig & config, graph_access & G) {
    internal_k_partitions.clear();
    internal_k_partitions.resize(config.number_of_k_partitions);
    for (unsigned int i = 0; i < k_partitions_size; ++i) {
        diversifier div;
        div.diversify(config);
        // Add randomization to the imbalance
        config.imbalance = random_functions::nextDouble(0.05, 0.75);
        // Randomize the number of blocks
        int k = random_functions::nextInt(2, config.multiway_blocks);
        partition part; 
        forall_nodes(G, node) {
            G.setPartitionIndex(node, 0);
        } endfor
        create_partition(config, G, part, k);
        part.id = i;
        internal_k_partitions[i] = part;
    }
}

void separator_pool::generate_k_separators(MISConfig & config, graph_access & G) {
    internal_k_separators.clear();
    internal_k_separators.resize(config.number_of_k_separators);
    for (unsigned int i = 0; i < k_separators_size; ++i) {
        diversifier div;
        div.diversify(config);
        // Add randomization to the imbalance
        config.imbalance = random_functions::nextDouble(0.05, 0.75);
        // Randomize the number of blocks
        int k = random_functions::nextInt(2, config.multiway_blocks);
        separator sep; 
        forall_nodes(G, node) {
            G.setPartitionIndex(node, 0);
        } endfor
        create_separator(config, G, sep, k);
        sep.id = i;
        internal_k_separators[i] = sep;
    }
}

void separator_pool::create_separator(MISConfig & config, graph_access & G, separator & sep, int number_of_blocks) {
    int n = G.number_of_nodes();
    int k = number_of_blocks;
    int *separator = NULL;
    int separator_size;
    int* part = new int [G.number_of_nodes()];

    node_separator( &n, 
                    NULL, 
                    xadj,
                    NULL, 
                    adjncy,
                    &k, 
                    &config.imbalance, 
                    false, 
                    config.seed, 
                    config.kahip_mode, 
                    &separator_size, 
                    &separator,
                    part);


    std::cout << "Separator [Seed: " << config.seed << ", K: " << k << ", Size: " << separator_size << "]" << std::endl;

    forall_nodes(G, node) {
        G.setPartitionIndex(node, part[node]);
    } endfor

    for( int i = 0; i < separator_size; i++) {
        G.setPartitionIndex(separator[i], k);
    } 

    sep.k = k;
    sep.separator_size = separator_size;
    sep.separator_map = new int [G.number_of_nodes()];
    forall_nodes(G, node) {
        sep.separator_map[node] = G.getPartitionIndex(node);
    } endfor

    delete [] separator;
    delete [] part;
    part = NULL;
    separator = NULL;
}

void separator_pool::separator_bfs(MISConfig & config, graph_access & G, unsigned int k) {
    for (unsigned int i = 0; i < k - 1; ++i) {
        std::queue<NodeID> node_list;
        NodeID start = random_functions::nextInt(0, G.number_of_nodes() - 1);
        while (G.getPartitionIndex(start) != 0) start = (start + 1) % G.number_of_nodes();
        G.setPartitionIndex(start, i + 1);
        node_list.push(start);
        while (!node_list.empty()) {
            NodeID node = node_list.front();
            node_list.pop();
            forall_out_edges(G, edge, node) {
                NodeID target = G.getEdgeTarget(edge);
                if (G.getPartitionIndex(target) != 0) continue;
                else {
                    G.setPartitionIndex(target, i + 1);
                    node_list.push(target);
                }
            } endfor
        }
    }
}

void separator_pool::create_partition(MISConfig & config, graph_access & G, partition & part, int number_of_blocks) {
    int edgecut = 0;
    int n = G.number_of_nodes();
    int k = number_of_blocks;
    part.k = k;
    part.partition_map = new int [G.number_of_nodes()];

    kaffpa( &n, 
            NULL,
            xadj,
            NULL,
            adjncy,
            &k, 
            &config.imbalance,
            false,
            config.seed,
            config.kahip_mode,
            &edgecut,
            part.partition_map);

    std::cout << "Partition [Seed: " << config.seed << ", K: " << k << "]" << std::endl;
}

void separator_pool::calculate_partition_scores(MISConfig & config, graph_access & G, population_mis & pop) {
    for (unsigned int i = 0; i < config.population_size; ++i) {
        individuum_mis ind;
        pop.get_individuum(i, ind); 
        update_scores_for_individuum(config, G, ind);
    } 
    // print_scores(config);
}

void separator_pool::calculate_separator_scores(MISConfig & config, graph_access & G, population_mis & pop) {
    for (unsigned int i = 0; i < config.population_size; ++i) {
        individuum_mis ind;
        pop.get_individuum(i, ind); 
        update_scores_for_individuum(config, G, ind);
    } 
}

void separator_pool::update_scores_for_individuum(MISConfig & config, graph_access & G, individuum_mis & ind) {
    if (config.use_multiway_vc) {
        for (unsigned int j = 0; j < k_partitions_size; ++j) {
            for (unsigned k = 0; k < config.multiway_blocks; ++k) {
                unsigned int index = (config.population_size * config.number_of_k_partitions * k) 
                                        + config.number_of_k_partitions * ind.id
                                        + j;
                scores[index] = 0;
            }
        }
        for (unsigned int j = 0; j < k_partitions_size; ++j) {
            partition part = internal_k_partitions[j];
            apply_partition(config, G, part);
            calculate_score(config, G, ind, j);
        }
    } 
    else {
        for (unsigned int j = 0; j < k_separators_size; ++j) {
            for (unsigned k = 0; k < config.multiway_blocks; ++k) {
                unsigned int index = (config.population_size * config.number_of_k_separators * k) 
                                        + config.number_of_k_separators * ind.id
                                        + j;
                scores[index] = 0;
            }
        }
        for (unsigned int j = 0; j < k_separators_size; ++j) {
            separator sep = internal_k_separators[j];
            apply_separator(config, G, sep);
            calculate_score(config, G, ind, j);
        }
    }
}

void separator_pool::calculate_score(MISConfig & config, graph_access & G, individuum_mis & ind, unsigned int id) {
    if (config.use_multiway_vc) {
        forall_nodes(G, node) {
            unsigned int index = (config.population_size * config.number_of_k_partitions * G.getPartitionIndex(node)) 
                                    + config.number_of_k_partitions * ind.id 
                                    + id;
            if (ind.solution[node]) scores[index]++;
        } endfor
    }
    else { 
        forall_nodes(G, node) {
            if (G.getPartitionIndex(node) == G.get_partition_count()) continue;
            unsigned int index = (config.population_size * config.number_of_k_separators * G.getPartitionIndex(node)) 
                                    + config.number_of_k_separators * ind.id 
                                    + id;
            if (ind.solution[node]) scores[index]++;
        } endfor
    }
}

void separator_pool::clear_scores(MISConfig & config) {
    std::fill(scores.begin(), scores.end(), 0);
}

void separator_pool::apply_partition(MISConfig & config, graph_access & G, partition & part) {
    G.set_partition_count(part.k);
    forall_nodes(G, node) {
        G.setPartitionIndex(node, part.partition_map[node]);
    } endfor
}

void separator_pool::apply_separator(MISConfig & config, graph_access & G, separator & sep) {
    G.set_partition_count(sep.k);
    forall_nodes(G, node) {
        G.setPartitionIndex(node, sep.separator_map[node]);
    } endfor
}

void separator_pool::get_random_separator(separator & sep) {
    unsigned int rand = random_functions::nextInt(0, separators_size - 1); 
    sep = internal_separators[rand];
}

void separator_pool::get_random_k_separator(separator & sep) {
    unsigned int rand = random_functions::nextInt(0, k_separators_size - 1); 
    sep = internal_k_separators[rand];
}

void separator_pool::get_random_partition(partition & part) {
    unsigned int rand = random_functions::nextInt(0, partitions_size - 1);    
    part = internal_partitions[rand];
}

void separator_pool::get_random_k_partition(partition & part) {
    unsigned int rand = random_functions::nextInt(0, k_partitions_size - 1);    
    part = internal_k_partitions[rand];
}

void separator_pool::get_best_candidates(MISConfig & config, graph_access & G, population_mis & pop, std::vector<individuum_mis> & candidates, partition & part) {
    // print_scores(config);
    for (int i = 0; i < part.k; ++i) {
        unsigned int max_score = 0;
        for (unsigned int j = 0; j < config.population_size; ++j) {
            individuum_mis ind;
            pop.get_individuum(j, ind);
            unsigned int index = (config.population_size * config.number_of_k_partitions * i) 
                                    + config.number_of_k_partitions * j
                                    + part.id;
            unsigned int score = scores[index];
            if (score > max_score) {
                max_score = score;
                candidates[i] = ind;
            }
        }
        if (max_score == 0) {
            individuum_mis ind;
            pop.get_individuum(random_functions::nextInt(0, config.population_size - 1), ind);
            candidates[i] = ind;
        }
    }
}

void separator_pool::get_best_candidates(MISConfig & config, graph_access & G, population_mis & pop, std::vector<individuum_mis> & candidates, separator & sep) {
    // print_scores(config);
    for (int i = 0; i < sep.k; ++i) {
        unsigned int max_score = 0;
        for (unsigned int j = 0; j < config.population_size; ++j) {
            individuum_mis ind;
            pop.get_individuum(j, ind);
            unsigned int index = (config.population_size * config.number_of_k_separators * i) 
                                    + config.number_of_k_separators * j
                                    + sep.id;
            unsigned int score = scores[index];
            if (score > max_score) {
                max_score = score;
                candidates[i] = ind;
            }
        }
        if (max_score == 0) {
            individuum_mis ind;
            pop.get_individuum(random_functions::nextInt(0, config.population_size - 1), ind);
            candidates[i] = ind;
        }
    }
}

void separator_pool::print_separator_sizes() {
    for (unsigned int i = 0; i < separators_size; ++i) {
        printf("Separator size:\t%d\n", internal_separators[i].separator_size);
    }
}

void separator_pool::print_scores(MISConfig & config) {
    if (config.use_multiway_vc) {
        for (unsigned int i = 0; i < config.population_size; ++i) {
            printf("Ind: %d\n", i);
            for (unsigned int j = 0; j < k_partitions_size; ++j) {
                printf("|-Part (k=%d): %d\n", internal_k_partitions[j].k, j);
                for (unsigned int k = 0; k < config.multiway_blocks; ++k) {
                    unsigned int index = (config.population_size * config.number_of_k_partitions * k) 
                                            + config.number_of_k_partitions * i
                                            + j;
                    unsigned int score = scores[index];
                    printf("%d\t", score);
                }
                printf("\n");
            }
        }
    }
    else {
        for (unsigned int i = 0; i < config.population_size; ++i) {
            printf("Ind: %d\n", i);
            for (unsigned int j = 0; j < k_separators_size; ++j) {
                printf("|-Sep (k=%d): %d\n", internal_k_separators[j].k, j);
                for (unsigned int k = 0; k < config.multiway_blocks; ++k) {
                    unsigned int index = (config.population_size * config.number_of_k_separators * k) 
                                            + config.number_of_k_separators * i
                                            + j;
                    unsigned int score = scores[index];
                    printf("%d\t", score);
                }
                printf("\n");
            }
        }
    }
}

