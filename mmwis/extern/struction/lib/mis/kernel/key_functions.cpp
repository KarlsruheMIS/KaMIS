//
// Created by alex on 31.01.20.
//
#include "key_functions.h"
#include "struction_branch_and_reduce_algorithm.h"

using namespace struction;

Gain ApproximateIncreaseKey::key(branch_and_reduce_algorithm* br_alg, NodeID n) {
    auto &status = br_alg->status;
    auto &neighbors = status.graph[n];
    auto &set = br_alg->set_1;

    size_t new_nodes_estimate = 0;
    for (NodeID u : status.graph[n]) {
        new_nodes_estimate += status.weights[n] < status.weights[u];
    }
    set.resize(status.graph.size());
    for (NodeID u : neighbors) {
        set.clear();
        for (NodeID v : status.graph[u]) {
            set.add(v);
        }
        for (NodeID v : neighbors) {
            new_nodes_estimate += u > v && !set.get(v) && status.weights[n] < status.weights[u] + status.weights[v];
        }
    }

    return key_by_set_estimate(br_alg, n, new_nodes_estimate);
}

Gain ApproximateIncreaseKey::key_by_set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, size_t set_estimate) {
    return (int) br_alg->deg(n) - (int) set_estimate;
}

size_t ApproximateIncreaseKey::set_limit(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return key_reinsert_factor * set_estimate(br_alg, n, key) + 1;
}

size_t ApproximateIncreaseKey::set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return br_alg->deg(n) - key;
}



Gain DegreeKey::key(branch_and_reduce_algorithm* br_alg, NodeID n) {
    return -br_alg->deg(n);
}

Gain DegreeKey::key_by_set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, size_t set_estimate) {
    return std::numeric_limits<Gain>::min();
}

size_t DegreeKey::set_limit(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return std::numeric_limits<size_t>::max() - 1;
}

size_t DegreeKey::set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return 0;
}


Gain IncreaseKey::key(branch_and_reduce_algorithm* br_alg, NodeID n) {
    auto &status = br_alg->status;
    auto &neighbors = status.graph[n];
    finder.findAllMWIS<false>(br_alg, neighbors, status.weights[n], 1000);
    size_t new_nodes_estimate = finder.get_sets().size();

    return key_by_set_estimate(br_alg, n, new_nodes_estimate);
}

Gain IncreaseKey::key_by_set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, size_t set_estimate) {
    return (int) br_alg->deg(n) - (int) set_estimate;
}

size_t IncreaseKey::set_limit(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return set_estimate(br_alg, n, key) + 1;
}

size_t IncreaseKey::set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return br_alg->deg(n) - key;
}


Gain RandomKey::key(branch_and_reduce_algorithm* br_alg, NodeID n) {
    return distribution(generator);
}

Gain RandomKey::key_by_set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, size_t set_estimate) {
    return std::numeric_limits<Gain>::min();
}

size_t RandomKey::set_limit(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return std::numeric_limits<size_t>::max() - 1;
}

size_t RandomKey::set_estimate(branch_and_reduce_algorithm* br_alg, NodeID n, Gain key) {
    return 0;
}
