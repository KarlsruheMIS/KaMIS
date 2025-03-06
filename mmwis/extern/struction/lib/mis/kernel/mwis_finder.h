//
// Created by alex on 14.12.19.
//

#ifndef COMPONENTS_MWIS_FINDER_H
#define COMPONENTS_MWIS_FINDER_H


#include <vector>
#include <data_structure/sized_vector.h>
#include "definitions.h"
#include "struction_dynamic_graph.h"


namespace struction {

using neighbor_list =struction_dynamic_graph::neighbor_list;
class branch_and_reduce_algorithm;

struct mwis {
    std::vector<NodeID> nodes;
    NodeWeight weight = 0;

    void push_back(NodeID n, NodeWeight w) {
        weight += w;
        nodes.push_back(n);
    }

    void pop_back(NodeWeight w) {
        weight -= w;
        nodes.pop_back();
    }
};

class mwis_finder {

public:

    std::vector<mwis> &get_sets();
    template<bool just_greater=false>
    bool findAllMWIS(branch_and_reduce_algorithm* br_alg, neighbor_list &nodes, NodeWeight min_weight, size_t max_sets = std::numeric_limits<size_t>::max() - 1);

private:
    static constexpr size_t NOT_INCLUDABLE = 0;
    static constexpr size_t INCLUDABLE = 1;

    std::vector<mwis> independent_sets;
    mwis tmp_mwis;

    branch_and_reduce_algorithm *br_alg;
    NodeWeight min_weight;
    size_t max_sets;


    template<bool just_greater=false>
    void findAllMWISRec(neighbor_list &nodes, size_t cur_node, NodeWeight remaining_neighbor_weight);

};
}
#endif //COMPONENTS_MWIS_FINDER_H
