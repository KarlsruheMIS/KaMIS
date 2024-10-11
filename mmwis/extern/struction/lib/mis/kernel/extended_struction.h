//
// Created by alex on 14.12.19.
//

#ifndef COMPONENTS_EXTENDED_STRUCTION_H
#define COMPONENTS_EXTENDED_STRUCTION_H

#include "mwis_finder.h"

namespace struction {
class branch_and_reduce_algorithm;

template<bool reduced = false>
class extended_struction {
public:
    bool reduce(branch_and_reduce_algorithm* br_alg, NodeID n, size_t max_nodes = std::numeric_limits<size_t>::max() - 1);
    void apply(branch_and_reduce_algorithm* br_alg);
    void restore(branch_and_reduce_algorithm* br_alg);
    size_t removed_vertices(branch_and_reduce_algorithm* br_alg, NodeID n);

private:
    branch_and_reduce_algorithm* br_alg;
    mwis_finder finder;
    struct restore_data {
        explicit restore_data(NodeID main, size_t set_nodes) : main(main), set_nodes(set_nodes) {}
        NodeID main;
        size_t set_nodes;
    };
    std::vector<restore_data> restore_vec;
    neighbor_list includable_neighbors;

    bool findAdditionalNodes(NodeID n, const std::vector<mwis> &independent_sets,
                             sized_vector<NodeID> &neighbor_ids, sized_vector<NodeID> &layer_prefix_count,
                             sized_vector<NodeID> &additional_nodes, size_t max_nodes);
};
}

#endif //COMPONENTS_EXTENDED_STRUCTION_H
