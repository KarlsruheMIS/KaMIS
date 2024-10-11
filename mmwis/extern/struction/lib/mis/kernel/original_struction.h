//
// Created by alex on 16.01.20.
//

#ifndef COMPONENTS_ORIGINAL_STRUCTION_H
#define COMPONENTS_ORIGINAL_STRUCTION_H

namespace struction {
class branch_and_reduce_algorithm;
template<bool modified = false>
class original_struction {

public:
    bool reduce(branch_and_reduce_algorithm* br_alg, NodeID n, size_t max_nodes);
    void apply(branch_and_reduce_algorithm* br_alg) {}
    void restore(branch_and_reduce_algorithm* br_alg) {}
    size_t removed_vertices(branch_and_reduce_algorithm* br_alg, NodeID n) {return 1;}

private:
    branch_and_reduce_algorithm* br_alg;
    struct restore_data {
        explicit restore_data(NodeID main, size_t set_nodes) : main(main), set_nodes(set_nodes) {}
        NodeID main;
        size_t set_nodes;
    };
    std::vector<restore_data> restore_vec;

    bool find_additional_nodes(sized_vector<NodeID> &additional_nodes, neighbor_list &neighbors, sized_vector<NodeID> &layer_prefix_count, size_t max_nodes) const;
};

}

#endif //COMPONENTS_ORIGINAL_STRUCTION_H
