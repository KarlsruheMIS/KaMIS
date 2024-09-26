/******************************************************************************
* struction_reductions.h
*
* Copyright (C) 2015-2018 Robert Williger
*
* This program is free software: you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the Free
* Software Foundation, either version 2 of the License, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
* more details.
*
* You should have received a copy of the GNU General Public License along with
* this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _struction_REDUCTIONS_H
#define _struction_REDUCTIONS_H

// local includes
#include "definitions.h"
#include "fast_set.h"
#include "sized_vector.h"
#include "struction_dynamic_graph.h"
#include "extended_struction.h"
#include "original_struction.h"
#include "struction_vertex_marker.h"
#include "key_functions.h"
#include "maxNodeHeap.h"
#include "mmwis_config.h"

// system includes
#include <vector>
#include <memory>
#include <array>
#include <data_structure/priority_queues/MaxHeap.h>

namespace struction {

class branch_and_reduce_algorithm;
using Struction_Type = ::mmwis::MISConfig::Struction_Type;
using Key_Type = ::mmwis::MISConfig::Key_Type;

enum reduction_type { struction_neighborhood, struction_fold2, struction_clique, struction_clique_neighborhood_fast, struction_clique_neighborhood, struction_twin,
        struction_domination, struction_generalized_neighborhood, path, struction_decrease, struction_plateau, struction_blow, struction_critical_set, struction_generalized_fold};
constexpr size_t STRUCTION_REDUCTION_NUM = 14;

struct struction_general_reduction {
	struction_general_reduction(size_t n) : marker(n) {}
	virtual ~struction_general_reduction() {}
	virtual struction_general_reduction* clone() const = 0;

	virtual reduction_type get_reduction_type() const = 0;
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) = 0;
	virtual void restore(branch_and_reduce_algorithm* br_alg) {}
    virtual void apply(branch_and_reduce_algorithm* br_alg) {}
    virtual void reset(branch_and_reduce_algorithm* br_alg, size_t comp_size) {}

	bool has_run = false;
	struction_vertex_marker marker;

	template<typename F>
	void for_each_changed_vertex(branch_and_reduce_algorithm* br_alg, F f);
};

struct struction_neighborhood_reduction : public struction_general_reduction {
    struction_neighborhood_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_neighborhood_reduction() {}
    virtual struction_neighborhood_reduction* clone() const final { return new struction_neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_neighborhood; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct struction_clique_neighborhood_reduction_fast : public struction_general_reduction {
    struction_clique_neighborhood_reduction_fast(size_t n) : struction_general_reduction(n) {}
    ~struction_clique_neighborhood_reduction_fast() {}
    virtual struction_clique_neighborhood_reduction_fast* clone() const final { return new struction_clique_neighborhood_reduction_fast(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_clique_neighborhood_fast; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct struction_clique_neighborhood_reduction : public struction_general_reduction {
    struction_clique_neighborhood_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_clique_neighborhood_reduction() {}
    virtual struction_clique_neighborhood_reduction* clone() const final { return new struction_clique_neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_clique_neighborhood; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;

    bool partition_into_cliques(NodeID v);
    bool expand_clique(NodeID max_neighbor, sized_vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set);

    branch_and_reduce_algorithm* br_alg;
    NodeWeight target_weight;
    NodeWeight neighbor_weights;
};

struct struction_critical_set_reduction : public struction_general_reduction {
    struction_critical_set_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_critical_set_reduction() {}
    virtual struction_critical_set_reduction* clone() const final { return new struction_critical_set_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_critical_set; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct struction_fold2_reduction : public struction_general_reduction {
    struction_fold2_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_fold2_reduction() {}
    virtual struction_fold2_reduction* clone() const final { return new struction_fold2_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_fold2; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final;
    virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
    struct fold_nodes {
        NodeID main;
        std::array<NodeID, 2> rest;
    };

    struct restore_data {
        fold_nodes nodes;
        NodeWeight main_weight;
        struction_dynamic_graph::neighbor_list main_neighbor_list;
        std::array<std::vector<NodeID>, 2> node_vecs;
    };

    void fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);

    std::vector<restore_data> restore_vec;
};

struct struction_clique_reduction : public struction_general_reduction {
    struction_clique_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_clique_reduction() {}
    virtual struction_clique_reduction* clone() const final { return new struction_clique_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_clique; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final;
    virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
    struct weighted_node {
        NodeID node;
        NodeWeight weight;
    };

    struct restore_data {
        weighted_node isolated;
        std::vector<NodeID> non_isolated;

        restore_data() = default;
        restore_data(const weighted_node& isolated, std::vector<NodeID>&& non_isolated) : isolated(isolated), non_isolated(std::move(non_isolated)) {};
    };

    void fold(branch_and_reduce_algorithm* br_alg, const weighted_node& isolated, std::vector<NodeID>&& non_isolated);

    std::vector<restore_data> restore_vec;
};

struct struction_twin_reduction : public struction_general_reduction {
    struction_twin_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_twin_reduction() {}
    virtual struction_twin_reduction* clone() const final { return new struction_twin_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_twin; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final;
    virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
    struct restore_data {
        NodeID main;
        NodeID twin;
    };

    void fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin);

    std::vector<restore_data> restore_vec;
};

struct struction_domination_reduction : public struction_general_reduction {
    struction_domination_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_domination_reduction() {}
    virtual struction_domination_reduction* clone() const final { return new struction_domination_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_domination; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct struction_generalized_neighborhood_reduction : public struction_general_reduction {
    struction_generalized_neighborhood_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_generalized_neighborhood_reduction() {}
    virtual struction_generalized_neighborhood_reduction* clone() const final { return new struction_generalized_neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_generalized_neighborhood; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};


struct struction_generalized_fold_reduction : public struction_general_reduction {
    struction_generalized_fold_reduction(size_t n) : struction_general_reduction(n) {}
    ~struction_generalized_fold_reduction() {}
    virtual struction_generalized_fold_reduction* clone() const final { return new struction_generalized_fold_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::struction_generalized_fold; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final;
    virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
    struct fold_nodes {
        NodeID main;
        std::vector<NodeID> MWIS;
    };

    struct restore_data {
        fold_nodes nodes;
        NodeWeight main_weight;
        NodeWeight MWIS_weight;
        struction_dynamic_graph::neighbor_list main_neighbor_list;
        std::vector<std::vector<NodeID>> MWIS_node_vecs;
    };

    void fold(branch_and_reduce_algorithm* br_alg, NodeID main_node, fast_set& MWIS_set, NodeWeight MWIS_weight);

    std::vector<restore_data> restore_vec;
};

struct path_reduction : public struction_general_reduction {
    path_reduction(size_t n) : struction_general_reduction(n) {}
    ~path_reduction() {}
    virtual path_reduction* clone() const final { return new path_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::path; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final;
    virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
    branch_and_reduce_algorithm* br_alg;
    NodeWeight w_i_i, w_i_e, w_e_i, w_e_e; // c_{f}_{l} ^= weight of maximum weight IS on path (in/ex)cluding first node (f) and (in/ex)cluding last node (l)
    std::vector<std::pair<NodeID,NodeWeight>> node_weights; //tmp vector to store changed node weights
    struct restore_data {
        std::vector<std::pair<NodeID,NodeWeight>> node_weights;
        sized_vector<NodeID> path;
        size_t start;
        size_t end;
        size_t offset;
        bool relink = false;
    };

    std::vector<restore_data> restore_vec;

    template<bool add_global = true>
    void enqueue_node(NodeID n);
    bool dequeue_node(sized_vector<NodeID> &queue, NodeID &n, size_t degree);
    bool reduce_degree_one_node();
    bool reduce_path();

    void fold_path(size_t start, size_t end, restore_data &data);
    void reassign_weight(NodeID n, NodeWeight w, restore_data &data);
    void reconnect(NodeID v, NodeID w, restore_data &data);
    void add_reduction_offset(size_t reduction_offset, restore_data &data);
    bool are_connected(NodeID v, NodeID w) const;

    void find_max_deg_1_path(NodeID n, sized_vector<NodeID> &path);
    void find_max_path(NodeID n, sized_vector<NodeID> &path);
    bool next_node_on_path(NodeID current, NodeID last, NodeID first, NodeID &next);

    void find_MIS_on_path(sized_vector<NodeID> &path);

    template<bool track_choices = false>
    void find_MIS_on_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path);
    template<bool track_choices = false>
    void find_MIS_on_deg_1_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path);

    void apply(sized_vector<NodeID> &path);
};

template<typename struction_type, reduction_type type, int vertex_increase>
struct iterative_struction : public struction_general_reduction {
    iterative_struction(size_t n) : struction_general_reduction(n) {}
    ~iterative_struction() {}
    virtual iterative_struction* clone() const final { return new iterative_struction(*this); }

    virtual reduction_type get_reduction_type() const final { return type; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final {s.restore(br_alg);};
    virtual void apply(branch_and_reduce_algorithm* br_alg) final {s.apply(br_alg);};

private:
    struction_type s;
};

using reduction_struction = iterative_struction<extended_struction<false>, reduction_type ::struction_decrease, -1>;
using plateu_struction = iterative_struction<extended_struction<false>, reduction_type ::struction_plateau, 0>;

template<typename key_function = ApproximateIncreaseKey, typename struction_type = extended_struction<false>, reduction_type type = reduction_type ::struction_blow/**/>
struct blow_up_struction : public struction_general_reduction {
    blow_up_struction(const mmwis::MISConfig &config, size_t n) : struction_general_reduction(n), update_set(n), distribution(0, 1),
                                                     generator(config.seed), f(config) {}
    ~blow_up_struction() {}
    virtual blow_up_struction* clone() const final { return new blow_up_struction(*this); }

    virtual reduction_type get_reduction_type() const final { return type; }
    virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
    virtual void restore(branch_and_reduce_algorithm* br_alg) final;
    virtual void apply(branch_and_reduce_algorithm* br_alg) final {s.apply(br_alg);};
    virtual void reset(branch_and_reduce_algorithm* br_alg, size_t comp_size) final;

private:
    size_t phase_start_kernel;

    fast_set update_set;
    std::vector<std::pair<NodeID, Gain>> update_list;
    std::vector<NodeID> added_list;

    key_function f;
    struction_type s;
    bool restored = false;
    size_t blow_ups;
    branch_and_reduce_algorithm *br_alg;

    MaxHeap<float> queue;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution;

    void init_blow_up_phase();

    bool clean_up_queue();
    bool is_done();

    Gain denoise(float key) {
        return static_cast<Gain>(key);
    }

    float apply_noise(Gain key) {
        return key + distribution(generator);
    }

    void update_queue_by_key(NodeID n, Gain key) {
        if (queue.contains(n)) {
            if (blow_ups != 0 && update_set.add(n))
                update_list.emplace_back(n, denoise(queue.getKey(n)));
            queue.changeKey(n, apply_noise(key));
        } else {
            if (blow_ups != 0 && update_set.add(n))
                added_list.push_back(n);
            queue.insert(n, apply_noise(key));
        }
    }

    void update_queue(NodeID n);

};

struct reduction_ptr {
    struction_general_reduction* reduction = nullptr;

    reduction_ptr() = default;

    ~reduction_ptr() {
        release();
    }

    reduction_ptr(struction_general_reduction* reduction) : reduction(reduction) {};

    reduction_ptr(const reduction_ptr& other) : reduction(other.reduction->clone()) {};

    reduction_ptr& operator=(const reduction_ptr& other) {
        release();
        reduction = other.reduction->clone();
        return *this;
    };

    reduction_ptr(reduction_ptr&& other) : reduction(std::move(other.reduction)) {
        other.reduction = nullptr;
    };

    reduction_ptr& operator=(reduction_ptr&& other) {
        reduction = std::move(other.reduction);
        other.reduction = nullptr;
        return *this;
    };

    struction_general_reduction* operator->() const {
        return reduction;
    }

    void release() {
        if (reduction) {
            delete reduction;
            reduction = nullptr;
        }
    };
};

template<class Last>
void make_reduction_vector_helper(std::vector<reduction_ptr>& vec, size_t n) {
    vec.emplace_back(new Last(n));
};

template<class First, class Second, class ...Redus>
void make_reduction_vector_helper(std::vector<reduction_ptr>& vec, size_t n) {
    vec.emplace_back(new First(n));
    make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template<class ...Redus>
std::vector<reduction_ptr> make_reduction_vector(size_t n) {
    std::vector<reduction_ptr> vec;
    make_reduction_vector_helper<Redus...>(vec, n);
    return vec;
};
 
template<reduction_type type, int new_nodes>
reduction_ptr make_iterative_struction(const ::mmwis::MISConfig &config, size_t n);

reduction_ptr make_decreasing_struction(const ::mmwis::MISConfig &config, size_t n);

reduction_ptr make_plateau_struction(const ::mmwis::MISConfig &config, size_t n);

reduction_ptr make_increasing_struction(const ::mmwis::MISConfig &config, size_t n);

/*template<class ...Redus>
std::vector<std::unique_ptr<struction_general_reduction>> make_reduction_vector(size_t n) {
std::vector<std::unique_ptr<struction_general_reduction>> vec;
(vec.push_back(std::make_unique<Redus>(n)), ...);
return vec;
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class Last>
void make_reduction_vector_helper(std::vector<std::unique_ptr<struction_general_reduction>>& vec, size_t n) {
vec.push_back(make_unique<Last>(n));
};

template<class First, class Second, class ...Redus>
void make_reduction_vector_helper(std::vector<std::unique_ptr<struction_general_reduction>>& vec, size_t n) {
vec.push_back(make_unique<First>(n));
make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template<class ...Redus>
std::vector<std::unique_ptr<struction_general_reduction>> make_reduction_vector(size_t n) {
std::vector<std::unique_ptr<struction_general_reduction>> vec;
make_reduction_vector_helper<Redus...>(vec, n);
return vec;
};*/

} // end namespace struction
#endif //struction_REDUCTIONS_H
