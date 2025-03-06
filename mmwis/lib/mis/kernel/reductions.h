/******************************************************************************
* reductions.h
*
*****************************************************************************/

#pragma once

// local includes
#include "definitions.h"
#include "vertex_marker.h"
#include "fast_set.h"
#include "sized_vector.h"
#include "dynamic_graph.h"

// system includes
#include <cstddef>
#include <vector>
#include <memory>
#include <array>

namespace mmwis {

class branch_and_reduce_algorithm;

enum reduction_type {fold1, v_shape, v_shape_min, triangle, single_edge, extended_single_edge, neighborhood, clique, twin, heavy_set, clique_neighborhood_fast, clique_neighborhood, critical_set, generalized_neighborhood, generalized_fold};

constexpr size_t REDUCTION_NUM = 17;

struct general_reduction {
	general_reduction(size_t n) : marker(n) {}
	virtual ~general_reduction() {}
	virtual general_reduction* clone() const = 0;

	virtual reduction_type get_reduction_type() const = 0;
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) = 0;
	virtual void restore(branch_and_reduce_algorithm* br_alg) {}
	virtual void apply(branch_and_reduce_algorithm* br_alg) {}

	virtual void print_reduction_type() {};
	bool has_run = false;
	vertex_marker marker;
};

struct neighborhood_reduction : public general_reduction {
	neighborhood_reduction(size_t n) : general_reduction(n) {}
	~neighborhood_reduction() {}
	virtual neighborhood_reduction* clone() const final { return new neighborhood_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::neighborhood; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "neighborhood_reduction" << std::endl; }
};

struct clique_neighborhood_reduction_fast : public general_reduction {
	clique_neighborhood_reduction_fast(size_t n) : general_reduction(n) {}
	~clique_neighborhood_reduction_fast() {}
	virtual clique_neighborhood_reduction_fast* clone() const final { return new clique_neighborhood_reduction_fast(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood_fast; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "clique_neighborhood_reduction_fast" << std::endl; }
};

struct clique_neighborhood_reduction : public general_reduction {
	clique_neighborhood_reduction(size_t n) : general_reduction(n) {}
	~clique_neighborhood_reduction() {}
	virtual clique_neighborhood_reduction* clone() const final { return new clique_neighborhood_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "clique_neighborhood_reduction" << std::endl; }

	bool partition_into_cliques(NodeID v);
	bool expand_clique(NodeID max_neighbor, sized_vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set);

	branch_and_reduce_algorithm* br_alg;
	NodeWeight target_weight;
	NodeWeight neighbor_weights;
};

struct critical_set_reduction : public general_reduction {
	critical_set_reduction(size_t n) : general_reduction(n) {}
	~critical_set_reduction() {}
	virtual critical_set_reduction* clone() const final { return new critical_set_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::critical_set; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "critical_set_reduction" << std::endl; }
};


struct fold1_reduction : public general_reduction {
	fold1_reduction(size_t n) : general_reduction(n) {}
	~fold1_reduction() {}
	virtual fold1_reduction* clone() const final { return new fold1_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::fold1; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "fold1_reduction" << std::endl; }
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
	struct fold_nodes {
		NodeID deg1_node;
		NodeID fold_node;
	};

	struct restore_data {
		fold_nodes nodes;
		NodeWeight deg1_weight;
		dynamic_graph::neighbor_list old_neighbors;
	};

	void fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);

	std::vector<restore_data> restore_vec;
};


struct clique_reduction : public general_reduction {
	clique_reduction(size_t n) : general_reduction(n) {}
	~clique_reduction() {}
	virtual clique_reduction* clone() const final { return new clique_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::clique; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "clique_reduction" << std::endl; }

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

struct triangle_reduction : public general_reduction {
	triangle_reduction(size_t n) : general_reduction(n) {}
	~triangle_reduction() {}
	virtual triangle_reduction* clone() const final { return new triangle_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::triangle; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
	struct fold_nodes {
        NodeID deg2_node;
        NodeID bigger;
        NodeID smaller;
	};

	struct restore_data {
		fold_nodes nodes;
		NodeWeight deg2_weight;
        int fold_case; // 0->fold_min, 1->fold_mid
	};

	void fold_mid_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);
	void fold_min_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);

	std::vector<restore_data> restore_vec;
};

struct v_shape_reduction : public general_reduction {
	v_shape_reduction(size_t n) : general_reduction(n) {}
	~v_shape_reduction() {}
	virtual v_shape_reduction* clone() const final { return new v_shape_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::v_shape; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "v_shape_reduction" << std::endl; }

private:
	struct fold_nodes {
        NodeID deg2_node;
        std::vector<NodeID> neighbors;
	};

	struct restore_data {
		fold_nodes nodes;
		NodeWeight deg2_weight;
        int fold_case;
        dynamic_graph::neighbor_list neighbor_list;
        std::array<std::vector<NodeID>, 2> node_vecs;
	};

	void fold_mid_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);
	void fold_max_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);

	std::vector<restore_data> restore_vec;
};


struct v_shape_min_reduction : public general_reduction {
	v_shape_min_reduction(size_t n) : general_reduction(n) {}
	~v_shape_min_reduction() {}
	virtual v_shape_min_reduction* clone() const final { return new v_shape_min_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::v_shape_min; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "v_shape_min_reduction" << std::endl; }

private:
	struct fold_nodes {
        NodeID deg2_node;
        std::vector<NodeID> neighbors;
	};

	struct restore_data {
		fold_nodes nodes;
		NodeWeight deg2_weight;
        std::array<std::vector<NodeID>, 2> node_vecs;
	};

	void fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes);
	std::vector<restore_data> restore_vec;
};

struct single_edge_reduction: public general_reduction {
	single_edge_reduction(size_t n) : general_reduction(n) {}
	~single_edge_reduction() {}
	virtual single_edge_reduction* clone() const final { return new single_edge_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::single_edge; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "single_edge_reduction" << std::endl; }
};

struct extended_single_edge_reduction: public general_reduction {
	extended_single_edge_reduction(size_t n) : general_reduction(n) {}
	~extended_single_edge_reduction() {}
	virtual extended_single_edge_reduction* clone() const final { return new extended_single_edge_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::extended_single_edge; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "extended_single_edge_reduction" << std::endl; }
};


struct twin_reduction : public general_reduction {
	twin_reduction(size_t n) : general_reduction(n) {}
	~twin_reduction() {}
	virtual twin_reduction* clone() const final { return new twin_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::twin; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "twin_reduction" << std::endl; }

private:
	struct restore_data {
		NodeID main;
		NodeID twin;
	};

	void fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin);

	std::vector<restore_data> restore_vec;
};


struct heavy_set_reduction : public general_reduction {
	heavy_set_reduction(size_t n) : general_reduction(n) {}
	~heavy_set_reduction() {}
	virtual heavy_set_reduction* clone() const final { return new heavy_set_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_set; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "heavy_set_reduction" << std::endl; }
};

struct generalized_neighborhood_reduction : public general_reduction {
	generalized_neighborhood_reduction(size_t n) : general_reduction(n) {}
	~generalized_neighborhood_reduction() {}
	virtual generalized_neighborhood_reduction* clone() const final { return new generalized_neighborhood_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::generalized_neighborhood; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "generalized_neighborhood_reduction" << std::endl; }
};

struct generalized_fold_reduction : public general_reduction {
	generalized_fold_reduction(size_t n) : general_reduction(n) {}
	~generalized_fold_reduction() {}
	virtual generalized_fold_reduction* clone() const final { return new generalized_fold_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::generalized_fold; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;
	virtual void print_reduction_type() final {std::cout << "generalized_fold_reduction" << std::endl; }

private:
	struct fold_nodes {
		NodeID main;
		std::vector<NodeID> MWIS;
	};

	struct restore_data {
		fold_nodes nodes;
		NodeWeight main_weight;
		NodeWeight MWIS_weight;
		dynamic_graph::neighbor_list main_neighbor_list;
		std::vector<std::vector<NodeID>> MWIS_node_vecs;
	};

	void fold(branch_and_reduce_algorithm* br_alg, NodeID main_node, fast_set& MWIS_set, NodeWeight MWIS_weight);

	std::vector<restore_data> restore_vec;
};


struct reduction_ptr {
	general_reduction* reduction = nullptr;

	reduction_ptr() = default;

	~reduction_ptr() {
		release();
	}

	reduction_ptr(general_reduction* reduction) : reduction(reduction) {};

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

	general_reduction* operator->() const {
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


/*template<class ...Redus>
std::vector<std::unique_ptr<general_reduction>> make_reduction_vector(size_t n) {
std::vector<std::unique_ptr<general_reduction>> vec;
(vec.push_back(std::make_unique<Redus>(n)), ...);
return vec;
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class Last>
void make_reduction_vector_helper(std::vector<std::unique_ptr<general_reduction>>& vec, size_t n) {
vec.push_back(make_unique<Last>(n));
};

template<class First, class Second, class ...Redus>
void make_reduction_vector_helper(std::vector<std::unique_ptr<general_reduction>>& vec, size_t n) {
vec.push_back(make_unique<First>(n));
make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template<class ...Redus>
std::vector<std::unique_ptr<general_reduction>> make_reduction_vector(size_t n) {
std::vector<std::unique_ptr<general_reduction>> vec;
make_reduction_vector_helper<Redus...>(vec, n);
return vec;
};*/

}

