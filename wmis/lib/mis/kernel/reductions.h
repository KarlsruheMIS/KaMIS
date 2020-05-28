/******************************************************************************
* reductions.h
*
*****************************************************************************/

#ifndef REDUCTIONS_H
#define REDUCTIONS_H

// local includes
#include "definitions.h"
#include "fast_set.h"
#include "data_structure/sized_vector.h"
#include "data_structure/dynamic_graph.h"

// system includes
#include <vector>
#include <memory>
#include <array>

class branch_and_reduce_algorithm;

enum reduction_type { neighborhood, fold2, clique, critical_set, clique_neighborhood_fast, clique_neighborhood, twin, domination, generalized_neighborhood, generalized_fold };
constexpr size_t REDUCTION_NUM = 10;

class vertex_marker {
private:
	sized_vector<NodeID> current;
	sized_vector<NodeID> next;
	fast_set added_vertices;

public:
	vertex_marker(size_t size) : current(size), next(size), added_vertices(size) {};

	void add(NodeID vertex) {
		if (!added_vertices.get(vertex)) {
			next.push_back(vertex);
			added_vertices.add(vertex);
		}
	}

	void get_next() {
		if (next.size() != 0) {
			current.swap(next);
			clear_next();
		}
	}

	void clear_next() {
		next.clear();
		added_vertices.clear();
	}

	void fill_current_ascending(size_t n) {
		current.clear();
		for (size_t i = 0; i < n; i++) {
			current.push_back(i);
		}
	}

	NodeID current_vertex(size_t index) {
		return current[index];
	}

	size_t current_size() {
		return current.size();
	}
};

struct general_reduction {
	general_reduction(size_t n) : marker(n) {}
	virtual ~general_reduction() {}
	virtual general_reduction* clone() const = 0;

	virtual reduction_type get_reduction_type() const = 0;
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) = 0;
	virtual void restore(branch_and_reduce_algorithm* br_alg) {}
	virtual void apply(branch_and_reduce_algorithm* br_alg) {}

	bool has_run = false;
	vertex_marker marker;
};

struct neighborhood_reduction : public general_reduction {
	neighborhood_reduction(size_t n) : general_reduction(n) {}
	~neighborhood_reduction() {}
	virtual neighborhood_reduction* clone() const final { return new neighborhood_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::neighborhood; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct clique_neighborhood_reduction_fast : public general_reduction {
	clique_neighborhood_reduction_fast(size_t n) : general_reduction(n) {}
	~clique_neighborhood_reduction_fast() {}
	virtual clique_neighborhood_reduction_fast* clone() const final { return new clique_neighborhood_reduction_fast(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood_fast; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct clique_neighborhood_reduction : public general_reduction {
	clique_neighborhood_reduction(size_t n) : general_reduction(n) {}
	~clique_neighborhood_reduction() {}
	virtual clique_neighborhood_reduction* clone() const final { return new clique_neighborhood_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;

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
};

struct fold2_reduction : public general_reduction {
	fold2_reduction(size_t n) : general_reduction(n) {}
	~fold2_reduction() {}
	virtual fold2_reduction* clone() const final { return new fold2_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::fold2; }
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
		dynamic_graph::neighbor_list main_neighbor_list;
		std::array<std::vector<NodeID>, 2> node_vecs;
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

struct twin_reduction : public general_reduction {
	twin_reduction(size_t n) : general_reduction(n) {}
	~twin_reduction() {}
	virtual twin_reduction* clone() const final { return new twin_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::twin; }
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

struct domination_reduction : public general_reduction {
	domination_reduction(size_t n) : general_reduction(n) {}
	~domination_reduction() {}
	virtual domination_reduction* clone() const final { return new domination_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::domination; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct generalized_neighborhood_reduction : public general_reduction {
	generalized_neighborhood_reduction(size_t n) : general_reduction(n) {}
	~generalized_neighborhood_reduction() {}
	virtual generalized_neighborhood_reduction* clone() const final { return new generalized_neighborhood_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::generalized_neighborhood; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
};

struct generalized_fold_reduction : public general_reduction {
	generalized_fold_reduction(size_t n) : general_reduction(n) {}
	~generalized_fold_reduction() {}
	virtual generalized_fold_reduction* clone() const final { return new generalized_fold_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::generalized_fold; }
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


#endif //REDUCTIONS_H
