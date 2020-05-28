/******************************************************************************
* branch_and_reduce_algorithm.h
*
*****************************************************************************/

#ifndef BRANCH_AND_REDUCE_SOLVER_H
#define BRANCH_AND_REDUCE_SOLVER_H

// local includes
#include "fast_set.h"
#include "mis_config.h"
#include "timer.h"
#include "definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/sized_vector.h"
#include "data_structure/dynamic_graph.h"
#include "reductions.h"

// system includes
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <iostream>
#include <sstream>

class cout_handler {
private:
	static std::streambuf* cout_rdbuf_backup;
	static std::stringstream buffered_output;
	static int disable_count;

	cout_handler() = delete;

public:
	static void disable_cout() {
		disable_count++;

		if (disable_count == 1) {
			buffered_output.str(std::string());
			buffered_output.clear();
			std::cout.rdbuf(buffered_output.rdbuf());
		}
	}

	static void enable_cout() {
		if (disable_count == 0)
			return;

		disable_count--;

		if (disable_count == 0)
			std::cout.rdbuf(cout_rdbuf_backup);
	}
};


class branch_and_reduce_algorithm {
public:
	enum IS_status { not_set, included, excluded, folded };

private:
	friend general_reduction;
	friend neighborhood_reduction;
	friend clique_neighborhood_reduction;
	friend clique_neighborhood_reduction_fast;
	friend critical_set_reduction;
	friend fold2_reduction;
	friend clique_reduction;
	friend twin_reduction;
	friend domination_reduction;
	friend generalized_neighborhood_reduction;
	friend generalized_fold_reduction;



	struct node_pos {
		NodeID node;
		size_t pos;

		node_pos(NodeID node = 0, size_t pos = 0) : node(node), pos(pos) {}
	};

	struct graph_status {
		size_t n = 0;
		size_t remaining_nodes = 0;
		NodeWeight is_weight = 0;
		NodeWeight reduction_offset = 0;
		dynamic_graph graph;
		std::vector<NodeWeight> weights;
		std::vector<IS_status> node_status;
		std::vector<reduction_ptr> reductions;
		sized_vector<reduction_type> folded_queue;
		sized_vector<node_pos> branching_queue;
		sized_vector<NodeID> modified_queue;

		graph_status() = default;

		graph_status(graph_access& G) :
			n(G.number_of_nodes()), remaining_nodes(n), graph(G), weights(n, 0), node_status(n, IS_status::not_set),
			folded_queue(n), branching_queue(n), modified_queue(n + 1) {

			forall_nodes(G, node) {
				weights[node] = G.getNodeWeight(node);
			} endfor
		}
	};

	static constexpr NodeID BRANCHING_TOKEN = std::numeric_limits<NodeID>::max();
	//static constexpr NodeID INITIAL_REDUCTION_TOKEN = BRANCHING_TOKEN - 1;

	static bool is_token(NodeID node) {
		//return node == BRANCHING_TOKEN || node == INITIAL_REDUCTION_TOKEN;
		return node == BRANCHING_TOKEN;
	}

	// lower graph size limit for when to use ils pruning
	static constexpr size_t ILS_SIZE_LIMIT = 50;

	// min number of remaining nodes to split up connected components
	static constexpr size_t SPLIT_CC_LIMIT = 100;

	// max number of neighbors to recurse on induced neighborhood graph in reductions
	static constexpr size_t REDU_RECURSION_LIMIT = 150;

	// max number of nodes to use all reductions during branch reduce
	static constexpr size_t FULL_REDUCTIONS_RECURSION_LIMIT = 50;

	MISConfig config;
	graph_status best_solution_status;
	NodeWeight best_weight = 0;
	timer t;
	bool is_ils_best_solution = false;
	size_t active_reduction_index;
	bool timeout = false;

	graph_status global_status;
	graph_access global_graph;
	std::vector<NodeID> global_mapping;
	std::vector<size_t> global_reduction_map;
	size_t total_ils_node_count;

	graph_status status;
	graph_access* local_graph;
	std::vector<NodeID> local_mapping;
	std::vector<size_t> local_reduction_map;
	std::vector<reduction_ptr> local_reductions;

	graph_access recursive_graph;
	std::vector<NodeID> recursive_mapping;
	std::vector<int> recursive_comp_map;
	std::vector<NodeID> recursive_local_mapping;

	std::function<void()> set_local_reductions;

	// presized onjects for temporary use
	fast_set set_1;
	fast_set set_2;
	fast_set double_set;
	sized_vector<sized_vector<NodeID>> buffers;


	size_t deg(NodeID node) const;
	void set(NodeID node, IS_status status, bool push_modified = true);
	void unset(NodeID node, bool restore = true);
	void flip_include_exclude(NodeID node);

	void fill_global_greedy();
	void compute_ils_pruning_bound();
	NodeWeight compute_cover_pruning_bound();

	void init_reduction_step();
	void add_next_level_node(NodeID node);
	void add_next_level_nodes(const std::vector<NodeID>& nodes);
	void add_next_level_nodes(const sized_vector<NodeID>& nodes);
	void add_next_level_neighborhood(NodeID node);
	void add_next_level_neighborhood(const std::vector<NodeID>& nodes);

	void reduce_graph_internal();
	bool branch_reduce_recursive();
	void branch_reduce_single_component();
	void initial_reduce();

	void update_best_solution();
	void reverse_branching();
	void restore_best_local_solution();
	void restore_best_global_solution();

	void build_global_graph_access();
	void build_induced_neighborhood_subgraph(graph_access& G, NodeID source_node);
	void build_induced_subgraph(graph_access& G, const sized_vector<NodeID>& nodes, const fast_set& nodes_set, sized_vector<NodeID>& reverse_mapping);

	void disable_cout();
	void enable_cout();

public:
	branch_and_reduce_algorithm(graph_access& G, const MISConfig& config, bool called_from_fold = false);

	void reduce_graph();
	bool run_branch_reduce();

	static size_t run_ils(const MISConfig& config, graph_access& G, sized_vector<NodeID>& tmp_buffer, size_t max_swaps);
	static void greedy_initial_is(graph_access& G, sized_vector<NodeID>& tmp_buffer);

	NodeWeight get_current_is_weight() const;
	void reverse_reduction(graph_access & G, graph_access & reduced_G, std::vector<NodeID> & reverse_mapping);
	void apply_branch_reduce_solution(graph_access & G);

	void build_graph_access(graph_access& G, std::vector<NodeID>& reverse_mapping) const;
};

#endif //BRANCH_AND_REDUCE_SOLVER_H
