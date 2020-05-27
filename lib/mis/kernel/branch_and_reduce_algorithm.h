 /******************************************************************************
 * branch_and_reduce_algorithm.h
 *
  *****************************************************************************/

#ifndef BRANCH_AND_REDUCE_SOLVER_H
#define BRANCH_AND_REDUCE_SOLVER_H

// local includes
#include "fast_set.h"
#include "modified.h"
#include "mis_config.h"

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "timer.h"

// system includes
#include <vector>
#include <string>
#include <memory>

class branch_and_reduce_algorithm
{
friend class modified;
friend class fold;
friend class alternative;

public:
	static int REDUCTION;
	static int LOWER_BOUND;
	static int BRANCHING;
	static bool outputLP;
	
	std::vector<std::vector<int>> adj;
	static long nBranchings;
	static int debug;
	double SHRINK;
	int    depth;
    int    maxDepth;
    int    rootDepth;
	int    n;
    int    N;

	/**
	 * current best solution
	 */
	int opt;
    std::vector<int> y;

	/**
	 * current solution (-1: not determined, 0: not in the vc, 1: in the vc, 2: removed by foldings)
	 */
	int crt;
    std::vector<int> x;

	/**
	 * #remaining vertices
	 */
	int rn;
	
	/**
	 * max flow
	 */
	std::vector<int> in;
    std::vector<int> out;
	
	/**
	 * lower bound
	 */
	int lb;

	std::vector<int> que, level, iter;

	std::vector<int> modTmp;

	std::vector<std::shared_ptr<modified>> modifieds;
	int modifiedN;
	
	/**
	 * Packing constraints
	 */
////	std::list<std::vector<int>> packing;
	std::vector<std::vector<int>> packing;
	
	std::vector<int> vRestore;

    int reductionSnapshotSize;
    std::vector<int> snapshotX;
	
	branch_and_reduce_algorithm(std::vector<std::vector<int>> &_adj, int const _N);

	int deg(int v);
	void set(int v, int a);

	fast_set used;

    // helpers for modifying the graph
	void compute_fold(std::vector<int> const &S, std::vector<int> const &NS);
	void compute_alternative(std::vector<int> const &A, std::vector<int> const &B);
	void restore(int n);
	void reverse();

    // helpers for lpReduction
	bool dinicDFS(int v);
	void updateLP();

    // reduction methods
	bool lpReduction();
	bool deg1Reduction();
	bool dominateReduction();
	bool fold2Reduction();
	bool twinReduction();
	bool funnelReduction();
	bool deskReduction();
	bool unconfinedReduction();
	int  packingReduction();

    // lower bounds for pruning
	int lpLowerBound();
	int cycleLowerBound();
	int cliqueLowerBound();
	int lowerBound();

    // recursive methods
	void branching(timer & t, double time_limit);
	bool decompose(timer & t, double time_limit);
	bool reduce();
	void rec(timer & t, double time_limit);

    // vestiges of original Java code
#if 0
	
	void debug(String str, Object...os) {
		StringBuilder sb = new StringBuilder();
		Calendar c = Calendar.getInstance();
		sb.append(String.format("%02d:%02d:%02d  ", c.get(Calendar.HOUR_OF_DAY), c.get(Calendar.MINUTE), c.get(Calendar.SECOND)));
		for (int i = 0; i < depth && i <= maxDepth; i++) sb.append(' ');
		System.err.print(sb);
		System.err.printf(str, os);
	}
#endif // 0
	
	int solve(timer & t, double time_limit);

    void initial_reduce_graph();
    void reduce_graph();

    void restore_to_snapshot();

    std::string debugString() const;
    void PrintState() const;

    size_t get_current_is_size() const;
    size_t get_current_is_size_with_folds() const;
    bool   folded_vertices_exist() const;
    std::vector<int> compute_maximal_is();
    size_t compute_alternative_maximal_is_size();
    size_t number_of_nodes_remaining() const;
    void force_into_independent_set(std::vector<NodeID> const &nodes);
    void extend_finer_is(std::vector<bool> &independent_set);
    void get_solved_is(std::vector<bool> &independent_set);

    void convert_adj_lists(graph_access & G, std::vector<NodeID> & reverse_mapping) const;

#if 0
}
#endif // 0
};

#endif //BRANCH_AND_REDUCE_SOLVER_H
