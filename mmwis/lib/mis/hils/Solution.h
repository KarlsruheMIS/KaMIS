/*
 *
 * Solution: solution data structre
 *
 *  Created on: 16/09/2015
 *      Author: bruno
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_
#include <assert.h>
#include <vector>
#include <string>
#include "mmwis_graph_access.h"


class Solution
{

public:

	template <typename graph>
	Solution(graph *G);

	// add a vertex to the solution

	void addVertex(const int v);

	// remove a vertex from the solution

	void removeVertex(const int v);

	// randomly add a free vertex to the solution

	void addRandomVertex();

	// find a (1,2)-swap improvment

	bool twoImprovement();
	bool candTwoImprovement(NodeID candidate);

	// find a (omega,1)-swap improvment

	bool omegaImprovement();
	bool candOmegaImprovement(NodeID candidate);

	// randomly insert k vertices into the solution

	void force(int k);

	void force_candidate(NodeID cand);

	// check if the solution is maximal

	bool isMaximal()
	{
		return free_size_ == 0;
	}

	// perform integrity check on the solution. returns true if
	// the solution is ok, and false, otherwise. for testing purposes
	// only.

	bool integrityCheck() const;

	// return the current solution weight

	int weight() const
	{
		return weight_;
	}

	// return the current size of the independent set
	
	int size() const
	{
		return solution_size_;
	}	

	// return a vector indicating the state of each vertex
	
	std::vector<int> solution() const
	{
		std::vector<int> sol;
		for (int idx = 0; idx < G->number_of_nodes(); idx++) {
			if(position_[idx] < solution_size_) {
				sol.push_back(1);
			} else {
				sol.push_back(0);
			}
		}
		return sol;
	}

	// return a vector with the vertices in the solution

	std::vector<int> i_set() const
	{
		std::vector<int> iset;
		for (int idx = 0; idx < G->number_of_nodes(); idx++) {
			if(position_[idx] < solution_size_) 
				iset.push_back(idx);
		}
		return iset;
	}


private:

	// problem instance

    mmwis::graph_access *G;

	// the solution_ vector is partitioned into three blocks: first vertices in the solution, then 
	// the free vertices (i.e., vertices that are not adjacent to any vertex in the solution), and 
	// finally the non-solution vertices that are not free

	std::vector<int> solution_;

	// size of the solution verticies partition

	int solution_size_;

	// size of the free vertices partition

	int free_size_;

	// for each vertex, the number of adjacent vertices that are on the solution

	std::vector<int> tightness_;

	// position of each vertex in the solution_ vector

	std::vector<int> position_;

	// weight of each vertex i minus the sum of the weights of its neighbors that
	// are in the independent set

	std::vector<int> mu_;
	
	// current independent vertex weight
	
	int weight_;

	// move a vertex from the free partition to solution partition

	void moveFreeToSolutionPartition(const int v);

	// move a vertex from the free patition to non free partition

	void moveFreeToNonFreePartition(const int v);

	// move a vertex from the solution partition to free partition

	void moveSolutionToFreePartition(const int v);

	// move a vertex from the non free partition to free partition

	void moveNonFreeToFreePartition(const int v);

}; // class Solution

template <typename graph>
inline Solution::Solution(graph *G) :
	G(G),
	solution_(G->number_of_nodes()),
	solution_size_(0),
	free_size_(G->number_of_nodes()),
	tightness_(G->number_of_nodes(), 0),
	position_(G->number_of_nodes()),
	mu_(G->number_of_nodes()),
	weight_(0) {
	for (int idx = 0; idx < G->number_of_nodes(); idx++) {
		position_[idx] = idx;
		solution_[idx] = idx;
		mu_[idx] = G->getNodeWeight(idx);
	}
} // Solution::Solution(const Graph *g)

#endif // #ifndef SOLUTION_H_
