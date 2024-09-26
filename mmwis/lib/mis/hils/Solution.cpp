/*
 *
 * Solution.cpp
 *
 *  Created on: 16/09/2015
 *      Author: bruno
 */

#include "random_functions.h"
#include "Solution.h"

#include <algorithm>
#include <random>
#include <vector>

/* extern std::mt19937 number_generator; // Mersenne Twister 19937 generator */

Solution::Solution(graph_access *G) :
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

void Solution::moveFreeToSolutionPartition(const int v) {
	assert(v < G->number_of_nodes());

	// current position of v in the solution_ vector
	int pos_v = position_[v];

	// new position of v in the solution_ vector
	int new_pos_v = solution_size_;

	// first vertex of the second partition
	int j = solution_[solution_size_];

	// ensures v is in the free partition of the solution vector
	assert((solution_size_ <= pos_v) && (solution_size_ + free_size_ > pos_v));

	// swap v with the first vertex of the second partition
    std::swap(solution_[pos_v], solution_[new_pos_v]);
	position_[v] = new_pos_v;
	position_[j] = pos_v;

	// change the boundary between the blocks to make v the last vertex of the
	// first partition
	solution_size_++;
	free_size_--;
} // void Solution::moveFreeToSolutionPartition(const int v)

void Solution::moveFreeToNonFreePartition(const int v)
{
	assert(v < G->number_of_nodes());

	// current position of v in the solution vector
	int pos_v = position_[v];

	// new position of v in the solution vector
	int new_pos_v = solution_size_ + free_size_ - 1;

	// last vertex of the second partition
	int j = solution_[solution_size_ + free_size_ - 1];

	// ensures v is in the free partition of the solution vector
	assert((solution_size_ <= pos_v) && (solution_size_ + free_size_ > pos_v));

	// swap v with the last vertex of the second partition
	std::swap(solution_[pos_v], solution_[new_pos_v]);
	position_[v] = new_pos_v;
	position_[j] = pos_v;

	// change the boundary between the blocks to make v the last vertex of the
	// second partition
	free_size_--;
} // void Solution::moveFreeToNonFreePartition(const int v)

void Solution::moveSolutionToFreePartition(const int v)
{
	assert(v < G->number_of_nodes());

	// current position of v in the solution vector
	int pos_v = position_[v];

	// new position of v in the solution vector
	int new_pos_v = solution_size_ - 1;

	// last vertex of the first partition
	int j = solution_[solution_size_ - 1];

	// ensures v is in the solution partition of the solution vector
	assert(pos_v < solution_size_);

	// swap v with the last vertex of the second partition
    std::swap(solution_[pos_v], solution_[new_pos_v]);
	position_[v] = new_pos_v;
	position_[j] = pos_v;

	// change the boundary between the blocks to make v the first vertex of the
	// second partition
	solution_size_--;
	free_size_++;
} // void Solution::moveSolutionToFreePartition(const int v)

void Solution::moveNonFreeToFreePartition(const int v)
{
	assert(v < G->number_of_nodes());

	// current position of v in the solution vector
	int pos_v = position_[v];

	// new position of v in the solution vector
	int new_pos_v = solution_size_ + free_size_;

	// first vertex of the third partition
	int j = solution_[solution_size_ + free_size_];

	// ensures v is in the non free partition of the solution vector
	assert(pos_v >= solution_size_ + free_size_);

	// swap v with the last vertex of the second partition
    std::swap(solution_[pos_v], solution_[new_pos_v]);
	position_[v] = new_pos_v;
	position_[j] = pos_v;

	// change the boundary between the blocks to make v the last vertex of the
	// second partition
	free_size_++;
} // void Solution::moveNonFreeToFreePartition(const int v)

void Solution::addVertex(const int v)
{
	int weight_v = G->getNodeWeight(v);
	weight_ += weight_v;

	moveFreeToSolutionPartition(v);

	forall_out_edges((*G), e, v)
	    // increase the tighness of each neighbor by one
	    NodeID neighbor = G->getEdgeTarget(e);
		tightness_[neighbor]++;

		mu_[neighbor] -= weight_v;

		// if the neighbor is in the free partition, move to non free partition
		int neighbor_pos = position_[neighbor];
		if ((solution_size_ <= neighbor_pos) && (solution_size_ + free_size_ > neighbor_pos)) {
			moveFreeToNonFreePartition(neighbor);
		}
	endfor
} // void Solution::addVertex(const int v)

void Solution::removeVertex(const int v)
{
	int weight_v = G->getNodeWeight(v);
	weight_ -= weight_v;

	moveSolutionToFreePartition(v);

    forall_out_edges((*G), e, v)
        NodeID neighbor = G->getEdgeTarget(e);
		tightness_[neighbor]--;

		mu_[neighbor] += weight_v;

		// if the neighbor becomes free
		if (tightness_[neighbor] == 0) {
			moveNonFreeToFreePartition(neighbor);
		}
	endfor
} // void Solution::removeVertex(const int v)

bool Solution::integrityCheck() const
{
	for (int idx = 0; idx < solution_size_; idx++) {
		int vertex = solution_[idx];

		if (tightness_[vertex] > 0) {
			return false;
		}


        forall_out_edges((*G), e, vertex)
            NodeID neighbor = G->getEdgeTarget(e);
			if (find(solution_.begin(), solution_.begin() + solution_size_, neighbor)
			        != solution_.begin() + solution_size_) {
				return false;
			}
		endfor
	}

	for (int idx = solution_size_; idx < solution_size_ + free_size_; idx++) {
		int vertex = solution_[idx];
		if (tightness_[vertex] > 0) {
			return false;
		}
	}

	for (int idx = solution_size_ + free_size_; idx < G->number_of_nodes(); idx++) {
		int vertex = solution_[idx];
		if (tightness_[vertex] == 0) {
			return false;
		}
	}

	return true;
} // bool Solution::integrityCheck() const

void Solution::addRandomVertex()
{
	assert(!isMaximal());

	// generate a random number between [0, free_size_ - 1]
	int free_pos = random_functions::nextInt(0, free_size_ -1);

	int vertex = solution_[solution_size_ + free_pos];

	addVertex(vertex);
} // void Solution::addRandomVertex()

bool Solution::candOmegaImprovement(NodeID cand)
{
	int pos_cand = position_[cand];
    if (pos_cand<solution_size_) return false; //not in solution

	int v = solution_[pos_cand];
	if (mu_[v] > 0) {
        forall_out_edges((*G), e, v)
            NodeID neighbor = G->getEdgeTarget(e);
			if (position_[neighbor] < solution_size_) {
				removeVertex(neighbor);
			}
		endfor
		addVertex(v);
		return true;
    }

	return false;
} // bool Solution::candOmegaImprovement()



bool Solution::omegaImprovement()
{
	for (int idx = G->number_of_nodes() - 1; idx >= solution_size_; idx--) {
		int v = solution_[idx];
		if (mu_[v] > 0) {
            forall_out_edges((*G), e, v)
                NodeID neighbor = G->getEdgeTarget(e);
				if (position_[neighbor] < solution_size_) {
					removeVertex(neighbor);
				}
			endfor
			addVertex(v);
			return true;
		}
	}

	return false;
} // bool Solution::swapImprovement()

bool Solution::candTwoImprovement(NodeID cand)
{
	assert(isMaximal());
	int pos_cand = position_[cand];
    if (pos_cand>=solution_size_) return false; //not in solution

		// the candidate for removal
		int x = solution_[pos_cand];

		// sorted list of 1-tight nighbors of x
        std::vector<int> onetight_list;

		// build the list of 1-tight nighbors of x
        forall_out_edges((*G), e, x)
            NodeID neighbor = G->getEdgeTarget(e);
			if (tightness_[neighbor] == 1) {
				onetight_list.push_back(neighbor);
			}
		endfor
		assert(is_sorted(onetight_list.begin(), onetight_list.end()));

		// if x has fewer than two 1-tight neighors we are done with x
        if (onetight_list.size() < 2) return false;

		int x_weight = G->getNodeWeight(x);

		// attempt to find in onetight_list a pair {v, w} such that there
		// is no edge between v and w
		for (int v : onetight_list) {

			// stores the sorted list of v neighbors
			//vector<int> v_neighbors(g_->adj_l(v));
			//assert(is_sorted(v_neighbors.begin(), v_neighbors.end()));

			// check if there is a vertex w in onetight_list (besides v) that
			// does not belong to v_neighbors. since both onetight_list and v_neighbors
			// are sorted, this can be done by traversing both lists in tandem.
			size_t i_idx = G->get_first_edge(v), j_idx = 0;
			while (i_idx < G->get_first_invalid_edge(v)//v_neighbors.size()
			        && j_idx < onetight_list.size()) {
				if (onetight_list[j_idx] == v) {
					j_idx++;
					continue;
				} else if (G->getEdgeTarget(i_idx)/*v_neighbors[i_idx]*/ < onetight_list[j_idx]) {
					i_idx++;
					continue;
				}  else if (G->getEdgeTarget(i_idx)/*v_neighbors[i_idx]*/ == onetight_list[j_idx]) {
					i_idx++;
					j_idx++;
					continue;
				}

				// if this point is reached, this means we found the pair {v, w}
				// we were looking for !!
				int w = onetight_list[j_idx];

				int weight_v = G->getNodeWeight(v);
				int weight_w = G->getNodeWeight(w);

				if (x_weight >= weight_v + weight_w) {
					i_idx++;
					continue;
				}

				removeVertex(x);
				addVertex(v);
				addVertex(w);
				return true;
			}
		} // for(int v : onetight_list) {
	return false;
} // bool Solution::candTwoImprovment()


bool Solution::twoImprovement()
{
	assert(isMaximal());

	for (int idx = 0; idx < solution_size_; idx++) {
		// the candidate for removal
		int x = solution_[idx];

		// sorted list of 1-tight nighbors of x
        std::vector<int> onetight_list;

		// build the list of 1-tight nighbors of x
        forall_out_edges((*G), e, x)
            NodeID neighbor = G->getEdgeTarget(e);
			if (tightness_[neighbor] == 1) {
				onetight_list.push_back(neighbor);
			}
		endfor
		assert(is_sorted(onetight_list.begin(), onetight_list.end()));

		// if x has fewer than two 1-tight neighors we are done with x
		if (onetight_list.size() < 2) continue;

		int x_weight = G->getNodeWeight(x);

		// attempt to find in onetight_list a pair {v, w} such that there
		// is no edge between v and w
		for (int v : onetight_list) {

			// stores the sorted list of v neighbors
			//vector<int> v_neighbors(g_->adj_l(v));
			//assert(is_sorted(v_neighbors.begin(), v_neighbors.end()));

			// check if there is a vertex w in onetight_list (besides v) that
			// does not belong to v_neighbors. since both onetight_list and v_neighbors
			// are sorted, this can be done by traversing both lists in tandem.
			size_t i_idx = G->get_first_edge(v), j_idx = 0;
			while (i_idx < G->get_first_invalid_edge(v)//v_neighbors.size()
			        && j_idx < onetight_list.size()) {
				if (onetight_list[j_idx] == v) {
					j_idx++;
					continue;
				} else if (G->getEdgeTarget(i_idx)/*v_neighbors[i_idx]*/ < onetight_list[j_idx]) {
					i_idx++;
					continue;
				}  else if (G->getEdgeTarget(i_idx)/*v_neighbors[i_idx]*/ == onetight_list[j_idx]) {
					i_idx++;
					j_idx++;
					continue;
				}

				// if this point is reached, this means we found the pair {v, w}
				// we were looking for !!
				int w = onetight_list[j_idx];

				int weight_v = G->getNodeWeight(v);
				int weight_w = G->getNodeWeight(w);

				if (x_weight >= weight_v + weight_w) {
					i_idx++;
					continue;
				}

				removeVertex(x);
				addVertex(v);
				addVertex(w);
				return true;
			}
		} // for(int v : onetight_list) {
	} // for(int x : cadidate_list) {

	return false;
} // bool Solution::twoImprovment()

void Solution::force(int k)
{
	for(int i = 0; i < k; i++) {
		// select a non solution vertex to add
		int nonsolution_size = G->number_of_nodes() - (solution_size_ + free_size_);
		int nonsolution_pos = random_functions::nextInt(0, nonsolution_size -1);
		int vertex = solution_[solution_size_ + free_size_ + nonsolution_pos];

		// remove the neighboring vertices as necessary

        forall_out_edges((*G), e, vertex)
            NodeID neighbor = G->getEdgeTarget(e);
			if (position_[neighbor] < solution_size_) {
				removeVertex(neighbor);
			}
		endfor
		addVertex(vertex);
	}
} // void Solution::force()

void Solution::force_candidate(NodeID candidate)
{
	// remove the neighboring vertices as necessary

    forall_out_edges((*G), e, candidate)
        NodeID neighbor = G->getEdgeTarget(e);
		if (position_[neighbor] < solution_size_) {
			removeVertex(neighbor);
		}
	endfor
	addVertex(candidate);
} // void Solution::force_candidate()

