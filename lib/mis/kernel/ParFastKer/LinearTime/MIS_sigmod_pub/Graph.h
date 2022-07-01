/******************************************************************************
 * Copyright (C) 2019 Lijun Chang <ljchang.au@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *****************************************************************************/

#ifndef _GRAPH_H_
#define _GRAPH_H_

// #include "Utility.h"
#include "mis/kernel/ParFastKer/LinearTime/MIS_sigmod_pub/Utility.h"

namespace LinearTime {
struct Edge {
	int id, duplicate;
	int next;
};

class Graph {
private:
	std::string dir; //input graph directory
	ui n, m; //number of nodes and edges of the graph

	ui *pstart; //offset of neighbors of nodes
	ui *edges; //adjacent ids of edges

    char *is;
    char *fixed;
    ui *reverse_mapping;
	std::vector<std::pair<ui,ui> > S;
public:
	Graph() ;
	~Graph() ;

	void read_graph() ;
    void read_graph_metis() ;
    void inputGraphAdjList(std::vector<std::vector<int>> &adj);
    void inputEdgeList(std::vector<ui> &_pstart, std::vector<ui> &_edges);
    void write_graph_metis(const char *is, const int *degree, const ui *pend);
	void degree_one_kernal_and_remove_max_degree() ;
	void degree_two_kernal_and_remove_max_degree_with_contraction() ;
	int degree_two_kernal_and_remove_max_degree_without_contraction(std::vector<std::vector<int>> &out_kernel) ;
	void degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction() ;

    void UndoReductions(std::vector<bool> &in_out_is);

	void greedy() ;
	void greedy_dynamic() ;

private:
    void write_kernel_to_adj_list(const char *is, const int *degree, const ui *pend, std::vector<std::vector<int>> &out_kernel) ;
	int general_swap(char *is, char *fixed = NULL) ;
	void check_is(const char *is, int count) ;
	void compute_upperbound(const char *is, char *fixed = NULL) ;

	void get_two_neighbors(ui u, char *is, ui &u1, ui &u2) ;
	ui get_other_neighbor(ui u, char *is, ui u1) ;
	int exist_edge(ui u1, ui u2) ;
	void edge_rewire(ui u, ui u1, ui u2) ;

	int exist_edge(ui u, ui v, const ui *pend) ;
	int find_other_endpoint(ui u, ui v, char *is) ;
	ui edge_rewire(ui u, const ui *pend, ui v, ui w) ;

	int remove_degree_one_two(std::vector<ui> &degree_ones, std::vector<ui> &degree_twos, char *is, int *degree, std::vector<std::pair<ui,ui> > &S) ;

	int initial_dominance_and_degree_two_remove(std::vector<ui> &degree_ones, std::vector<ui> &degree_twos, char *is, int *degree, char *adj, std::vector<std::pair<ui,ui> > &S) ;

	int lp_reduction(ui *ids, ui ids_n, char *is, int *degree) ;

	void shrink(ui u, ui &end, const char *is) ;
	void shrink(ui u, ui &end, const char *is, ui *tri) ;
	void update_triangle(ui u1, ui u2, ui *pend, char *is, char *adj, ui *tri, int *degree, char *dominate, std::vector<ui> &dominated) ;
	int dominated_check(ui u, ui *pend, char *is, ui *tri, int *degree) ;
	int compute_triangle_counts(ui *tri, ui *pend, char *adj, char *is, int *degree, char *dominate, std::vector<ui> &dominated) ;
	void construct_degree_increase(ui *ids) ;

	int delete_vertex(ui v, char *is, int *degree, std::vector<ui> &degree_ones) ;
	int delete_vertex(ui v, char *is, int *degree, std::vector<ui> &degree_ones, std::vector<ui> &degree_twos) ;
	int delete_vertex(ui v, const ui *pend, char *is, int *degree, std::vector<ui> &degree_ones, std::vector<ui> &degree_twos) ;
	int delete_vertex(ui u, ui *pend, char *is, std::vector<ui> &degree_twos, ui *tri, char *adj, int *degree, char *dominate, std::vector<ui> &dominated) ;
	int delete_vertex(ui v, char *is, int *degree, int *head, Edge *es, int *bin_head, int *bin_next, int *bin_pre, std::vector<ui> &degree_ones, std::vector<ui> &degree_twos) ;
};
}
#endif
