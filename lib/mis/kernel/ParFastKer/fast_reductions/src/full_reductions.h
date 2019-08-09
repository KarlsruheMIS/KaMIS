 /******************************************************************************
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
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

#ifndef FULL_REDUCTIONS_H
#define FULL_REDUCTIONS_H

#include <vector>
#include <memory>
#include "mis_config.h"
#include "parallel_reductions.h"
#include "mis/kernel/ParFastKer/LinearTime/MIS_sigmod_pub/Graph.h"


class full_reductions
{
public:
	full_reductions(std::vector<std::vector<int>> &_adj, int _N);
	void reduce_graph();
	size_t get_current_is_size_with_folds();
	size_t number_of_nodes_remaining();
    void force_into_independent_set(std::vector<NodeID> &nodes_to_force);
    void convert_adj_lists(graph_access &G, std::vector<NodeID> &reverse_mapping);
    void extend_finer_is(std::vector<bool> &independent_set);

	static std::vector<int> partitions;
private:
	std::vector<std::vector<int>> &adj;
	std::unique_ptr<parallel_reductions> parallel_reducer;
    std::vector<unsigned int> forced_vertices;
    int linearTimeOffset;
    LinearTime::Graph LineartimeKernelizer;
    std::vector<int> after_force_mapping;
    std::vector<int> after_force_reverse_mapping;
};


#endif // FULL_REDUCTIONS_H
