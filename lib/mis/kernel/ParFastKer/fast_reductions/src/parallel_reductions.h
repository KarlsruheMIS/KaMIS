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

#ifndef PARALLEL_REDUCTIONS_H
#define PARALLEL_REDUCTIONS_H

// #include "Set.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "Reduction.h"
#include "SimpleSet.h"
#include "fast_set.h"
#include "MaximumMatching.h"
#include "data_structure/graph_access.h"

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <ctime>
#include <string>
#include <atomic>
#include <mutex>

#define TIMERS

class parallel_reductions
{
public:
    parallel_reductions(std::vector<std::vector<int>> const &adjacencyArray, std::vector<int> const &vertexPartitions);
    ~parallel_reductions();

    void reduce_graph_parallel(std::vector<unsigned int> forced_vertices);
    void reduce_graph_sequential();
    void reduce_graph_sequential_reduction_wise();

    void ApplyReductions(int const partition, std::vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining, std::vector<int> &tempInt1, std::vector<int> &tempInt2, fast_set &fastSet, std::vector<int> &tempIntDoubleSize, double &time, int &isolatedCliqueCount, int &foldedVertexCount, int &removedTwinCount, int &foldedTwinCount, int &removedUnconfinedVerticesCount, int &numDiamondReductions);
    void UndoReductions(std::vector<Reduction> const &vReductions);
    std::vector<std::vector<int>> getKernel();
    void applyKernelSolution(std::vector<int> kernel_solution);
    void ApplyKernelSolutionToReductions(std::vector<Reduction> const &vReductions);

    std::vector<int> independent_set;

    size_t size() const { return inGraph.Size(); }

    std::vector<SparseArraySet> const& Neighbors()  const { return neighbors;  }

    void SetAllowVertexFolds(bool const allow) { m_bAllowVertexFolds = allow; }


    void force_into_independent_set(std::vector<unsigned int> &nodes_to_force);
	size_t get_current_is_size_with_folds();
    void getGraphAccess(graph_access &G, std::vector<unsigned int> &reverse_mapping);
    void ExtendPartialSolution(std::vector<bool> &independent_set);

protected: // methods
    bool removeUnconfined(int const partition, int const vertex, ArraySet &remaining, fast_set &closedNeighborhood, std::vector<int> &neighborhood, std::vector<int> &numNeighborsInS, std::vector<int> &neighborsInS, int &removedUnconfinedVerticesCount, int &numDiamondReductions);
    bool removeTwin(int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, std::vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount);
    bool RemoveIsolatedClique    (int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, std::vector<bool> &vMarkedVertices, int &isolatedCliqueCount);
    bool FoldVertex(int const partition, int const vertex, std::vector<Reduction> &vReductions, ArraySet &remaining, int &foldedVertexCount);
    bool LPReduction(std::vector<ArraySet> &remainingPerPartition, std::vector<std::vector<int>> &bufferPerPartition, int &numLPReductions);
    void initReducableVertices(int numPartitions);
    bool isTwoNeighborhoodInSamePartition(int const vertex, int const partition, ArraySet &remaining);
    void UpdateRemaining(std::vector<ArraySet> &remainingPerPartition, std::vector<std::vector<int>> &bufferPerPartition);

    bool removeAllUnconfined(int const partition, ArraySet *remainingInsert, fast_set &closedNeighborhood, std::vector<int> &neighborhood, std::vector<int> &numNeighborsInS, std::vector<int> &neighborsInS, int &removedUnconfinedVerticesCount, int &numDiamondReductions);
    bool removeAllTwin(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, std::vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount);
    bool RemoveAllIsolatedClique(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, std::vector<bool> &vMarkedVertices, int &isolatedCliqueCount);
    bool FoldAllVertices(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, int &foldedVertexCount);

    int degree(int const vertex);
    bool isBoundaryVertex(const int vertex);

    void updateDependencyCheckingEstimation(int partition);
    void initDependencyCheckingEstimation(int partition);
    bool shouldStopDependencyCheckingReductions(int partition);

    void initGlobalBurstEstimator();
    double finishThreadAndGetEstimatedBurstLength(int tid);
    bool isLastFinishedThread(int tid);
    bool shouldTerminate();
    void terminateOtherThreads();

    // Just for testing
    bool checkDegrees();

protected: // members
    std::vector<int> graph_to_kernel_map;
    std::vector<int> kernel_solution;
    std::vector<std::vector<std::vector<Reduction>>> AllReductions;
    std::vector<std::vector<int>> const m_AdjacencyArray;
    std::vector<SparseArraySet>     neighbors;
    SimpleSet inGraph;
    SimpleSet neighborhoodChanged;
    std::vector<int> partitions;
    std::vector<std::vector<int>> partition_nodes;
    std::vector<ArraySet> inGraphPerPartition;
    MaximumMatching maximumMatching;
    std::vector<std::atomic_int> vertexDegree;
    std::vector<std::atomic_int> numCutEdges;
    std::vector<double> dependecy_checking_burst_estimation;
    std::vector<double> dependency_checking_times;
    double global_burst_timer;
    double global_burst_estimation;
    std::mutex global_burst_timer_mutex;
    int lastFinishedThread;
    bool terminationFlag;
    int firstFinished;
    int is_offset;
#ifdef TIMERS
    clock_t replaceTimer;
    #endif // TIMERS
    bool m_bAllowVertexFolds;
};

#endif //PARALLEL_REDUCTIONS_H
