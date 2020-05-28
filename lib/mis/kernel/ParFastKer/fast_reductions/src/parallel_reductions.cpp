 /******************************************************************************
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
 *
  *****************************************************************************/

#include "parallel_reductions.h"
#include "ArraySet.h"
#include "SparseArraySet.h"
#include "ProfilingHelper.h"

#include <vector>
#include <set>
#include <iostream>
#include <ctime>
#include <numeric>
#include <cassert>
#include <climits>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "omp.h"
#include <assert.h>
#include <limits.h>
#include <functional>

#define ISOLATED_CLIQUE_MAX_NEIGHBORS 2
#define MAX_SIZE_UNCONFINED 6
#define DEPENDENCY_CHECKING_THRESHOLD_MULTIPLIER 3.0
#define DEPENDENCYCHECKING_BURST_ESTIMATION_ALPHA 0.5
#define GLOBAL_BURST_ESTIMATION_ALPHA 0.5
#define GLOBAL_BURST_THRESHOLD_MULTIPLIER 6

#define INSERT_REMAINING(partition, remaining, v) if(partitions[v] == partition) remaining.Insert(v);
// Remove vertex from inGraph first!
#define REMOVE_NEIGHBOR(partition, neighbor, vertex) {assert(!inGraph.Contains(vertex)); assert(partitions[vertex] == partition); vertexDegree[neighbor]--; if(partition != partitions[neighbor]) {numCutEdges[neighbor]--; neighborhoodChanged.Insert(neighbor);}}
#define REMOVE_VERTEX(partition, vertex) {inGraph.Remove(vertex); inGraphPerPartition[partition].Remove(vertex);}

using namespace std;

ProfilingHelper_t profilingHelper;

parallel_reductions::parallel_reductions(vector<vector<int>> const &adjacencyArray, vector<int> const &vertexPartitions)
 : m_AdjacencyArray(adjacencyArray)
 , neighbors(adjacencyArray.size())
 , inGraph(adjacencyArray.size(), true)
 , neighborhoodChanged(adjacencyArray.size(), false)
 , partitions(vertexPartitions)
 , independent_set(adjacencyArray.size(), -1)
 , maximumMatching(adjacencyArray)
 , vertexDegree(adjacencyArray.size())
 , numCutEdges(adjacencyArray.size())
#ifdef TIMERS
 , replaceTimer(0)
 #endif // TIMERS
 , m_bAllowVertexFolds(true)
{
    int N = adjacencyArray.size();
    for (size_t u=0; u < (size_t)N; ++u) {
        neighbors[u].InitializeFromAdjacencyArray(m_AdjacencyArray, u);
    }
    auto numPartitions_ptr = max_element(partitions.begin(), partitions.end());
    int numPartitions = numPartitions_ptr!=partitions.end() ? *numPartitions_ptr + 1 : 0;
    if(numPartitions == 1) {
        numPartitions = 2;
    }
    partition_nodes = std::vector<std::vector<int>>(numPartitions);
    for(int node = 0; node < N; ++node) {
        assert(partitions[node] >= 0);
        assert(partitions[node] < numPartitions);
        partition_nodes[partitions[node]].push_back(node);
        numCutEdges[node] = 0;
    }

    for(int node = 0; node < N; ++node) {
        vertexDegree[node] = neighbors[node].Size();
        for(auto neighbor: neighbors[node]) {
            if(partitions[neighbor] != partitions[node]) {
                numCutEdges[node]++;
            }
        }
    }
}

parallel_reductions::~parallel_reductions()
{

#ifdef TIMERS
    //cout << "Total time spent undoing  reductions  : " << (replaceTimer/(double)CLOCKS_PER_SEC) << endl;
#endif // TIMERS
}

std::vector<std::vector<int>> parallel_reductions::getKernel() {
  double startclock = omp_get_wtime();
    graph_to_kernel_map = std::vector<int> (m_AdjacencyArray.size());
    int nodecount = 0;
    for(int node = 0; node < (int)m_AdjacencyArray.size(); node++) {
        if(inGraph.Contains(node)) {
            assert(independent_set[node] == -1);
            graph_to_kernel_map[node] = nodecount++;
        }
    }

    std::vector<std::vector<int>> kernel_adj(nodecount);

    // Build adjacency vectors
    #pragma omp parallel for
    for(int node = 0; node < (int)m_AdjacencyArray.size(); node++) {
        if(inGraph.Contains(node)) {
            kernel_adj[graph_to_kernel_map[node]].reserve(neighbors[node].Size());
            for(auto neighbor : neighbors[node]) {
                if(inGraph.Contains(neighbor)) {
                    kernel_adj[graph_to_kernel_map[node]].push_back(graph_to_kernel_map[neighbor]);
                }
            }
            std::sort(kernel_adj[graph_to_kernel_map[node]].begin(), kernel_adj[graph_to_kernel_map[node]].end());
        }
    }
    double endclock = omp_get_wtime();
    std::cout << "getKernel took " << (endclock - startclock) << std::endl;
    return kernel_adj;
}

void parallel_reductions::applyKernelSolution(std::vector<int> kernel_solution){
    for(int node = 0; node < (int)m_AdjacencyArray.size(); ++node) {
        if(inGraph.Contains(node)) {
            independent_set[node] = kernel_solution[graph_to_kernel_map[node]];
        }
    }
    for(int i = AllReductions.size(); i > 0; i--) {
        for(auto Reductions: AllReductions[i - 1]) {
            ApplyKernelSolutionToReductions(Reductions);
        }
    }
}

int parallel_reductions::degree(int const vertex) {
    return vertexDegree[vertex];
}

bool parallel_reductions::isBoundaryVertex(const int vertex) {
    return numCutEdges[vertex] > 0;
}

bool parallel_reductions::RemoveIsolatedClique(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &isolatedCliqueCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    if(isBoundaryVertex(vertex)) {
        profilingAddTimeUnsuccessfulIsolatedCliquePartition(&profilingHelper, partition);
        return false;
    }

    int const deg = degree(vertex);

    if(deg > ISOLATED_CLIQUE_MAX_NEIGHBORS)
        return false;

    for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        assert(partitions[neighbor] == partition);
        if (degree(neighbor) < deg) {
            profilingAddTimeUnsuccessfulIsolatedCliqueDegree(&profilingHelper, partition);
            return false;
        }
    }

    bool superSet(true);

    for (int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        // TODO Should be possible faster
        for (int const nNeighbor : neighbors[neighbor]) if(inGraph.Contains(nNeighbor)) {
            vMarkedVertices[nNeighbor] = true;
        }
        vMarkedVertices[neighbor] = true;

        for (int const neighbor2 : neighbors[vertex]) if(inGraph.Contains(neighbor2)) {
            superSet = superSet && vMarkedVertices[neighbor2];
        }

        for (int const nNeighbor : neighbors[neighbor]) if(vMarkedVertices[nNeighbor]) {
            vMarkedVertices[nNeighbor] = false;
        }
        vMarkedVertices[neighbor] = false;

        if (!superSet) {
            profilingAddTimeUnsuccessfulIsolatedCliqueNoClique(&profilingHelper, partition);
            return false;
        }
    }
    if (superSet) {
        // Reduction reduction(ISOLATED_VERTEX);
        // reduction.SetVertex(vertex);
        independent_set[vertex] = 0;
        REMOVE_VERTEX(partition, vertex);
        for (const int neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            assert(partitions[neighbor] == partition);
            REMOVE_VERTEX(partition, neighbor);
            remaining.Remove(neighbor);
            independent_set[neighbor] = 1;
            for (const int nNeighbor : neighbors[neighbor]) if(inGraph.Contains(nNeighbor)) {
                // assert(neighbors[nNeighbor].Contains(neighbor));
                REMOVE_NEIGHBOR(partition, nNeighbor, neighbor);
                INSERT_REMAINING(partition, remaining, nNeighbor);
            }
            neighbors[neighbor].Clear();
            assert(!inGraph.Contains(neighbor));
        }
        isolatedCliqueCount += deg + 1;
        neighbors[vertex].Clear();

        // vReductions.emplace_back(std::move(reduction));

        profilingAddTimeSuccessfulIsolatedClique(&profilingHelper, partition);
        return true;
    }
    assert(false);
    return false;
}

bool parallel_reductions::isTwoNeighborhoodInSamePartition(int const vertex, int const partition, ArraySet &remaining) {
    if(isBoundaryVertex(vertex)) {
        return false;
    }
    for(int neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        if(isBoundaryVertex(neighbor)) {
            return false;
        } else if(neighborhoodChanged.Contains(neighbor)) {
            neighborhoodChanged.Remove(neighbor);
            remaining.Insert(neighbor);
        }
    }
    return true;
}

bool parallel_reductions::removeUnconfined(int const partition, int const vertex, ArraySet &remaining, fast_set &closedNeighborhood, vector<int> &neighborhood, vector<int> &numNeighborsInS, vector<int> &neighborsInS, int &removedUnconfinedVerticesCount, int &numDiamondReductions) {
    assert(neighborhood.size() >= neighbors.size());
    assert(numNeighborsInS.size() >= neighbors.size());
    assert(neighborsInS.size() >= 2 * neighbors.size());
    closedNeighborhood.clear();
    closedNeighborhood.add(vertex);
    int sizeS = 1, sizeNeighborhood = 0;
    for (int u : neighbors[vertex]) if(inGraph.Contains(u)) {
        closedNeighborhood.add(u);
        if(partitions[u] == partition) {
            neighborhood[sizeNeighborhood++] = u;
            numNeighborsInS[u] = 1;
        }
    }
    bool vertexAddedToS = true;

    while (vertexAddedToS) {
        vertexAddedToS = false;
        for (int i = 0; i < sizeNeighborhood; i++) {
            int const u = neighborhood[i];
            if (numNeighborsInS[u] != 1)  {
                continue;
            } 
            int neighborToAdd = -1;
            for (int const w : neighbors[u]) if(inGraph.Contains(w) && !closedNeighborhood.get(w)) {
                if (neighborToAdd >= 0) {
                    // There is more than 1 neighbor outside of N[S]
                    neighborToAdd = -2;
                    break;
                }
                neighborToAdd = w;
            }
            if (neighborToAdd == -1) {
                // There is a vertex in N(u) that doesn't have any neighbors outside of N[S]
                // Input vertex is unconfined
                independent_set[vertex] = 1;
                inGraph.Remove(vertex);
                for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                    REMOVE_NEIGHBOR(partition, neighbor, vertex);
                    INSERT_REMAINING(partition, remaining, neighbor);
                }
                neighbors[vertex].Clear();
                remaining.Remove(vertex);
                ++removedUnconfinedVerticesCount;
                return true;
            } else if (neighborToAdd >= 0) {
                // if(sizeS >= MAX_SIZE_UNCONFINED) {
                //     continue;
                // }
                // there is a vertex in N(u) that has exactly one neighbor outside of N[S]
                // that vertex has to be added to S
                if(partitions[neighborToAdd] == partition) {
                    vertexAddedToS = true;
                    closedNeighborhood.add(neighborToAdd);
                    sizeS++;
                    for (int w : neighbors[neighborToAdd]) if(inGraph.Contains(w)) {
                        if (closedNeighborhood.add(w)) {
                            if(partitions[w] == partition) {
                                neighborhood[sizeNeighborhood++] = w;
                                numNeighborsInS[w] = 1;
                            }
                        } else {
                            if(partitions[w] == partition)
                                numNeighborsInS[w]++;
                        }
                    }
                }
            }
        }
    }
    
    if (sizeS >= 2) {
        closedNeighborhood.clear();
        int N = neighbors.size();
        for (int i = 0; i < sizeNeighborhood; i++) closedNeighborhood.add(neighborhood[i]);
        for (int i = 0; i < sizeNeighborhood; i++) {
            neighborsInS[i] = neighborsInS[N + i] = -1;
            int u = neighborhood[i];
            if (numNeighborsInS[u] != 2) continue;
            int v1 = -1, v2 = -1;
            // numNeighborsInS[u] == 2 assures that there are exactly two neighbors in S
            // !closedNeighborhood.get(w) can only cause the loop to find more vertices, not less
            // => only vertices with exactly 2 neighbors outside of N(S) are found
            for (int w : neighbors[u]) if (inGraph.Contains(w) && !closedNeighborhood.get(w)) {
                if (v1 < 0) v1 = w;
                else if (v2 < 0) v2 = w;
                else {
                    v1 = v2 = -1;
                    break;
                }
            }
            if (v1 > v2) {
                int t = v1;
                v1 = v2;
                v2 = t;
            }
            neighborsInS[i] = v1;
            neighborsInS[N + i] = v2;
        }
        for (int i = 0; i < sizeNeighborhood; i++) if (neighborsInS[i] >= 0 && neighborsInS[N + i] >= 0) {
            int u = neighborhood[i];
            closedNeighborhood.clear();
            // TODO
            for (int w : neighbors[u]) if(inGraph.Contains(w)) closedNeighborhood.add(w);
            for (int j = i + 1; j < sizeNeighborhood; j++) if (neighborsInS[i] == neighborsInS[j] && neighborsInS[N + i] == neighborsInS[N + j] && !closedNeighborhood.get(neighborhood[j])) {
                // Vertex is unconfined
                independent_set[vertex] = 1;
                inGraph.Remove(vertex);
                for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                    REMOVE_NEIGHBOR(partition, neighbor, vertex);
                    INSERT_REMAINING(partition, remaining, neighbor);              
                }
                neighbors[vertex].Clear();
                remaining.Remove(vertex);
                ++numDiamondReductions;
                return true;
            }
        }
    }
    return false;
}

bool parallel_reductions::removeTwin(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount)
{
    assert(partitions[vertex] == partition);
    assert(vMarkedVertices.size() == neighbors.size());
    // This takes really long (it's O(n))
    // assert(std::accumulate(vMarkedVertices.begin(), vMarkedVertices.end(), false, std::logical_or<bool>()) == false);
    if(isBoundaryVertex(vertex))
        return false;

    if(degree(vertex) != 3)
        return false;

    int twinNeighbors[3];
    int numNeighbors = 0;


    for(int const neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        assert(numNeighbors < 3);
        twinNeighbors[numNeighbors++] = neighbor;
    }


    int smallestDegreeNeighbor = -1;
    int smallesDegreeNeighborDegree = INT_MAX;
    for(int i = 0; i < 3; ++i) {
        int const neighbor = twinNeighbors[i];
        assert(partitions[neighbor] == partition);
        assert(neighbor != vertex);

        vMarkedVertices[neighbor] = true;
        int const neighborDegree = degree(neighbor);
        if(neighborDegree < smallesDegreeNeighborDegree) {
            smallesDegreeNeighborDegree = neighborDegree;
            smallestDegreeNeighbor = neighbor;
        }
    }
    assert(smallestDegreeNeighbor != -1);

    int twin = -1;
    for(int possibleTwin: neighbors[smallestDegreeNeighbor]) if(inGraph.Contains(possibleTwin)) {
        if(possibleTwin == vertex) continue;
        if(partitions[possibleTwin] != partition) continue;
        if(vMarkedVertices[possibleTwin]) continue;
        if(degree(possibleTwin) != 3) continue;
        assert(partitions[possibleTwin] == partitions[vertex]);
        bool isTwin = true;
        int neighborCount(0);
        for(int twinNeighbor: neighbors[possibleTwin]) if(inGraph.Contains(twinNeighbor)) {
            if(!vMarkedVertices[twinNeighbor]) {
                isTwin = false;
                break;
            }
            ++neighborCount;
        }
        if(isTwin && neighborCount == 3) {
            twin = possibleTwin;
            break;
        }
    }
    if(twin == -1) {
        for(int i = 0; i < 3; ++i) {
            int const markedVertex = twinNeighbors[i];
            vMarkedVertices[markedVertex] = false;
        }
        return false;
    }
    assert(twin >= 0);
    assert(partitions[twin] == partitions[vertex]);
    assert(!isBoundaryVertex(twin));
    assert(!neighbors[vertex].Contains(twin));

    bool isNeighborhoodAdjacent = false;
    for(int i = 0; i < 3; ++i) {
        int const neighbor1 = twinNeighbors[i];
        for(int neighbor2: neighbors[neighbor1]) {
            if(vMarkedVertices[neighbor2]) {
                isNeighborhoodAdjacent = true;
                goto afterNeighborhoodCheck;
            }
        }
    }
afterNeighborhoodCheck:

    bool reduced = false;
    if(isNeighborhoodAdjacent) {
        // Case where all vertices get removed from the graph and the twins get inserted into the independent set
        REMOVE_VERTEX(partition, vertex);
        independent_set[vertex] = 0;
        remaining.Remove(vertex);
        REMOVE_VERTEX(partition, twin);
        independent_set[twin] = 0;
        remaining.Remove(twin);
        for(int i = 0; i < 3; ++i) {
            int const neighbor1 = twinNeighbors[i];
            assert(partitions[neighbor1] == partition);
            REMOVE_VERTEX(partition, neighbor1);
            remaining.Remove(neighbor1);
            independent_set[neighbor1] = 1;
            for(int const neighbor2: neighbors[neighbor1]) if(inGraph.Contains(neighbor2)) {
                REMOVE_NEIGHBOR(partition, neighbor2, neighbor1);
                INSERT_REMAINING(partition, remaining, neighbor2);
            }
            neighbors[neighbor1].Clear();
        }
        neighbors[vertex].Clear();
        neighbors[twin].Clear();
        removedTwinCount += 5;
        reduced = true;

    } else {
        // Case where the vertices get folded
        bool twoNeighborHoodInSamePartition = isTwoNeighborhoodInSamePartition(vertex, partition, remaining);
        if(!twoNeighborHoodInSamePartition) {
            reduced = false;
        } else {
            Reduction reduction(FOLDED_TWINS);
            reduction.SetVertex(vertex);
            reduction.SetTwin(twin);
            for(int i = 0; i < 3; ++i) {
                int const neighbor = twinNeighbors[i];
                reduction.AddNeighbor(neighbor);
            }

            int neighborHoodSize(0);
            for(int neighbor: reduction.GetNeighbors()) {
                assert(partitions[neighbor] == partitions[twin]);
                // TODO
                neighbors[neighbor].Remove(twin);
                // TODO
                neighbors[neighbor].Remove(vertex);
                neighborHoodSize += degree(neighbor);
            }
            neighbors[twin].Clear();
            neighbors[vertex].Clear();
            neighbors[vertex].Resize(neighborHoodSize);
            for(int neighbor1: reduction.GetNeighbors()) {
                assert(!isBoundaryVertex(neighbor1));
                for(int neighbor2: neighbors[neighbor1]) if(inGraph.Contains(neighbor2)) {
                    assert(partitions[neighbor2] == partitions[neighbor1]);
                    assert(neighbor2 != vertex);
                    assert(neighbor2 != twin);
                    assert(!vMarkedVertices[neighbor2]);
                    // TODO
                    neighbors[neighbor2].Remove(neighbor1);
                    vertexDegree[neighbor2]--;
                    neighbors[vertex].Insert(neighbor2);
                    if(!neighbors[neighbor2].Contains(vertex))
                        vertexDegree[neighbor2]++;
                    neighbors[neighbor2].Insert(vertex);
                    remaining.Insert(neighbor2);
                }
                neighbors[neighbor1].Clear();
                REMOVE_VERTEX(partition, neighbor1);
                remaining.Remove(neighbor1);
            }
            vertexDegree[vertex] = neighbors[vertex].Size();
            REMOVE_VERTEX(partition, twin);
            assert(!isBoundaryVertex(twin));
            remaining.Remove(twin);
            remaining.Insert(vertex);
            vReductions.push_back(reduction);
            assert(inGraph.Contains(vertex));
            assert(!inGraph.Contains(reduction.GetNeighbors()[0]));
            assert(!inGraph.Contains(reduction.GetNeighbors()[1]));
            assert(!inGraph.Contains(reduction.GetNeighbors()[2]));
            foldedTwinCount += 4;
            reduced = true;
        }
    }

    for(int i = 0; i < 3; ++i) {
        int const markedVertex = twinNeighbors[i];
        vMarkedVertices[markedVertex] = false;
    }
    return reduced;
}

bool parallel_reductions::FoldVertex(int const partition, int const vertex, vector<Reduction> &vReductions, ArraySet &remaining, int &foldedVertexCount)
{
    assert(partitions[vertex] == partition);

    profilingStartClock(&profilingHelper, partition, vertex);

    if(degree(vertex) != 2) {
        profilingAddTimeUnsuccessfulFoldDegree(&profilingHelper, partition);
        return false;
    }

    if(isBoundaryVertex(vertex)) {
        return false;
    }

    int neighbor1 = -1;
    int neighbor2 = -1;
    for(int const neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
        if(neighbor1 == -1)
            neighbor1 = neighbor;
        else if (neighbor2 == -1) {
            neighbor2 = neighbor;
            break;
        }
        else {
            assert(false);
        }
    }
    assert(neighbor2 != -1);

    int const vertex1(neighbor1);
    int const vertex2(neighbor2);
    int const vertex1degree = degree(vertex1);
    int const vertex2degree = degree(vertex2);
    int smallDegreeNeighbor = vertex1degree > vertex2degree ? vertex2 : vertex1;
    int highDegreeNeighbor = vertex1degree > vertex2degree ? vertex1 : vertex2;

    assert(smallDegreeNeighbor != highDegreeNeighbor);
    assert(smallDegreeNeighbor < m_AdjacencyArray.size());
    assert(highDegreeNeighbor < m_AdjacencyArray.size());
    assert(degree(highDegreeNeighbor) + degree(smallDegreeNeighbor) == vertex1degree + vertex2degree);

    if(isBoundaryVertex(smallDegreeNeighbor)) {
        return false;
    }

    for(int const neighbor2 : neighbors[smallDegreeNeighbor]) if(inGraph.Contains(neighbor2)) {
            if(neighbor2 == highDegreeNeighbor) {
                return false; // neighbors must not be adjacent.
            }
    }

    foldedVertexCount += 2;

    assert(partitions[highDegreeNeighbor] == partition);
    assert(partitions[smallDegreeNeighbor] == partition);

    Reduction reduction(FOLDED_VERTEX);
    reduction.SetVertex(vertex);
    reduction.AddNeighbor(highDegreeNeighbor);
    reduction.AddNeighbor(smallDegreeNeighbor);
    reduction.SetKeptVertex(highDegreeNeighbor);

    // if(degree(highDegreeNeighbor) > 10000) {
    //     std::cout << "Vertex degree " << degree(vertex) << std::endl;
    //     std::cout << "Neighbor1 degree: " << degree(highDegreeNeighbor) << std::endl;
    //     std::cout << "Neighbor2 degree: " << degree(smallDegreeNeighbor) << std::endl;
    // }

    // neighbors[vertex].Clear();
    neighbors[highDegreeNeighbor].Resize(neighbors[highDegreeNeighbor].Size() + degree(smallDegreeNeighbor));
    // neighbors[vertex1].Remove(vertex);
    // neighbors[vertex2].Remove(vertex);

    for (int const neighbor2 : neighbors[smallDegreeNeighbor]) if(inGraph.Contains(neighbor2)) {
        assert(partitions[neighbor2] == partition);
        if (neighbor2 == vertex) continue;
        assert(neighbor2 != highDegreeNeighbor);
        vertexDegree[neighbor2]--;
        assert(partitions[highDegreeNeighbor] == partitions[neighbor2]);
        if(!neighbors[neighbor2].Contains(highDegreeNeighbor)) {
            assert(!neighbors[highDegreeNeighbor].Contains(vertex2));
            neighbors[highDegreeNeighbor].Insert(neighbor2);
            vertexDegree[neighbor2]++;
            vertexDegree[highDegreeNeighbor]++;
            assert(neighbors[neighbor2].Contains( smallDegreeNeighbor ));
            neighbors[neighbor2].Remove(smallDegreeNeighbor);
            neighbors[neighbor2].Insert(highDegreeNeighbor);
        }
        INSERT_REMAINING(partition, remaining, neighbor2);
        // remaining.Insert(neighbor2);
    }


    INSERT_REMAINING(partition, remaining, highDegreeNeighbor);
    // remaining.Insert(vertex);

    vReductions.emplace_back(std::move(reduction));

    remaining.Remove(smallDegreeNeighbor);
    remaining.Remove(vertex);
    vertexDegree[highDegreeNeighbor]--;
    REMOVE_VERTEX(partition, smallDegreeNeighbor);
    REMOVE_VERTEX(partition, vertex);

    profilingAddTimeSuccessfulFold(&profilingHelper, partition);
    // std::cout << "Degree after fold: " << degree(highDegreeNeighbor) << std::endl;
    return true;
}


void parallel_reductions::UpdateRemaining(vector<ArraySet> &remainingPerPartition, vector<vector<int>> &bufferPerPartition) {
    int numPartitions = remainingPerPartition.size();
#pragma omp parallel for schedule(dynamic)
    for(int partition = 0; partition < numPartitions; ++partition) {
      auto tid = omp_get_thread_num();
      vector<int> &buffer = bufferPerPartition[tid];
        ArraySet &remaining = remainingPerPartition[partition];
        int numVerticesRemoved = 0;
        for(int vertex: inGraphPerPartition[partition]) {
            assert(partitions[vertex] == partition);
            if(!inGraph.Contains(vertex)) {
                neighborhoodChanged.Remove(vertex);
                remaining.Remove(vertex);
                buffer[numVerticesRemoved++] = vertex;
            }
            else if(neighborhoodChanged.Contains(vertex)) {
                remaining.Insert(vertex);
            }
        }
        for(int i = 0; i < numVerticesRemoved; ++i) {
            int vertex = buffer[i];
            inGraphPerPartition[partition].Remove(vertex);
        }
    }
}

bool parallel_reductions::LPReduction(vector<ArraySet> &remainingPerPartition, vector<vector<int>> &bufferPerPartition, int &numLPReductions) {
    int sizeBefore = inGraph.Size();
    int N = neighbors.size();
    UpdateRemaining(remainingPerPartition, bufferPerPartition);
        bool changed = false;
#pragma omp parallel for
    for(int vertex = 0; vertex < N; ++vertex) {
        if(!inGraph.Contains(vertex))
            continue;
        if(maximumMatching.reachableVertices[vertex] == 0 && maximumMatching.reachableVertices[vertex + N] > 0) {
            changed = true;
            // vertex is in the vertex cover
            independent_set[vertex] = 1;
            inGraph.Remove(vertex);
            for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                vertexDegree[neighbor]--;
                if(partitions[vertex] != partitions[neighbor])
                    numCutEdges[neighbor]--;
                neighborhoodChanged.Insert(neighbor);
            }
            neighbors[vertex].Clear();
        } else if (maximumMatching.reachableVertices[vertex] > 0 && maximumMatching.reachableVertices[vertex + N] == 0) {
            changed = true;
            // vertex is in the independent set
            // Nothing to to for the neighbors because they get removed too (two vertices on the same edge can't be 0)
            independent_set[vertex] = 0;
            inGraph.Remove(vertex);
            neighbors[vertex].Clear();
        }
        // else: We don't know it
    }
    UpdateRemaining(remainingPerPartition, bufferPerPartition);

    int sizeAfter = inGraph.Size();
    numLPReductions += sizeBefore - sizeAfter;
    return changed;
}


void parallel_reductions::initGlobalBurstEstimator() {
    global_burst_timer = omp_get_wtime();
    global_burst_estimation = -1.0;
    lastFinishedThread = -1;
    terminationFlag = false;;
    firstFinished = -1;
}
double parallel_reductions::finishThreadAndGetEstimatedBurstLength(int tid) {
    global_burst_timer_mutex.lock();
    double current_time = omp_get_wtime();
    bool firstFinisher = lastFinishedThread == -1;
    lastFinishedThread = tid;
    double lastBurstLength = current_time - global_burst_timer;
    if (!firstFinisher) {
        global_burst_estimation = global_burst_estimation <= 0.0 ? lastBurstLength : GLOBAL_BURST_ESTIMATION_ALPHA * lastBurstLength + (1 - GLOBAL_BURST_ESTIMATION_ALPHA) * global_burst_estimation;
    } else {
        firstFinished = tid;
    }
    global_burst_timer = current_time;
    global_burst_timer_mutex.unlock();
    if(firstFinisher) {
        return 0.0;
    } else {
        return global_burst_estimation;
    }
}
bool parallel_reductions::isLastFinishedThread(int tid) {
    return lastFinishedThread == tid && tid != firstFinished;
}
bool parallel_reductions::shouldTerminate() {
    return terminationFlag;
}
void parallel_reductions::terminateOtherThreads() {
    terminationFlag = true;
}

size_t parallel_reductions::get_current_is_size_with_folds(){
    return is_offset;
}

void parallel_reductions::getGraphAccess(graph_access &G, std::vector<unsigned int> &reverse_mapping) {
    // Number of nodes
    unsigned int const node_count = inGraph.Size();
    // Number of edges
    int m = 0;

    // Nodes -> Range
    // std::vector<NodeID> mapping(neighbors.size(), UINT_MAX);
    graph_to_kernel_map = std::vector<int> (m_AdjacencyArray.size());

    // Get number of edges and reorder nodes
    unsigned int node_counter = 0;
    for (int node = 0; node < (int)neighbors.size(); ++node) if (inGraph.Contains(node)) {
        for (int const neighbor : neighbors[node]) if (inGraph.Contains(neighbor)) m++;
        graph_to_kernel_map[node] = node_counter;
        reverse_mapping[node_counter] = node;
        node_counter++;
    }

    // Create the adjacency array
    std::vector<int> xadj(node_count + 1);
    std::vector<int> adjncy(m);
    unsigned int adjncy_counter = 0;
    for (unsigned int i = 0; i < node_count; ++i) {
        xadj[i] = adjncy_counter;
        for (int const neighbor : neighbors[reverse_mapping[i]]) if(inGraph.Contains(neighbor)){
            // if (graph_to_kernel_map[neighbor] == i) continue;
            // if (graph_to_kernel_map[neighbor] == UINT_MAX) continue;
            adjncy[adjncy_counter++] = graph_to_kernel_map[neighbor];
        }
        std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
    }
    xadj[node_count] = adjncy_counter;

    // Build the graph
    G.build_from_metis(node_count, &xadj[0], &adjncy[0]);
    // std::cout << "m_adj: " << m_AdjacencyArray.size() << "; neighbors: " << neighbors.size() << std::endl;
}

void parallel_reductions::ExtendPartialSolution(std::vector<bool> &in_out_independent_set){
    // std::cout << "calling Extend" << std::endl;
    // std::cout << "m_adj: " << m_AdjacencyArray.size() << "; neighbors: " << neighbors.size() << std::endl;
    // assert(in_out_independent_set.size() == neighbors.size());
    // if(in_out_independent_set.size() != neighbors.size()){
    //     std::cout << "in out is has wrong size (neighbors)" << std::endl;
    // }
    // assert(in_out_independent_set.size() == independent_set.size());
    // if(in_out_independent_set.size() != independent_set.size()){
    //     std::cout << "in out is has wrong size (independent_set)" << std::endl;
    // }
    // assert(in_out_independent_set.size() == m_AdjacencyArray.size());
    // if(in_out_independent_set.size() != m_AdjacencyArray.size()) {
    //     std::cout << "in out is has wrong size (m_AdjacencyArray)" << std::endl;
    //     // std::cout << "in out is has wrong size (m_AdjacencyArray)" << in_out_independent_set.size() << " : " << m_AdjacencyArray.size() << std::endl;
    // }


    int size_before = is_offset;
    for(bool is : in_out_independent_set) {
        if(is)
            size_before++;
    }

    int size_without_folds = 0;
    for(bool is : independent_set) {
        if(is == 0)
            size_without_folds++;
    }

    // std::cout << "gathering is" << std::endl;
    // #pragma omp parallel for
    for(int node = 0; node < (int)neighbors.size(); ++node) {
        // std::cout << "gathering is for " << node << std::endl;
        if(inGraph.Contains(node)) {
            if( independent_set[node] != -1 ) {
                std::cout << "Vertex already removed or added" << std::endl;
                exit(1);
            }
            independent_set[node] = in_out_independent_set[node] == true ? 0 : 1;
            // independent_set[node] = 1;
        } else if(in_out_independent_set[node] == true) {
            std::cout << "Error. Non contained vertex is in is" << std::endl;
            exit(1);
        }
    }
    // std::cout << "undoing reductions" << std::endl;
    // #pragma omp parallel for
    for(int i = AllReductions.size(); i > 0; i--) {
        for(auto Reductions: AllReductions[i - 1]) {
            ApplyKernelSolutionToReductions(Reductions);
        }
    }

    // std::cout << "updating output" << std::endl;
    // update full independent set
    for (unsigned int i = 0; i < neighbors.size(); ++i) {
        if (independent_set[i] == 0) {
            in_out_independent_set[i] = true;
////        in_out_independent_set_size++;
        }else if(independent_set[i] == 1) {
            in_out_independent_set[i] = false;
        } else if(independent_set[i] == -1) {
            std::cout << "Error! Vertex still -1!" << std::endl;
        }
    }

    int size_after = 0;
    for(bool is : in_out_independent_set) {
        if(is)
            size_after++;
    }

    if(size_after == size_before) {
        //std::cout << "Same size" << std::endl;
    } else {
        std::cout << "ERROR: Before: " << size_before << "; after: " << size_after << "; without folds: " << size_without_folds <<  std::endl;
    }
////    if (in_out_independent_set_size != new_in_out_independent_set_vertices + get_current_is_size()) {
////        cout << "ERROR: incorrect original count for independent set size with reductions!" << endl << flush;
////    }
    in_out_independent_set.resize(neighbors.size());
}

void parallel_reductions::force_into_independent_set(std::vector<unsigned int> &nodes_to_force) {
    // std::cout << "calling force" << std::endl;
    for(int vertex: nodes_to_force) {
        // std::cout << "testing " << vertex << std::endl;
        // std::cout << "inGraph.Contains(vertex): " <<inGraph.Contains(vertex) << std::endl;
        // std::cout << "neighbors.size(): " << neighbors.size() << std::endl;
        assert(vertex < neighbors.size());
        if(vertex >= (int)neighbors.size()) {
            std::cout << "Trying to force invalid vertex" << std::endl;
        }
        assert(inGraph.Contains(vertex));
        if(!inGraph.Contains(vertex)) {
            std::cout << "Trying to force vertex that has already been removed" << std::endl;
        }
        // std::cout << "Forcing " << vertex << std::endl;
        is_offset++;
        inGraph.Remove(vertex);
        independent_set[vertex] = 0;
        for(int neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
                independent_set[neighbor] = 1;
                inGraph.Remove(neighbor);
                int neighborPartition = partitions[neighbor];
                for(int nextNeighbor : neighbors[neighbor]) {
                    vertexDegree[nextNeighbor]--;
                    if(partitions[nextNeighbor] != neighborPartition) {
                        numCutEdges[nextNeighbor]--;
                    }
                }
            REMOVE_NEIGHBOR(partitions[vertex], neighbor, vertex);
        }
    }
}

void parallel_reductions::reduce_graph_parallel(std::vector<unsigned int> forced_vertices) {

    int numPartitions = partition_nodes.size();
    auto numPartitionsRealPointer = max_element(partitions.begin(), partitions.end());
    int numPartitionsReal = numPartitionsRealPointer!=partitions.end() ? *numPartitionsRealPointer + 1 : 0;

    omp_set_num_threads(numPartitions);
    long numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    //std::cout << "num threads: " << numThreads << std::endl;

    profilingInit(&profilingHelper, &neighbors, numPartitions);

    inGraphPerPartition = vector<ArraySet>(numPartitions);

    vector<vector<bool>> vMarkedVerticesPerTid(numThreads);
    vector<vector<int>> tempInt1PerTid(numThreads);
    vector<vector<int>> tempInt2PerTid(numThreads);
    vector<fast_set> fastSetPerTid(numThreads, fast_set(0));
    vector<vector<int>> tempIntDoubleSizePerTid(numThreads);
    vector<ArraySet> remainingPerPartition(numPartitions);
    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(numPartitions);
    for(int partition = 0; partition < numPartitions; partition++) {
      ArraySet remaining = ArraySet(m_AdjacencyArray.size());
      remainingPerPartition[partition] = remaining;
      inGraphPerPartition[partition] = ArraySet(m_AdjacencyArray.size());
      for (int const vertex : partition_nodes[partition]) {
        if(inGraph.Contains(vertex)) {
          assert(partitions[vertex] == partition);
          inGraphPerPartition[partition].Insert(vertex);
        }
      };
    }
#pragma omp parallel for
    for(int tid = 0; tid < numThreads; tid++) {
        // std::cout << "Start allocating memory for block " << tid << std::endl;
        vMarkedVerticesPerTid[tid] = std::vector<bool>(m_AdjacencyArray.size(), false);
        tempInt1PerTid[tid] = vector<int>(m_AdjacencyArray.size());
        tempInt2PerTid[tid] = vector<int>(m_AdjacencyArray.size());
        fastSetPerTid[tid] = fast_set(m_AdjacencyArray.size());
        tempIntDoubleSizePerTid[tid] = vector<int>(m_AdjacencyArray.size() * 2);
    }
    // std::cout << "Finished allocating memory" << std::endl;

    vector<double> partitionTimes(numPartitions);
    vector<double> partitionFinishTimes(numPartitions);
    vector<int> partitionFinishSizes(numPartitions);
    vector<int> numIsolatedCliqueReductions(numPartitions, 0);
    vector<int> numVertexFoldReductions(numPartitions, 0);
    vector<int> numTwinReductionsRemoved(numPartitions, 0);
    vector<int> numTwinReductionsFolded(numPartitions, 0);
    vector<int> removedUnconfinedVerticesCount(numPartitions, 0);
    vector<int> numDiamondReductions(numPartitions, 0);
    dependecy_checking_burst_estimation = std::vector<double>(numPartitions, -1.0);
    dependency_checking_times = std::vector<double>(numPartitions, 0.0);

    global_burst_timer = 0.0;
    global_burst_estimation = -1.0;
    lastFinishedThread = -1;
    terminationFlag = false;

    int numLPReductions = 0;
    
    assert(checkDegrees());

    double startClock = omp_get_wtime();
    double tmpClock;
    double LPTime = 0;
    double restTime = 0;
    // std::cout << "Filling remaining vertices" << std::endl;
    // if(forced_vertices.size() == 0) {
#pragma omp parallel for schedule(dynamic)
        for(int partition = 0; partition < numPartitions; ++partition) {
            remainingPerPartition[partition].Clear();
            for (int const vertex : partition_nodes[partition]) {
                if(inGraph.Contains(vertex)) {
                    assert(partitions[vertex] == partition);
                    remainingPerPartition[partition].Insert(vertex);
                }
            }
        }
//     } else {
// #pragma omp parallel for schedule(dynamic)
//         for(int partition = 0; partition < numPartitions; ++partition) {
//             remainingPerPartition[partition].Clear();
//         }
//         for(int vertex : forced_vertices) {
//             for(int neighbor : neighbors[vertex]) if(inGraph.Contains(neighbor)) {
//                 remainingPerPartition[partitions[neighbor]].Insert(neighbor);
//             }
//         }
//     }


    // std::cout << "Start LP reduction" << std::endl;
    tmpClock = omp_get_wtime();

    // LPReduction(remainingPerPartition, tempInt1PerTid, numLPReductions);
    // LPTime += omp_get_wtime() - tmpClock;

    // std::cout << "done with LP reduction" << std::endl;

    bool changed = true;
    int numIterations = 0;
    double terminationTime = 0.0;
    while(changed) {
        terminationTime = -1.0;
        initGlobalBurstEstimator();
        // std::cout << "Iteration " << numIterations << " starts at " << omp_get_wtime() - startClock << " current size: " << inGraph.Size() << std::endl;
        // int sizeBefore = inGraph.Size();
        // std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
      // std::cout << "starting new iteration: " << numIterations << std::endl;
      // std::cout << "Current rest time: " << restTime << std::endl;
      auto partitionTimesCopy(partitionTimes);
      //auto old_graphsize = inGraph.Size();
        int graphSizeLastSample = 0;
        // for(int partition = 0; partition < numPartitions; ++partition) {
        //     graphSizeLastSample += inGraphPerPartition[partition].Size();
        // }
        graphSizeLastSample = inGraph.Size();
        double timeLastSample = omp_get_wtime();
        double timeBefore = timeLastSample;
        int graphSizeBefore = graphSizeLastSample;

        atomic_int finishedThreads(0);
        tmpClock = omp_get_wtime();
#pragma omp parallel for schedule(dynamic,1)
        for(int partition = 0; partition < numPartitions; partition++) {
          auto tid = omp_get_thread_num();
          // std::cout << "partition " << partition << " on tid " << tid << std::endl;
          //std::cout << partition << ": starting new iteration" << std::endl;
            ApplyReductions(partition, ReductionsPerPartition[partition], vMarkedVerticesPerTid[tid], remainingPerPartition[partition], tempInt1PerTid[tid], tempInt2PerTid[tid], fastSetPerTid[tid], tempIntDoubleSizePerTid[tid], partitionTimes[partition], numIsolatedCliqueReductions[partition], numVertexFoldReductions[partition], numTwinReductionsRemoved[partition], numTwinReductionsFolded[partition], removedUnconfinedVerticesCount[partition], numDiamondReductions[partition]);
            partitionFinishTimes[partition] = omp_get_wtime() - startClock;
            partitionFinishSizes[partition] = inGraph.Size();
            // if(!shouldTerminate()) {
            //     double sleeptime = finishThreadAndGetEstimatedBurstLength(tid) * GLOBAL_BURST_THRESHOLD_MULTIPLIER;
            //     double startTime = omp_get_wtime();
            //     std::cout << tid << " sleeping for " << sleeptime << std::endl;
            //     while( (omp_get_wtime() - startTime) < sleeptime);
            //     if(isLastFinishedThread(tid)) {
            //         std::cout << tid << " Terminating others" << std::endl;
            //         terminationTime = omp_get_wtime() - startClock;
            //         terminateOtherThreads();
            //     }
            // }
            finishedThreads++;
            global_burst_timer_mutex.lock();
            bool isFirstFinisher = firstFinished == -1;
            firstFinished = 1;
            global_burst_timer_mutex.unlock();
            if(isFirstFinisher) {
                int graphSizeCurrentSample = 0;
                // for(int partition = 0; partition < numPartitions; ++partition) {
                //     graphSizeCurrentSample += inGraphPerPartition[partition].Size();
                // }
                double current_time = omp_get_wtime();
                // double last_delta = (graphSizeLastSample - graphSizeCurrentSample) / (current_time - timeLastSample);
                //double last_delta = 0;
                graphSizeLastSample = graphSizeBefore;
                timeLastSample = timeBefore;
                while(true) {
                    double startTime = omp_get_wtime();
                    while( (omp_get_wtime() - startTime) < 1 && finishedThreads < numPartitions);
                    if(finishedThreads == numPartitions)
                        break;
                    graphSizeCurrentSample = 0;
                    graphSizeCurrentSample = inGraph.Size();
                    current_time = omp_get_wtime();
                    double current_delta = (graphSizeLastSample - graphSizeCurrentSample) / (current_time - timeLastSample);
                    double global_delta = (graphSizeBefore - graphSizeCurrentSample) / (current_time - timeBefore);
                    graphSizeLastSample = graphSizeCurrentSample;
                    timeLastSample = current_time;
                    if(current_delta <= 0.05 * global_delta) {
                        terminationTime = omp_get_wtime() - startClock;
                        terminateOtherThreads();
                        break;
                    }
                }
            }
        }
        restTime += omp_get_wtime() - tmpClock;
        tmpClock = omp_get_wtime();

        if(numPartitionsReal != numPartitions) {
            omp_set_num_threads(numPartitionsReal);
        }
        changed = LPReduction(remainingPerPartition, tempInt1PerTid, numLPReductions);
        if(numPartitionsReal != numPartitions) {
            omp_set_num_threads(numPartitions);
        }
        LPTime += omp_get_wtime() - tmpClock;
        numIterations++;
    }

    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    int sum_vertex_fold = std::accumulate(numVertexFoldReductions.begin(), numVertexFoldReductions.end(), 0);
    int sum_twin_folded = std::accumulate(numTwinReductionsFolded.begin(), numTwinReductionsFolded.end(), 0);
    assert(checkDegrees());

    is_offset = 0;
    is_offset += sum_vertex_fold / 2;
    is_offset += sum_twin_folded / 2;
    for(int vertex = 0; vertex < (int)m_AdjacencyArray.size(); vertex++) {
      if(independent_set[vertex] == 0) {
        is_offset++;
      }
    }
    //std::cout << "Independent set offset: " << is_offset << std::endl;
}

bool parallel_reductions::checkDegrees() {
    for(int vertex = 0; vertex < neighbors.size(); ++vertex) if(inGraph.Contains(vertex)) {
        int deg = 0;
        int cutEdges = 0;
        for(int const neighbor: neighbors[vertex]) if(inGraph.Contains(neighbor)) {
            deg++;
            if(partitions[vertex] != partitions[neighbor])
                cutEdges++;
        }
        if(vertexDegree[vertex] != deg) {
            std::cout << "Vertex " << vertex << " has degree " << deg << " but vertexDegree[vertex] has value " << vertexDegree[vertex] << " number of elements in neighbors[vertex]: " << neighbors[vertex].Size() << std::endl;
            return false;
        }
        if(numCutEdges[vertex] != cutEdges){
            std::cout << "Vertex " << vertex << " has " << cutEdges << " cut edges but numCutEdges[vertex] has value " << numCutEdges[vertex] << std::endl;
            return false;
        }
    }
    return true;
}

void parallel_reductions::reduce_graph_sequential() {
    // return;
    long numThreads;

    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    omp_set_num_threads(1);
    // std::cout << "numThreads: " << numThreads << std::endl;
    profilingInit(&profilingHelper, &neighbors, 1);

    int N = m_AdjacencyArray.size();

    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(1);
    std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
    vector<int> tempInt1 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempInt2 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempIntDoubleSize = vector<int>(m_AdjacencyArray.size() * 2);
    fast_set fastSet(m_AdjacencyArray.size());
    ArraySet remaining = ArraySet(N);

    partition_nodes.resize(1);
    partition_nodes[0].clear();
    inGraphPerPartition = vector<ArraySet>(1);
    inGraphPerPartition[0] = ArraySet(N);
    for(int vertex = 0; vertex < N; ++vertex) {
        numCutEdges[vertex] = 0;
        partitions[vertex] = 0;
        partition_nodes[0].push_back(vertex);
        if(inGraph.Contains(vertex)) {
            inGraphPerPartition[0].Insert(vertex);
        }
    }

    double time(0);
    int numIsolatedCliqueReductions(0);
    int numVertexFoldReductions(0);
    int numTwinReductionsRemoved(0);
    int numTwinReductionsFolded(0);
    int removedUnconfinedVerticesCount(0);
    int numLPReductions(0);
    int numDiamondReductions(0);
    
    double startClock = omp_get_wtime();

    vector<ArraySet> remainingPerPartition = {remaining};
    vector<vector<int>> tempInt1PerPartition = {tempInt1};

    std::cout << "Filling remaining vertices set..." << std::endl;
    remaining.Clear();
    for(int vertex = 0; vertex < N; ++vertex)  {
        if(inGraph.Contains(vertex)) {
            assert(partitions[vertex] == 0);
            remaining.Insert(vertex);
        }
    }

    bool changed = true;
    int numIterations = 0;
    initGlobalBurstEstimator();
    while(changed) {
        numIterations++;
        ApplyReductions(0, ReductionsPerPartition[0], vMarkedVertices, remaining, tempInt1, tempInt2, fastSet, tempIntDoubleSize, time, numIsolatedCliqueReductions, numVertexFoldReductions, numTwinReductionsRemoved, numTwinReductionsFolded, removedUnconfinedVerticesCount, numDiamondReductions);
        changed = LPReduction(remainingPerPartition, tempInt1PerPartition, numLPReductions);
    }

    double endClock = omp_get_wtime();
    std::cout << "Num iterations: " << numIterations << std::endl;
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
    cout << "Number of isolated clique reductions: " << numIsolatedCliqueReductions << endl;
    cout << "Number of vertex fold reductions: " << numVertexFoldReductions << endl;
    cout << "Number of twin reductions (removed): " << numTwinReductionsRemoved << endl;
    cout << "Number of twin reductions (folded): " << numTwinReductionsFolded << endl;
    cout << "Number of unconfined vertices removed: " << removedUnconfinedVerticesCount << endl;
    cout << "Number of diamond reductions: " << numDiamondReductions << endl;
    cout << "Number of vertices removed by LP reduction: " << numLPReductions << endl;
    omp_set_num_threads(numThreads);
}

void parallel_reductions::initDependencyCheckingEstimation(int partition) {
    dependency_checking_times[partition] = omp_get_wtime();
}

void parallel_reductions::updateDependencyCheckingEstimation(int partition) {
    double current_time = omp_get_wtime();
    double timeSinceLastReduction = current_time - dependency_checking_times[partition];
    dependecy_checking_burst_estimation[partition] = dependecy_checking_burst_estimation[partition] <= 0.0 ? timeSinceLastReduction : DEPENDENCYCHECKING_BURST_ESTIMATION_ALPHA * timeSinceLastReduction + (1 - DEPENDENCYCHECKING_BURST_ESTIMATION_ALPHA) * dependecy_checking_burst_estimation[partition];
    dependency_checking_times[partition] = current_time;
}

bool parallel_reductions::shouldStopDependencyCheckingReductions(int partition) {
    return dependecy_checking_burst_estimation[partition] <= 0.0 ? false : ((omp_get_wtime() - dependency_checking_times[partition]) > dependecy_checking_burst_estimation[partition] * DEPENDENCY_CHECKING_THRESHOLD_MULTIPLIER);
}

void parallel_reductions::ApplyReductions(int const partition, vector<Reduction> &vReductions, std::vector<bool> &vMarkedVertices, ArraySet &remaining, vector<int> &tempInt1, vector<int> &tempInt2, fast_set &fastSet, vector<int> &tempIntDoubleSize, double &time, int &isolatedCliqueCount, int &foldedVertexCount, int &removedTwinCount, int &foldedTwinCount, int &removedUnconfinedVerticesCount, int &numDiamondReductions)
{
    double startClock = omp_get_wtime();
    bool changed = true;
    // Only do this for singe threaded and only in first iteration
    bool finishedOneIteration = inGraphPerPartition[partition].Size() != neighbors.size();
    while (changed) {
        changed = false;
        initDependencyCheckingEstimation(partition);
        while (!remaining.Empty()) {
            if(shouldTerminate() && finishedOneIteration) {
                break;
            }
            int const vertex = *(remaining.begin());
            remaining.Remove(vertex);
            assert(inGraph.Contains(vertex));
            assert(partitions[vertex] == partition);
            assert(independent_set[vertex] == -1);
            bool reduction = RemoveIsolatedClique(partition, vertex, vReductions, remaining, vMarkedVertices, isolatedCliqueCount);
            if (!reduction && m_bAllowVertexFolds) {
                reduction = FoldVertex(partition, vertex, vReductions, remaining, foldedVertexCount);
		}
            if (!reduction) {
                reduction = removeTwin(partition, vertex, vReductions, remaining, vMarkedVertices, removedTwinCount, foldedTwinCount);
		}
        }
        std::vector<int> verticesToRemove;
        for (int const vertex : inGraphPerPartition[partition]) {
            if(shouldTerminate() && finishedOneIteration) {
                break;
            }
            assert(partitions[vertex] == partition);
            if(inGraph.Contains(vertex)) {
                bool reduction = removeUnconfined(partition, vertex, remaining, fastSet, tempInt1, tempInt2, tempIntDoubleSize, removedUnconfinedVerticesCount, numDiamondReductions);
                if(reduction) {
                    changed = true;
                    verticesToRemove.push_back(vertex);
                }
            } else {
                verticesToRemove.push_back(vertex);
            }
        }
        for(int vertex : verticesToRemove) {
            inGraphPerPartition[partition].Remove(vertex);
        }
        finishedOneIteration = true;
        if(shouldTerminate()) {
            changed = false;
        }
	}
    double endClock = omp_get_wtime();
    time += (endClock - startClock);
}

void parallel_reductions::reduce_graph_sequential_reduction_wise() {
    long numThreads;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
    omp_set_num_threads(1);
    profilingInit(&profilingHelper, &neighbors, 1);

    int N = m_AdjacencyArray.size();

    vector<vector<Reduction>> ReductionsPerPartition = vector<vector<Reduction>>(1);
    std::vector<bool> vMarkedVertices = std::vector<bool>(m_AdjacencyArray.size(), false);
    vector<int> tempInt1 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempInt2 = vector<int>(m_AdjacencyArray.size());
    vector<int> tempIntDoubleSize = vector<int>(m_AdjacencyArray.size() * 2);
    fast_set fastSet(m_AdjacencyArray.size());
    ArraySet remainingUse = ArraySet(N);
    ArraySet remaining2 = ArraySet(N);

    partition_nodes.resize(1);
    partition_nodes[0].clear();
    inGraphPerPartition = vector<ArraySet>(1);
    inGraphPerPartition[0] = ArraySet(N);
    for(int vertex = 0; vertex < N; ++vertex) {
        numCutEdges[vertex] = 0;
        partitions[vertex] = 0;
        partition_nodes[0].push_back(vertex);
        if(inGraph.Contains(vertex)) {
            inGraphPerPartition[0].Insert(vertex);
        }
    }

    double time(0);
    int numIsolatedCliqueReductions(0);
    int numVertexFoldReductions(0);
    int numTwinReductionsRemoved(0);
    int numTwinReductionsFolded(0);
    int removedUnconfinedVerticesCount(0);
    int numLPReductions(0);
    int numDiamondReductions(0);
    
    double startClock = omp_get_wtime();

    vector<vector<int>> tempInt1PerPartition = {tempInt1};

    remainingUse.Clear();
    for(int vertex = 0; vertex < N; ++vertex)  {
        if(inGraph.Contains(vertex)) {
            assert(partitions[vertex] == 0);
            remainingUse.Insert(vertex);
        }
    }

    ArraySet *remainingUseptr = &remainingUse;
    ArraySet *remainingInsertptr = &remaining2;
    double isolated_clique_time = 0.0;
    double unconfined_time = 0.0;
    double lp_time = 0.0;
    double vertex_fold_time = 0.0;
    double twin_time = 0.0;
    while(true) {
        bool changed = false;
        ArraySet *temp;
        double start_time;

        start_time = omp_get_wtime();
        changed = RemoveAllIsolatedClique(0, ReductionsPerPartition[0], remainingUseptr, remainingInsertptr, vMarkedVertices, numIsolatedCliqueReductions);
        isolated_clique_time += omp_get_wtime() - start_time;
        temp = remainingInsertptr;
        remainingInsertptr = remainingUseptr;
        remainingUseptr = temp;
        // if(changed) continue;

        start_time = omp_get_wtime();
        changed = removeAllUnconfined(0, remainingUseptr, fastSet, tempInt1, tempInt2, tempIntDoubleSize, removedUnconfinedVerticesCount, numDiamondReductions);
        unconfined_time += omp_get_wtime() - start_time;
        if(changed) continue;

        vector<ArraySet> remainingPerPartition = {*remainingUseptr};
        start_time = omp_get_wtime();
        changed = LPReduction(remainingPerPartition, tempInt1PerPartition, numLPReductions);
        lp_time += omp_get_wtime() - start_time;
        if(changed) continue;

        start_time = omp_get_wtime();
        changed = FoldAllVertices(0, ReductionsPerPartition[0], remainingUseptr, remainingInsertptr, numVertexFoldReductions);
        vertex_fold_time += omp_get_wtime() - start_time;
        temp = remainingInsertptr;
        remainingInsertptr = remainingUseptr;
        remainingUseptr = temp;
        if(changed) continue;

        start_time = omp_get_wtime();
        changed = removeAllTwin(0, ReductionsPerPartition[0], remainingUseptr, remainingInsertptr, vMarkedVertices, numTwinReductionsRemoved, numTwinReductionsFolded);
        twin_time += omp_get_wtime() - start_time;
        temp = remainingInsertptr;
        remainingInsertptr = remainingUseptr;
        remainingUseptr = temp;
        if(changed) continue;

        break;
    }

    double endClock = omp_get_wtime();
    AllReductions.push_back(ReductionsPerPartition);
    profilingPrint(&profilingHelper);

    cout << "Total time spent applying reductions  : " << (endClock - startClock) << endl;
    cout << "Number of isolated clique reductions: " << numIsolatedCliqueReductions << endl;
    cout << "Number of vertex fold reductions: " << numVertexFoldReductions << endl;
    cout << "Number of twin reductions (removed): " << numTwinReductionsRemoved << endl;
    cout << "Number of twin reductions (folded): " << numTwinReductionsFolded << endl;
    cout << "Number of unconfined vertices removed: " << removedUnconfinedVerticesCount << endl;
    cout << "Number of diamond reductions: " << numDiamondReductions << endl;
    cout << "Number of vertices removed by LP reduction: " << numLPReductions << endl;

    cout << "Isolated clique time: " << isolated_clique_time << std::endl;
    cout << "Unconfined time: " << unconfined_time << std::endl;
    cout << "LP time: " << lp_time << std::endl;
    cout << "Vertex fold time: " << vertex_fold_time << std::endl;
    cout << "Twin time: " << twin_time << std::endl;
    omp_set_num_threads(numThreads);
}

bool parallel_reductions::removeAllUnconfined(int const partition, ArraySet *remainingInsert, fast_set &closedNeighborhood, std::vector<int> &neighborhood, std::vector<int> &numNeighborsInS, std::vector<int> &neighborsInS, int &removedUnconfinedVerticesCount, int &numDiamondReductions) {
    std::vector<int> verticesToRemove;
    bool reduced = false;
    for (int const vertex : inGraphPerPartition[partition]) {
        assert(partitions[vertex] == partition);
        if(inGraph.Contains(vertex)) {
            bool reduction = removeUnconfined(partition, vertex, *remainingInsert, closedNeighborhood, neighborhood, numNeighborsInS, neighborsInS, removedUnconfinedVerticesCount, numDiamondReductions);
            if(reduction) {
                reduced = true;
                verticesToRemove.push_back(vertex);
            }
        }
    }
    for(int vertex : verticesToRemove) {
        inGraphPerPartition[partition].Remove(vertex);
    }
    return reduced;
}
bool parallel_reductions::removeAllTwin(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, std::vector<bool> &vMarkedVertices, int &removedTwinCount, int &foldedTwinCount) {
    bool reduced = false;
    while (!(*remainingUse).Empty()) {
        int const vertex = *((*remainingUse).begin());
        (*remainingUse).Remove(vertex);
        if(!inGraph.Contains(vertex))
            continue;
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        bool reduction = removeTwin(partition, vertex, vReductions, *remainingInsert, vMarkedVertices, removedTwinCount, foldedTwinCount);
        if(!reduction) {
            remainingInsert->Insert(vertex);
        }
        reduced |=reduction;
    }
    return reduced;
}

bool parallel_reductions::RemoveAllIsolatedClique(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, std::vector<bool> &vMarkedVertices, int &isolatedCliqueCount){
    bool reduced = false;
    while (!(*remainingUse).Empty()) {
        int const vertex = *((*remainingUse).begin());
        (*remainingUse).Remove(vertex);
        if(!inGraph.Contains(vertex))
            continue;
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        bool reduction = RemoveIsolatedClique(partition, vertex, vReductions, *remainingInsert, vMarkedVertices, isolatedCliqueCount);
        if(!reduction) {
            remainingInsert->Insert(vertex);
        }
        reduced |=reduction;
    }
    return reduced;
}

bool parallel_reductions::FoldAllVertices(int const partition, std::vector<Reduction> &vReductions, ArraySet *remainingUse, ArraySet *remainingInsert, int &foldedVertexCount) {
    bool reduced = false;
    while (!(*remainingUse).Empty()) {
        int const vertex = *((*remainingUse).begin());
        (*remainingUse).Remove(vertex);
        if(!inGraph.Contains(vertex))
            continue;
        assert(partitions[vertex] == partition);
        assert(independent_set[vertex] == -1);
        bool reduction = FoldVertex(partition, vertex, vReductions, *remainingInsert, foldedVertexCount);
        if(!reduction) {
            remainingInsert->Insert(vertex);
        }
        reduced |=reduction;
    }
    return reduced;
}

void parallel_reductions::UndoReductions(vector<Reduction> const &vReductions)
{
#ifdef TIMERS
    clock_t startClock = clock();
#endif //TIMERS
    for (size_t index = vReductions.size(); index > 0; index--) {
        Reduction const &reduction(vReductions[index-1]);
        switch(reduction.GetType()) {
            case ISOLATED_VERTEX:

                inGraph.Insert(reduction.GetVertex());
                independent_set[reduction.GetVertex()] = 0;
                for (int const neighbor : reduction.GetNeighbors()) {
                    inGraph.Insert(neighbor);
                    independent_set[neighbor] = 1;
                }
                for (pair<int,int> const &edge : reduction.GetRemovedEdges()) {
                    neighbors[edge.first].Insert(edge.second);
                }
            break;
            case FOLDED_VERTEX:
                inGraph.Insert(reduction.GetNeighbors()[0]);
                inGraph.Insert(reduction.GetNeighbors()[1]);
                // first remove all added edges
                for (int const neighbor : neighbors[reduction.GetVertex()]) {
                    neighbors[neighbor].Remove(reduction.GetVertex());
                }
                neighbors[reduction.GetVertex()].Clear();
                // then replace all removed edges
                for (pair<int,int> const &edge : reduction.GetRemovedEdges()) {
                    neighbors[edge.first].Insert(edge.second);
                }
                if(independent_set[reduction.GetVertex()] == 0) {
                    independent_set[reduction.GetNeighbors()[0]] = 0;
                    independent_set[reduction.GetNeighbors()[1]] = 0;
                    independent_set[reduction.GetVertex()] = 1;
                } else {
                    independent_set[reduction.GetNeighbors()[0]] = 1;
                    independent_set[reduction.GetNeighbors()[1]] = 1;
                    independent_set[reduction.GetVertex()] = 0;
                }
            break;
            default:
                cout << "ERROR!: Unsupported reduction type used..." << endl << flush;
                exit(1);
            break;
        };
    }

#ifdef TIMERS
    clock_t endClock = clock();
    replaceTimer += (endClock - startClock);
    #endif // TIMERS
}

void parallel_reductions::ApplyKernelSolutionToReductions(vector<Reduction> const &vReductions)
{
#ifdef TIMERS
    clock_t startClock = clock();
#endif //TIMERS
    for (size_t index = vReductions.size(); index > 0; index--) {
        Reduction const &reduction(vReductions[index-1]);
        switch(reduction.GetType()) {
            /*case ISOLATED_VERTEX:
                independent_set[reduction.GetVertex()] = 0;
                for (int const neighbor : reduction.GetNeighbors()) {
                    independent_set[neighbor] = 1;
                }
            break;*/
            case FOLDED_VERTEX:
                assert(independent_set[reduction.GetNeighbors()[0]] == -1 || reduction.GetNeighbors()[0] == reduction.GetKeptVertex());
                assert(independent_set[reduction.GetNeighbors()[1]] == -1 || reduction.GetNeighbors()[1] == reduction.GetKeptVertex());
                assert(independent_set[reduction.GetVertex()] == -1);
                if(independent_set[reduction.GetKeptVertex()] == 0) {
                    independent_set[reduction.GetNeighbors()[0]] = 0;
                    independent_set[reduction.GetNeighbors()[1]] = 0;
                    independent_set[reduction.GetVertex()] = 1;
                } else {
                    independent_set[reduction.GetNeighbors()[0]] = 1;
                    independent_set[reduction.GetNeighbors()[1]] = 1;
                    independent_set[reduction.GetVertex()] = 0;
                }
                break;
            case FOLDED_TWINS:
                assert(reduction.GetNeighbors().size() == 3);
                assert(independent_set[reduction.GetNeighbors()[0]] == -1);
                assert(independent_set[reduction.GetNeighbors()[1]] == -1);
                assert(independent_set[reduction.GetNeighbors()[2]] == -1);
                assert(independent_set[reduction.GetTwin()] == -1);
                if(independent_set[reduction.GetVertex()] == 0) {
                    independent_set[reduction.GetVertex()] = 1;
                    independent_set[reduction.GetTwin()] = 1;
                    independent_set[reduction.GetNeighbors()[0]] = 0;
                    independent_set[reduction.GetNeighbors()[1]] = 0;
                    independent_set[reduction.GetNeighbors()[2]] = 0;
                } else {
                    independent_set[reduction.GetVertex()] = 0;
                    independent_set[reduction.GetTwin()] = 0;
                    independent_set[reduction.GetNeighbors()[0]] = 1;
                    independent_set[reduction.GetNeighbors()[1]] = 1;
                    independent_set[reduction.GetNeighbors()[2]] = 1;
                }
            break;
            default:
                cout << "ERROR!: Unsupported reduction type used..." << endl << flush;
                exit(1);
            break;
        };
    }

#ifdef TIMERS
    clock_t endClock = clock();
    replaceTimer += (endClock - startClock);
    #endif // TIMERS
}
