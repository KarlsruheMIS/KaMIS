/* 
 * This is a modified version of code written by Ariful Azad (azad@lbl.gov) and Aydin Buluc (abuluc@lbl.gov), released under the following license with the following copyright notice: 
 * 
 * ----------------------------- 
 *   *** Copyright Notice *** 
 * "Parallel Maximum Cardinality Matchings via Tree-Grafting" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved. 
 * 
 * If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov. 
 * 
 * NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so. 
 * ---------------------------------- 
 *  *** License Agreement *** 
 *  
 * ""Parallel Maximum Cardinality Matchings via Tree-Grafting" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved." 
 *  
 *  
 *  
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 
 *  
 *  
 *  
 * (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
 *  
 *  
 *  
 * (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 
 *  
 *  
 *  
 * (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
 *  
 *  
 *  
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *  
 *  
 * You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form. 
 *  
 * ------------------------------------- 
 * 
 * Adaptation to C++, addition of functionality for finding a minimum vertex cover and additional functionality for integration into ParFastKer is 
 * 
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu> 
 * */ 
 

#include <parallel/numeric>
#include "MaximumMatching.h"
#include <sys/resource.h>

void increaseStackLimit(unsigned const size) {

    const rlim_t kStackSize = size * 1024 * 1024;   // size = min stack size in MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
        {
            if (rl.rlim_cur < kStackSize)
                {
                    rl.rlim_cur = kStackSize;
                    result = setrlimit(RLIMIT_STACK, &rl);
                    if (result != 0)
                        {
                            fprintf(stderr, "setrlimit returned result = %d\n", result);
                            exit(1);
                        }
                }
        }
}
MaximumMatching::MaximumMatching(std::vector<std::vector<int>> const &adjacencyArray) {
    increaseStackLimit(16);
	G = (graph *) malloc(sizeof(graph));
	G->weight = NULL;
	long numVertices = adjacencyArray.size();
	G->vtx_pointer = new long[numVertices * 2 + 1];
	G->nrows = numVertices;
	G->n = 2 * numVertices;
	degree = (long*) malloc(sizeof(long) * (G->n));
	#pragma omp parallel for
	for(int i = 0; i < numVertices; ++i ) {
		degree[i] = adjacencyArray[i].size();
		degree[i + numVertices] = adjacencyArray[i].size();
	}
	auto end_ptr = __gnu_parallel::partial_sum(degree, degree + (G->n), (G->vtx_pointer) + 1);
	assert(end_ptr == &(G->vtx_pointer[2 * numVertices]) + 1);
	G->vtx_pointer[0] = 0;
	long numEdges = G->vtx_pointer[2 * numVertices];
	assert(G->vtx_pointer[2 * numVertices] == G->vtx_pointer[2 * numVertices - 1] + degree[2 * numVertices - 1]);
	G->m = numEdges;
	G->endV = new long[numEdges];
	#pragma omp parallel for
	for(int i = 0; i < numVertices; ++i) {
		long index = G->vtx_pointer[i];
		for(int neighbor: adjacencyArray[i]) {
			G->endV[index++] = neighbor + numVertices;
		}
		index = G->vtx_pointer[i + numVertices];
		for(int neighbor: adjacencyArray[i]) {
			G->endV[index++] = neighbor;
		}
	}


	const long NV = G->n;
	const long nrows = G->nrows; 
	QF = (long*) malloc(NV * sizeof(long));
    QFnext = (long*) malloc(NV * sizeof(long));
	flag = (long*) malloc(NV * sizeof(long));
	parent = (long*) malloc(NV * sizeof(long));
	leaf = (long*) malloc(NV * sizeof(long));
    root = (long*) malloc(NV * sizeof(long));
	mate = (long*) malloc(NV * sizeof(long));
    unmatchedU = (long*) malloc(nrows * sizeof(long));
    nextUnmatchedU = (long*) malloc(nrows * sizeof(long));

	degree1Vtx = (long*) malloc(sizeof(long) * nrows);

	reachableVertices = std::vector<int>(G->n);

	long numThreads;
#pragma omp parallel
	{
		numThreads = omp_get_num_threads();
	}

	stacks = std::vector<std::vector<long>>(numThreads);
	for(int i = 0; i < numThreads; i++)
	{
		stacks[i] = std::vector<long>(G->n);
	}

    firstKarpSipser = true;
}

int MaximumMatching::VertexDegree(const int vertex, std::vector<SparseArraySet> &neighbors, SimpleSet &inGraph, std::vector<std::atomic_int> &vertexDegree) {
    if(!inGraph.Contains(vertex))
        return 0;
    return vertexDegree[vertex];
}

void MaximumMatching::LoadGraph(std::vector<SparseArraySet> &neighbors, SimpleSet &inGraph, std::vector<std::atomic_int> &vertexDegree) {
	assert(neighbors.size() == G->nrows);
	#pragma omp parallel for
	for(int i = 0; i < G->nrows; ++i ) {
        int deg = VertexDegree(i, neighbors, inGraph, vertexDegree);
		degree[i] = deg;
		degree[i + G->nrows] = deg;
	}
	auto end_ptr = __gnu_parallel::partial_sum(degree, degree + (G->n), (G->vtx_pointer) + 1);
	assert(end_ptr == &(G->vtx_pointer[G->n]) + 1);
	G->vtx_pointer[0] = 0;
	long numEdges = G->vtx_pointer[G->n];
	G->vtx_pointer[G->n] = numEdges;
	assert(G->vtx_pointer[G->n] == G->vtx_pointer[G->n - 1] + degree[G->n - 1]);
	G->m = numEdges;

	#pragma omp parallel for
	for(int i = 0; i < G->nrows; ++i) if(degree[i] > 0) {
		long index = G->vtx_pointer[i];
        long rhsOffset = G->vtx_pointer[i + G->nrows] - index;
		for(int neighbor: neighbors[i]) if(inGraph.Contains(neighbor)) {
            // assert(index - G->vtx_pointer[i] < VertexDegree(i, neighbors, inGraph));
            G->endV[rhsOffset + index] = neighbor;
			G->endV[index++] = neighbor + G->nrows;
		}
	}
}

MaximumMatching::~MaximumMatching() {
	free(degree1Vtx);
	free(degree);

	free(flag);
	free(QF);
    free(QFnext);
	free(parent);
	free(leaf);
	free(root);
    free(unmatchedU);
    free(nextUnmatchedU);
    free(mate);

	free_graph(G);
	free(G);
}

void free_graph( graph* bGraph)
{
	delete [] bGraph->vtx_pointer;
	delete [] bGraph->endV;
    if(bGraph->weight) delete [] bGraph->weight;
}

long* MaximumMatching::MS_BFS_Graft ()
{
	const long NE = G->m; // number of edges
	const long NV = G->n; // numver of vertices in both sides
    const long nrows = G->nrows; // number of vertices in the left side
	long * __restrict__ endVertex = G->endV; // adjacency
	long * __restrict__ vtx_pointer = G->vtx_pointer; // adjacency pointer
  
    double time_start = omp_get_wtime();
    
    #define THREAD_BUF_LEN 16384
    long numUnmatchedU = 0;
    
    // identify unmatched and non-isolated vertices from where search will begin
#pragma omp parallel
    {
        long kbuf=0, nbuf[THREAD_BUF_LEN];
#pragma omp for
        for(long u=0; u<nrows; u++)
        {
            if(mate[u] == -1 && (vtx_pointer[u+1] > vtx_pointer[u]))
            {
                if (kbuf < THREAD_BUF_LEN)
                {
                    nbuf[kbuf++] = u;
                }
                else
                {
                    long voff = __sync_fetch_and_add (&numUnmatchedU, THREAD_BUF_LEN);
                    for (long vk = 0; vk < THREAD_BUF_LEN; ++vk)
                        unmatchedU[voff + vk] = nbuf[vk];
                    nbuf[0] = u;
                    kbuf = 1;
                }
                root[u] = u;
            }
            else
                root[u] = -1;
                
                
        }
        if(kbuf>0)
        {
            long voff = __sync_fetch_and_add (&numUnmatchedU, kbuf);
            for (long vk = 0; vk < kbuf; ++vk)
                unmatchedU[voff + vk] = nbuf[vk];
        }
        
    }
    
    
#pragma omp parallel for schedule(static)
    for(long i=0; i<nrows; i++)
    {
        parent[i] = -1;
        leaf[i] = -1;
        flag[i] = 0;
    }
    // I separated them out so that root information for row vertices can be set earlier
#pragma omp parallel for schedule(static)
    for(long i=nrows; i<NV; i++)
    {
        parent[i] = -1;
        leaf[i] = -1;
        root[i] = -1;
        flag[i] = 0;
    }

    // prepare frontier for the first iteration.
#pragma omp parallel for schedule(static)
	for(long i=0; i<numUnmatchedU; i++) // &&
	{
		long u  = unmatchedU[i];
        QF[i] = u;
        unmatchedU[i] = u;
	}
    
    
    long QFsize = numUnmatchedU;
    long QRsize = NV-nrows;
    long total_aug_path_len=0, total_aug_path_count=0;
    long edgeVisited = 0;
    long eFwdAll=0, eRevAll = 0, eGraftAll=0;
    double timeFwdAll=0, timeRevAll=0, timeGraftAll = 0, timeAugmentAll=0, timeStatAll=0;
    
	long iteration = 1;
	long matched = 1;
	while(matched)
	{
        double time_phase = omp_get_wtime();
        double timeFwd=0, timeRev=0;
        long  phaseEdgeVisited = 0;
        long curLayer = 0;
        long QFnextSize = 0;
        long eFwd=0, eRev = 0;
        long eFwdFrontier = 0;
        double tsLayer;
        
        // Step 1: BFS
#pragma omp parallel
        {
            long kbuf, nbuf[THREAD_BUF_LEN]; // temporary thread private Queue buffer
            while(QFsize > 0)
            {
                bool isTopdownBFS = true;
                double alpha=5;
                if(QFsize > QRsize/alpha)
                    isTopdownBFS=false;

#pragma omp single nowait
                {
                    tsLayer = omp_get_wtime();
                }
                
                
                
                kbuf=0;
                #pragma omp barrier
                //isTopdownBFS=false;
                if(isTopdownBFS) // top-down BFS
                {
                    long teFwd = 0;
#pragma omp for
                    for(long i=0; i<QFsize; i++)
                    {
                        long u = QF[i]; // fairness in Queue acess does not help... tested
                        long curRoot = root[u]; // for unmatched U root[u] = u;
                        if( leaf[curRoot] == -1) // without this test this algorithm is still correct, but continues expanding a dead tree
                        {
                            long j;
                            //for(long j=vtx_pointer[u+1]-1; j>=vtx_pointer[u]; j--)
                            for(j=vtx_pointer[u]; j<vtx_pointer[u+1]; j++)
                            {
                                long v = endVertex[j]; // fairness in accessing adjacenty is not helping. why??: no matter how we access, every neighbor will be in the same tree (in serial case). Hence it does not change #iteration. In DFS this may discover shorter augmenting path, but in BFS it does not help.
                                
                                if(flag[v]==0) // avoids unnecessary __sync_fetch_and_or
                                {
                                    if( __sync_fetch_and_or(&flag[v], 1) == 0 )
                                    {
                                        root[v] = curRoot;
                                        parent[v] = u;
                                        
                                        if(mate[v] == -1)
                                        {
                                            leaf[curRoot] = v; //benign race
                                            break;
                                        }
                                        else
                                        {
                                            long next_u = mate[v];
                                            root[next_u] = curRoot;
                                            
                                            if (kbuf < THREAD_BUF_LEN)
                                            {
                                                nbuf[kbuf++] = next_u;
                                            }
                                            else
                                            {
                                                long voff = __sync_fetch_and_add (&QFnextSize, THREAD_BUF_LEN);
                                                for (long vk = 0; vk < THREAD_BUF_LEN; ++vk)
                                                    QFnext[voff + vk] = nbuf[vk];
                                                nbuf[0] = next_u;
                                                kbuf = 1;
                                            }
                                        }
                                    }
                                }
                            }
                            teFwd += j - vtx_pointer[u];
                            
                        }
                    }
                    __sync_fetch_and_add(&eFwd,teFwd);
                    
                }
                else // bottom up BFS
                {
                    long teRev=0;
#pragma omp for
                    for(long v=nrows; v<NV; v++)
                    {
                        if(flag[v]==0)
                        {
                            long j;
                            for(j=vtx_pointer[v+1]-1; j>=vtx_pointer[v]; j--)  // fairness here is important
                                //for(j=vtx_pointer[v]; j<vtx_pointer[v+1]; j++)
                            {
                                long u= endVertex[j];
                                // u must be in the current layer or current layer+1, both cases are fine
                                // if u in layer+1 we are taking a super step in graph traversal (meaning parent and children can be in Queue)
                                if(root[u]!= -1 && leaf[root[u]] == -1) // u is in an active tree
                                {
                                    // Obtaining a parent in the lowest layer gives the shorter augmenting paths
                                    // But it does not reduce size of frontier at any point
                                    // It requires travesing whole adjaency (without the break below), thus costlier
                                    root[v] = root[u];
                                    parent[v] = u;
                                    flag[v] = 1;
                                    if(mate[v] == -1)
                                    {
                                        leaf[root[v]] = v;  // possible benign race
                                    }
                                    else
                                    {
                                        long next_u = mate[v];
                                        root[next_u] = root[v];
                                        
                                        if (kbuf < THREAD_BUF_LEN)
                                        {
                                            nbuf[kbuf++] = next_u;
                                        }
                                        else
                                        {
                                            long voff = __sync_fetch_and_add (&QFnextSize, THREAD_BUF_LEN);
                                            for (long vk = 0; vk < THREAD_BUF_LEN; ++vk)
                                                QFnext[voff + vk] = nbuf[vk];
                                            nbuf[0] = next_u;
                                            kbuf = 1;
                                        }
                                    }
                                    break;
                                    
                                }
                            }
                            //teRev += j - vtx_pointer[v];
                            teRev += vtx_pointer[v+1] - 1 -j;
                        }
                    }
                    __sync_fetch_and_add(&eRev, teRev);
                    
                }
                if(kbuf>0)
                {
                    int64_t voff = __sync_fetch_and_add (&QFnextSize, kbuf);
                    for (long vk = 0; vk < kbuf; ++vk)
                        QFnext[voff + vk] = nbuf[vk];
                }
                
#pragma omp barrier
#pragma omp single
                {
                    long* t;
                    t = QF;
                    QF = QFnext;
                    QFnext = t;
                    QFsize = QFnextSize;
                    QFnextSize = 0;
                    QRsize = QRsize - QFsize; // underestimate
                    if(isTopdownBFS)
                    {
                        timeFwd += omp_get_wtime() - tsLayer;
                        timeFwdAll += omp_get_wtime() - tsLayer;
                    }
                    else
                    {
                        timeRev += omp_get_wtime() - tsLayer;
                        timeRevAll += omp_get_wtime() - tsLayer;
                    }
                    curLayer ++;
                }
            }
        
        }
        double timeBFS = omp_get_wtime() - time_phase;
        
        // ============================================================
        // ---------- Step2: Augment Matching -------------------------
        // ============================================================
        
        double timeAugment_start = omp_get_wtime();
		long nextUnmatchedSize = 0;
#pragma omp parallel
        {
            long kbuf=0, nbuf[THREAD_BUF_LEN];
            long taug_path_len = 0, taug_path_count = 0;
#pragma omp for
            for(long i=0; i<numUnmatchedU; i++)
            {
                long first_u = unmatchedU[i];
                long last_v = leaf[first_u];
                if(last_v != -1)
                {
                    long v = last_v;
                    taug_path_count++;
                    while(v != - 1)
                    {
                        long u = parent[v];
                        long next_v = mate[u];
                        mate[v] = u;
                        mate[u]=v;
                        v = next_v;
                        taug_path_len += 2;
                    }
                }
                else
                {
                    if (kbuf < THREAD_BUF_LEN)
                    {
                        nbuf[kbuf++] = first_u;
                    }
                    else
                    {
                        long voff = __sync_fetch_and_add (&nextUnmatchedSize, THREAD_BUF_LEN);
                        for (long vk = 0; vk < THREAD_BUF_LEN; ++vk)
                            nextUnmatchedU[voff + vk] = nbuf[vk];
                        nbuf[0] = first_u;
                        kbuf = 1;
                    }
                    //nextUnmatchedU[__sync_fetch_and_add(&nextUnmatchedSize,1)] = first_u;
                }
            }
            if(kbuf>0)
            {
                long voff = __sync_fetch_and_add (&nextUnmatchedSize, kbuf);
                for (long vk = 0; vk < kbuf; ++vk)
                    nextUnmatchedU[voff + vk] = nbuf[vk];
            }
            __sync_fetch_and_add(&total_aug_path_len, taug_path_len);
            __sync_fetch_and_add(&total_aug_path_count, taug_path_count);

            
        }
        
        matched = numUnmatchedU - nextUnmatchedSize; // number of new vertices matched in this phase
       
        long* t = unmatchedU;
		unmatchedU = nextUnmatchedU;
        nextUnmatchedU = t;
        long tNumUnmatchedU = numUnmatchedU;
        numUnmatchedU = nextUnmatchedSize;
		timeAugmentAll  += omp_get_wtime() - timeAugment_start ;
        
        
        
        // ===========================================================================
        // statistics: active, inactive & renewable vertices
        // This statistics are use to decide whether to apply tree grafting mechanism
        // ===========================================================================
        
        double timeStat_start = omp_get_wtime();
        long ActiveVtx = 0, InactiveVtx=0, RenewableVtx = 0;
        long step = 100; // smpling, computing statistics for every 100th vertices
        
#pragma omp parallel
        {
            long tActiveVtx = 0, tInactiveVtx=0, tRenewableVtx = 0; //thread private variables
            
#pragma omp for
            for(long u=0; u<nrows; u+=step)
            {
                if(root[u]!=-1 && leaf[root[u]]==-1)
                    tActiveVtx++;
            }
            
#pragma omp for
            for(long v=nrows; v<NV; v+=step)
            {
                if(root[v]==-1)
                    tInactiveVtx ++;
                else if(leaf[root[v]]!=-1)
                    tRenewableVtx ++;
            }
            
            __sync_fetch_and_add(&ActiveVtx,tActiveVtx);
            __sync_fetch_and_add(&RenewableVtx, tRenewableVtx);
            __sync_fetch_and_add(&InactiveVtx,tInactiveVtx);
        }

        //cout << ActiveVtx << " " << RenewableVtx << " " << InactiveVtx << endl;
        timeStatAll+= omp_get_wtime() - timeStat_start;
        
        
        
        // ===========================================================================
        // Step3: reconstruct frontier for the next iteration
        // Either use tree-grafting or build from scratch
        // ===========================================================================
        
        // decide whther tree-grafting is beneficial or not
        bool isGrafting = true;
        double alpha1=5;
        if(ActiveVtx < RenewableVtx/alpha1)
            isGrafting=false;
        
        isGrafting = true;
        double timeGraft_start = omp_get_wtime(); // store the time to reconstruct frontiers
        long eGraft = 0;
        
        if(isGrafting) // use dead-alive scheme
        {
            
            QFsize = 0;
            QRsize = 0;
#pragma omp parallel
            {
                long teGraft = 0;
                long nbuf[THREAD_BUF_LEN];
                long kbuf = 0;
#pragma omp for
                for(long v = nrows; v < NV; v++)
                {
                    //if( root[v]==-1 || leaf[root[v]]!=-1) // consider both dead and unvisited vertices, enble for testing
                    if( root[v]!=-1 && leaf[root[v]]!=-1) // we will not consider unvisited vertices because they can not be part of frontier
                    {
                        // remove v from the dead tree
                        flag[v] = 0;
                        root[v] = -1;
                        
                        // to obtain best result we look for parents in the adjacenty from high to low indices (given that in forward BFS we traverse adjacenty from low to high indices)
                        long j;
                        for(j=vtx_pointer[v+1]-1; j>=vtx_pointer[v]; j--)
                        //for(j=vtx_pointer[v]; j<vtx_pointer[v+1]; j++)
                        {
                            long u = endVertex[j];
                            if(root[u]!= -1 && leaf[root[u]] == -1) // u in an active tree (no augmenting path found in the latest BFS)
                            {
                                if(mate[v] == -1)
                                {
                                    //  leaf[root[v]] = v;  // we can not allow this because then we will try to destroy a tree that has just found an augmenting path (but the augmentation is yet to perform)!
                                    QF[__sync_fetch_and_add(&QFsize,1)] = u; // insert into the frontier again so that we can discover v; we can insert next_u more than once?? No probem: bottom up does not use QF, top down will explore adjacency only once
                                }
                                else
                                {
                                    root[v] = root[u];
                                    parent[v] = u;
                                    flag[v] = 1;
                                    long next_u = mate[v];
                                    root[next_u] = root[v];
                                    //QF[__sync_fetch_and_add(&QFsize,1)] = next_u; // slow version, shared Queue
                                    if (kbuf < THREAD_BUF_LEN)
                                    {
                                        nbuf[kbuf++] = next_u;
                                    }
                                    else
                                    {
                                        long voff = __sync_fetch_and_add (&QFsize, THREAD_BUF_LEN);
                                        for (long vk = 0; vk < THREAD_BUF_LEN; ++vk)
                                            QF[voff + vk] = nbuf[vk];
                                        nbuf[0] = next_u;
                                        kbuf = 1;
                                    }
                                }
                                break;
                            }
                        }
                        //teGraft += j - vtx_pointer[v];
                        teGraft += vtx_pointer[v+1]-1 - j;
                    }
                    
                }
                
                if(kbuf>0)
                {
                    long voff = __sync_fetch_and_add (&QFsize, kbuf);
                    for (long vk = 0; vk < kbuf; ++vk)
                        QF[voff + vk] = nbuf[vk];
                }
                __sync_fetch_and_add(&eGraft, teGraft);
                
            }
            
            QRsize = NV- nrows - QFsize; // speculative reverse frontier size
            
        }
        else // constructs active trees from scratch
        {
#pragma omp parallel for schedule(static) default(shared)
            for(long v = nrows; v < NV; v++)
            {
                if( root[v]!=-1)
                {
                    flag[v]=0;
                    root[v] = -1;
                    //parent[v] = -1;
                }
            }
            
            
#pragma omp parallel for schedule(static) default(shared)
            for(long v = 0; v < nrows; v++)
            {
                if( root[v]!=-1 && leaf[root[v]]==-1)
                {
                    //flag[v]=0;
                    root[v] = -1; // we need this, otherwise in reverse BFS an Y vertex might attached to a zoombie X vertex
                    //parent[v] = -1;
                }
            }
            
            
            QFsize = numUnmatchedU;
#pragma omp parallel for default(shared)
            for(long i=0; i<QFsize; i++)
            {
                long u = unmatchedU[i];
                QF[i] = u;
                root[u] = u;
            }
            
            QRsize = NV-nrows;
        }
        
        
        
        eRevAll += eRev;
        eFwdAll += eFwd;
        eGraftAll += eGraft;
        phaseEdgeVisited += eRev + eFwd + eGraft;
        edgeVisited += eRev + eFwd + eGraft;
       

        
        
	}
    
    return(mate);
}

// helper function used in Karp-Sipser initialization 
void MaximumMatching::findMate(long u, graph* G, long* flag,long* mate, long* degree)
{
	if(__sync_fetch_and_add(&flag[u],1) != 0) return;
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	
	long neighbor_first = edgeStart[u];
	long neighbor_last = edgeStart[u+1];   
	for(long j=neighbor_first; j<neighbor_last; j++) 
	{
		long v = endVertex[j];
		if(__sync_fetch_and_add(&flag[v],1) == 0) // if I can lock then this v node is unmatched
		{
			mate[u] = v;
			mate[v] = u;
			// update degree
			long neighborFirstU = edgeStart[v];
			long neighborLastU = edgeStart[v+1];
			for(long k=neighborFirstU; k< neighborLastU; k++)
			{
				long nextU = endVertex[k];
				if( __sync_fetch_and_add(&degree[nextU],-1) == 2)
				{
					findMate(nextU,G,flag,mate,degree);
				}
				
			}
			break;
		} 
	}
}

long MaximumMatching::KarpSipserInit(SimpleSet &inGraph) {
    long result;
    if(firstKarpSipser) {
        result = KarpSipserInit1();
    } else {
        result = KarpSipserInit2(inGraph);
    }
    firstKarpSipser = false;
    return result;
}


// Multithreaded  Karp-Sipser maximal matching
long MaximumMatching::KarpSipserInit1()
{
	long nrows = G->nrows;
	
	long *endVertex = G->endV;
	long *edgeStart = G->vtx_pointer;
	long nrowsV = G->n;
	long numUnmatchedU = 0;
	
#pragma omp parallel for default(shared) schedule(static)
	for(long i=0; i< nrowsV; i++)
	{
		flag[i] = 0;
		mate[i] = -1;
	}
	
	long degree1Tail = 0;
	long degree1Count = 0;  
	
	
	
	// populate degree and degree1Vtx
#pragma omp parallel for default(shared) schedule(static)//schedule(dynamic)
	for(long u=0; u<nrows; u++)
	{      
		degree[u] = edgeStart[u+1] - edgeStart[u];
		if(degree[u] == 1)
		{
			degree1Vtx[__sync_fetch_and_add(&degree1Count,1)] = u;
			//flag[u] = 1; // means already taken 
		}
	}
	
	
	
#pragma omp parallel for default(shared) //schedule(dynamic,100)
	for(long u=0; u<degree1Count; u++)
	{
		//findMate1(degree1Vtx[u],G,flag,mate,degree,degree1Vtx,&degree1Head);
		findMate(degree1Vtx[u],G,flag,mate,degree);		  
	}
	
	
	// process other vertices 
#pragma omp parallel for default(shared) schedule(dynamic,100)//schedule(dynamic)
	for(long u=0; u<nrows; u++)
	{
		if(flag[u] == 0 && degree[u]>0)
			findMate(u,G,flag,mate,degree);	  
	}
	
#pragma omp parallel for default(shared) 
	for(long u=0; u<nrows; u++)
	{
		
		if(mate[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
		{
			unmatchedU[__sync_fetch_and_add(&numUnmatchedU, 1)] = u;
		}
	}
   	
	return numUnmatchedU;

}

// Multithreaded  Karp-Sipser maximal matching
long MaximumMatching::KarpSipserInit2(SimpleSet &inGraph)
{
    long nrows = G->nrows;
    
    long *endVertex = G->endV;
    long *edgeStart = G->vtx_pointer;
    long nrowsV = G->n;
    long numUnmatchedU = 0;
    
    
#pragma omp parallel for default(shared) schedule(static)
    for(long i=0; i< nrowsV; i++)
    {
        long vertexMate = mate[i];
        if(inGraph.Contains(i % nrows) && vertexMate >= 0 &&  inGraph.Contains(vertexMate % nrows)) {
            flag[i] = 1;
        } else {
            flag[i] = 0;
            mate[i] = -1;
        }
    }
    
    long degree1Tail = 0;
    long degree1Count = 0;  
    
    
    
    // populate degree and degree1Vtx
#pragma omp parallel for default(shared) schedule(static)//schedule(dynamic)
    for(long u=0; u<nrows; u++)
    {      
        if(mate[u] >= 0) continue;
        degree[u] = 0;
        for(long i = edgeStart[u+1]; i < edgeStart[u]; ++i) {
            if(mate[i] >= 0)
                ++degree[u];
        }
        if(degree[u] == 1)
        {
            degree1Vtx[__sync_fetch_and_add(&degree1Count,1)] = u;
            //flag[u] = 1; // means already taken 
        }
    }
    
    
    
#pragma omp parallel for default(shared) //schedule(dynamic,100)
    for(long u=0; u<degree1Count; u++)
    {
        //findMate1(degree1Vtx[u],G,flag,mate,degree,degree1Vtx,&degree1Head);
        findMate(degree1Vtx[u],G,flag,mate,degree);       
    }
    
    
    // process other vertices 
#pragma omp parallel for default(shared) schedule(dynamic,100)//schedule(dynamic)
    for(long u=0; u<nrows; u++)
    {
        if(flag[u] == 0 && degree[u]>0)
            findMate(u,G,flag,mate,degree);   
    }
    
#pragma omp parallel for default(shared) 
    for(long u=0; u<nrows; u++)
    {
        
        if(mate[u] == -1 && (edgeStart[u+1] > edgeStart[u]))
        {
            unmatchedU[__sync_fetch_and_add(&numUnmatchedU, 1)] = u;
        }
    }
    
    return numUnmatchedU;
}

void MaximumMatching::MarkReachableVertices() {
#pragma omp parallel for
	for(int vertex = 0; vertex < G->n; ++vertex) {
		reachableVertices[vertex] = 0;
	}
#pragma omp parallel for
	for(int startVertex = 0; startVertex < G->nrows; ++startVertex) {
		if(mate[startVertex] == -1) {
			assert(reachableVertices[startVertex] == 0);
			reachableVertices[startVertex] = 1;
			long threadId = omp_get_thread_num();
			std::vector<long> &stack = stacks[threadId];
			int top = -1;
			assert(stack.size() == G->n);
			assert(startVertex >= 0);
			assert(startVertex < G->nrows);
			stack[++top] = startVertex;
			
			while (top >= 0 )
			{
				int vertex = stack[top];
				top --;
				assert(vertex < G->nrows);


				for(long i = G->vtx_pointer[vertex]; i < G->vtx_pointer[vertex+1]; i++)
				{
					int neighbor = G->endV[i];
					int matchingMate = mate[neighbor];
					assert(matchingMate != -1);
					assert(matchingMate < G->nrows);
					if(neighbor == matchingMate)
						continue;

					if(__sync_fetch_and_add(&reachableVertices[neighbor],1) == 0) 
					{
						assert(reachableVertices[matchingMate] == 0);
						reachableVertices[matchingMate] += 1;
						assert(reachableVertices[matchingMate] == 1);
						stack[++top] = matchingMate;
						assert(top < G->n);
					}
					
				}				
			}
		}
	}
	assert(IsValidVertexCover());
	assert(CheckVertexCoverAndMatchingSize());
}

// Just for testing
bool MaximumMatching::IsValidVertexCover() {
	for(int vertex = 0; vertex < G->nrows; ++vertex) {
		if(reachableVertices[vertex] != 0) {
			for(long i = G->vtx_pointer[vertex]; i < G->vtx_pointer[vertex+1]; i++) { 
				long neighbor = G->endV[i];
				assert(neighbor >= G->nrows);
				if(reachableVertices[neighbor] == 0)
					return false;
			}
		}
	}
	return true;
}

// Just for testing
bool MaximumMatching::CheckVertexCoverAndMatchingSize() {
	long vertexCoverSize = 0;
	for(long vertex = 0; vertex < G->nrows; ++vertex) {
		if(reachableVertices[vertex] == 0) {
			++vertexCoverSize;
		}
	}

	for(long vertex = G->nrows; vertex < G->n; ++vertex) {
		if(reachableVertices[vertex] > 0) {
			++vertexCoverSize;
		}
	}

	long matchingSize = 0;
	for(long vertex = 0; vertex < G->nrows; ++vertex) {
		if(mate[vertex] != -1) {
			++matchingSize;
		}
	}

	return vertexCoverSize == matchingSize;
}
