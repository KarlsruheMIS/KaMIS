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

#ifndef PROFILING_HELPER_H
#define PROFILING_HELPER_H

#include "SparseArraySet.h"

#include <sys/time.h>
#include <unordered_map>

#ifdef PROFILING
struct ProfilingHelper_t {
	std::vector<struct timeval> startTimersPerPartition;
	std::vector<int> currentVertexPartition;
	std::vector<struct timeval> startTimersUpdateNeighborhoodPerPartition;
	std::vector<int> currentVertexUpdateNeighborhoodPartition;

	std::vector<SparseArraySet> *neighborsPtr;

	std::vector<std::unordered_map<int, std::vector<unsigned int>>> timesUpdateNeighborhoodsPerDegree; // partition, degree

    std::vector<std::vector<unsigned int>> timesUnsuccessfulFoldDegree; //partition
    std::vector<std::unordered_map<int, std::vector<unsigned int>>> timesUnsuccessfulFoldAdjacentPerNeighborDegree; // partition, first neighbor degree
    std::vector<std::vector<unsigned int>> timesUnsuccessfulFoldWrongPartition; // partition,
    std::vector<std::unordered_map<int, std::vector<unsigned int>>> timesSuccessfulFoldPerTwoNeighborhoodSize; // partition, size of 2-neighborhood

    std::vector<std::vector<unsigned int>> timesUnsuccessfulIsolatedCliquePartition;
    std::vector<std::unordered_map<int, std::vector<unsigned int>>> timesUnsuccessfulIsolatedCliqueDegreePerDegree; // partition, degree
    std::vector<std::unordered_map<int, std::vector<unsigned int>>> timesUnsuccessfulIsolatedCliqueNoCliquePerDegree; // partition, degree
    std::vector<std::unordered_map<int, std::vector<unsigned int>>> timesSuccessfulIsolatedCliquePerDegree; // partition, degree
};

namespace {

	int twoNeighbordhoodSize(ProfilingHelper_t *profilingHelper, int const vertex) {
		int size = (*profilingHelper->neighborsPtr)[vertex].Size();
		for(int neighbor : (*profilingHelper->neighborsPtr)[vertex]) {
	        size += (*profilingHelper->neighborsPtr)[neighbor].Size();
	    }
	    return size;
	}

	void addNewTime(std::unordered_map<int, std::vector<unsigned int>> &map, int const index, unsigned int const time) {
		map[index].push_back(time);
	}

	void printVectorMapVector(std::vector<std::unordered_map<int, std::vector<unsigned int>>> &vec, const char name[]) {
		std::cout << "#########################################################" << std::endl;
		for(int partition = 0; partition < vec.size(); partition++) {
			std::cout << "---------------------------------------------------------" << std::endl;
			for(std::pair<const int, std::vector<unsigned int>> element: vec[partition]) {
				int size = element.first;
				std::vector<unsigned int> times = element.second;
				std::cout << name << "[" <<partition << "]: " << size << " -";
				for(unsigned int time : times) {
					std::cout << " " << time;
				}
				std::cout << std::endl;
			}
		}
	}

	void printTimeVector(std::vector<std::vector<unsigned int>> &vec, const char name[]) {
		std::cout << "#########################################################" << std::endl;
		for(int partition = 0; partition < vec.size(); partition++) {
			std::cout << name << "[" <<partition << "]:";
			for(unsigned int time : vec[partition]) {
				std::cout << " " << time;
			}
			std::cout << std::endl;
		}
	}

	unsigned int getTime(struct timeval *tBefore) {
		struct timeval tAfter;
		gettimeofday(&tAfter, NULL);
		return ((tAfter.tv_sec - tBefore->tv_sec) * 1000000L + tAfter.tv_usec) - tBefore->tv_usec;
	}
}

void profilingInit(ProfilingHelper_t *profilingHelper, std::vector<SparseArraySet> *neighbors, int const numPartitions)  {
	profilingHelper->neighborsPtr = neighbors;
	profilingHelper->timesUpdateNeighborhoodsPerDegree = std::vector<std::unordered_map<int, std::vector<unsigned int>>>(numPartitions);
	profilingHelper->timesUnsuccessfulFoldDegree = std::vector<std::vector<unsigned int>>(numPartitions);
	profilingHelper->timesUnsuccessfulFoldAdjacentPerNeighborDegree = std::vector<std::unordered_map<int, std::vector<unsigned int>>>(numPartitions);
	profilingHelper->timesUnsuccessfulFoldWrongPartition = std::vector<std::vector<unsigned int>>(numPartitions);
	profilingHelper->timesSuccessfulFoldPerTwoNeighborhoodSize = std::vector<std::unordered_map<int, std::vector<unsigned int>>>(numPartitions);
	profilingHelper->timesUnsuccessfulIsolatedCliquePartition = std::vector<std::vector<unsigned int>>(numPartitions);
	profilingHelper->timesUnsuccessfulIsolatedCliqueDegreePerDegree = std::vector<std::unordered_map<int, std::vector<unsigned int>>>(numPartitions);
	profilingHelper->timesUnsuccessfulIsolatedCliqueNoCliquePerDegree = std::vector<std::unordered_map<int, std::vector<unsigned int>>>(numPartitions);
	profilingHelper->timesSuccessfulIsolatedCliquePerDegree = std::vector<std::unordered_map<int, std::vector<unsigned int>>>(numPartitions);
	profilingHelper->startTimersPerPartition = std::vector<struct timeval>(numPartitions);
	profilingHelper->currentVertexPartition = std::vector<int>(numPartitions);
	profilingHelper->startTimersUpdateNeighborhoodPerPartition = std::vector<struct timeval>(numPartitions);
	profilingHelper->currentVertexUpdateNeighborhoodPartition = std::vector<int>(numPartitions);
}

void profilingStartClock(ProfilingHelper_t *profilingHelper, int const partition, int const vertex) {
	profilingHelper->currentVertexPartition[partition] = vertex;
	gettimeofday(&(profilingHelper->startTimersPerPartition[partition]), NULL);
}

void profilingStartClockUpdateNeighborhood(ProfilingHelper_t *profilingHelper, int const partition, int const vertex) {
	profilingHelper->currentVertexUpdateNeighborhoodPartition[partition] = vertex;
	gettimeofday(&(profilingHelper->startTimersUpdateNeighborhoodPerPartition[partition]), NULL);
}

void profilingAddTimeUpdateNeighborhood(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersUpdateNeighborhoodPerPartition[partition]));
	int vertex = profilingHelper->currentVertexUpdateNeighborhoodPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesUpdateNeighborhoodsPerDegree[partition], degree, time);
}

void profilingAddTimeUnsuccessfulFoldDegree(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	profilingHelper->timesUnsuccessfulFoldDegree[partition].push_back(time);
}

void profilingAddTimeUnsuccessfulFoldAdjacent(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	int vertex = profilingHelper->currentVertexPartition[partition];
	int neighbor = (*profilingHelper->neighborsPtr)[vertex][0];
	int neighborDegree = (*profilingHelper->neighborsPtr)[neighbor].Size();
	addNewTime(profilingHelper->timesUnsuccessfulFoldAdjacentPerNeighborDegree[partition], neighborDegree, time);
}

void profilingAddTimeUnsuccessfulFoldWrongPartition(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	int vertex = profilingHelper->currentVertexPartition[partition];
	profilingHelper->timesUnsuccessfulFoldWrongPartition[partition].push_back(time);
}

void profilingAddTimeSuccessfulFold(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	int vertex = profilingHelper->currentVertexPartition[partition];
	int NeighbordhoodSize = twoNeighbordhoodSize(profilingHelper, vertex);
	addNewTime(profilingHelper->timesSuccessfulFoldPerTwoNeighborhoodSize[partition], NeighbordhoodSize, time);
}

void profilingAddTimeUnsuccessfulIsolatedCliqueDegree(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	int vertex = profilingHelper->currentVertexPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesUnsuccessfulIsolatedCliqueDegreePerDegree[partition], degree, time);
}

void profilingAddTimeUnsuccessfulIsolatedCliquePartition(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	profilingHelper->timesUnsuccessfulIsolatedCliquePartition[partition].push_back(time);
}

void profilingAddTimeUnsuccessfulIsolatedCliqueNoClique(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	int vertex = profilingHelper->currentVertexPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesUnsuccessfulIsolatedCliqueNoCliquePerDegree[partition], degree, time);
}

void profilingAddTimeSuccessfulIsolatedClique(ProfilingHelper_t *profilingHelper, int const partition) {
	unsigned int time = getTime(&(profilingHelper->startTimersPerPartition[partition]));
	int vertex = profilingHelper->currentVertexPartition[partition];
	int degree = (*profilingHelper->neighborsPtr)[vertex].Size();
	addNewTime(profilingHelper->timesSuccessfulIsolatedCliquePerDegree[partition], degree, time);
}

void profilingPrint(ProfilingHelper_t *profilingHelper) {
	std::cout << "#########################################################" << std::endl;
	std::cout << std::endl << "Profiling:" << std::endl;

	printVectorMapVector(profilingHelper->timesUpdateNeighborhoodsPerDegree, "UpdateNeighborhoodsPerDegree");

	printTimeVector(profilingHelper->timesUnsuccessfulFoldDegree, "UnsuccessfulFoldDegree");
	printVectorMapVector(profilingHelper->timesUnsuccessfulFoldAdjacentPerNeighborDegree, "UnsuccessfulFoldAdjacentPerNeighborDegree");
	printTimeVector(profilingHelper->timesUnsuccessfulFoldWrongPartition, "UnsuccessfulFoldWrongPartition");
	printVectorMapVector(profilingHelper->timesSuccessfulFoldPerTwoNeighborhoodSize, "SuccessfulFoldPerTwoNeighborhoodSize");

	printTimeVector(profilingHelper->timesUnsuccessfulIsolatedCliquePartition, "UnsuccessfulIsolatedCliquePartition");
	printVectorMapVector(profilingHelper->timesUnsuccessfulIsolatedCliqueDegreePerDegree, "UnsuccessfulIsolatedCliqueDegreePerDegree");
	printVectorMapVector(profilingHelper->timesUnsuccessfulIsolatedCliqueNoCliquePerDegree, "UnsuccessfulIsolatedCliqueNoCliquePerDegree");
	printVectorMapVector(profilingHelper->timesSuccessfulIsolatedCliquePerDegree, "SuccessfulIsolatedCliquePerDegree");

	std::cout << "#########################################################" << std::endl;
}
#else
struct ProfilingHelper_t {
	char nothing[0];
};
#define profilingInit (void)sizeof
#define profilingStartClock (void)sizeof
#define profilingStartClockUpdateNeighborhood (void)sizeof
#define profilingAddTimeUpdateNeighborhood (void)sizeof
#define profilingAddTimeUnsuccessfulFoldDegree (void)sizeof
#define profilingAddTimeUnsuccessfulFoldAdjacent (void)sizeof
#define profilingAddTimeUnsuccessfulFoldWrongPartition (void)sizeof
#define profilingAddTimeSuccessfulFold (void)sizeof
#define profilingAddTimeUnsuccessfulIsolatedCliqueDegree (void)sizeof
#define profilingAddTimeUnsuccessfulIsolatedCliquePartition (void)sizeof
#define profilingAddTimeUnsuccessfulIsolatedCliqueNoClique (void)sizeof
#define profilingAddTimeSuccessfulIsolatedClique (void)sizeof
#define profilingPrint (void)sizeof
#endif
#endif //PROFILING_HELPER_H
