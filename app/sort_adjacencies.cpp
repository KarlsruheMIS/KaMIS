/******************************************************************************
 * sort_adjacencies.cpp 
 *
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <algorithm>
#include <unordered_set>

using namespace std;

// This program implements the functions to sort the neighborhood of each node of an input graph
int main(int argn, char **argv) {
        if (argn != 2) {
            std::cout <<  "Usage: sort_adjancencies FILE"  << std::endl;
            exit(0);
        }

        std::string line;
        std::string filename(argv[1]);

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
            std::cerr << "Error opening " << filename << std::endl;
            return 1;
        }

        std::getline(in,line);
        //skip comments
        while (line[0] == '%') {
            std::getline(in, line);
        }

        long nmbNodes;
        long nmbEdges;
        long ew = 0;

        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;
        if(ew != 0) 
                std::cout <<  "sorting adjacencies only supported for unweighted graphs."  << std::endl;

        std::cout <<  nmbNodes << " " << nmbEdges << std::endl;
        while (std::getline(in, line)) {
            if (line[0] == '%') continue;

            std::stringstream ss(line);
            std::vector<long> adjacent_nodes;
                
            long target;
            while (ss >> target) {
                adjacent_nodes.push_back(target);
            }
            std::sort(adjacent_nodes.begin(), adjacent_nodes.end());
            for( unsigned long i = 0; i < adjacent_nodes.size(); i++) {
                    std::cout <<  adjacent_nodes[i] << " ";
            }
            std::cout <<  std::endl;

            if (in.eof()) break;
        }

        return 0;
}

