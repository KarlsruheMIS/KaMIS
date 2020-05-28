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

struct WeightedEdge {
        long target;
        long weight;
        bool operator < (const WeightedEdge& rhs) const {
                return target < rhs.target;
        }
};

// This program implements the functions to sort the neighborhood of each node of an input graph
int main(int argn, char **argv) {
        if (argn != 2) {
            std::cout <<  "Usage: sort_adjacencies FILE"  << std::endl;
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
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

        int ew = 0;
        long nmbNodes;
        long nmbEdges;
        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        std::cout <<  nmbNodes << " " << nmbEdges << " " << ew << std::endl;
        bool read_ew = false;
        bool read_nw = false;

        if(ew == 1) {
                read_ew = true;
        } else if (ew == 11) {
                read_ew = true;
                read_nw = true;
        } else if (ew == 10) {
                read_nw = true;
        }
        
        while(  std::getline(in, line)) {
       
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                std::stringstream ss(line);

                long weight = 1;
                if( read_nw ) {
                        ss >> weight;
                        std::cout <<  weight << " ";
                }

                long target;
                std::vector< WeightedEdge > edges_with_weight;
                while( ss >> target ) {
                        WeightedEdge e;
                        e.target = target;
                        e.weight = 1;
                        long edge_weight = 1;
                        if( read_ew ) {
                                ss >> edge_weight;
                                e.weight = edge_weight;
                        }
                        edges_with_weight.push_back(e);
                }
                std::sort(edges_with_weight.begin(), edges_with_weight.end());
                if(read_ew) {
                        for( unsigned long i = 0; i < edges_with_weight.size(); i++) {
                                std::cout <<  edges_with_weight[i].target << " " << edges_with_weight[i].weight << " " ;
                        }
                
                } else {
                        for( unsigned long i = 0; i < edges_with_weight.size(); i++) {
                                std::cout <<  edges_with_weight[i].target << " ";
                        }
                }
                std::cout << std::endl;

                if(in.eof()) {
                        break;
                }
        }

        return 0;
}

