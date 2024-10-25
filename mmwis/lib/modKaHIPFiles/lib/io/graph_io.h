/******************************************************************************
 * graph_io.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPHIO_H_
#define GRAPHIO_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "definitions.h"
#include "data_structure/graph_access.h"

class graph_io
{
public:
        graph_io();
        virtual ~graph_io();

        template <typename graph>
        static int readGraphWeighted(graph &G, const std::string &filename);

        template <typename graph>
        static int writeGraphWeighted(graph &G, const std::string &filename);

        template <typename graph>
        static int writeGraph(graph &G, const std::string &filename);

        template <typename graph>
        static int readPartition(graph &G, const std::string &filename);

        template <typename graph>
        static void writePartition(graph &G, const std::string &filename);

        template <typename vectortype>
        static void writeVector(std::vector<vectortype> &vec, const std::string &filename);

        template <typename vectortype>
        static void readVector(std::vector<vectortype> &vec, const std::string &filename);
};
#include "graph_io.tpp"
/*
graph_io::graph_io()
{
}

graph_io::~graph_io()
{
}

template <typename vectortype>
void graph_io::writeVector(std::vector<vectortype> &vec, const std::string &filename)
{
        std::ofstream f(filename.c_str());
        for (unsigned i = 0; i < vec.size(); ++i)
        {
                f << vec[i] << std::endl;
        }

        f.close();
}

template <typename vectortype>
void graph_io::readVector(std::vector<vectortype> &vec, const std::string &filename)
{

        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in)
        {
                std::cerr << "Error opening vectorfile" << filename << std::endl;
                return;
        }

        unsigned pos = 0;
        std::getline(in, line);
        while (!in.eof())
        {
                if (line[0] == '%')
                { // Comment
                        continue;
                }

                vectortype value = (vectortype)atof(line.c_str());
                vec[pos++] = value;
                std::getline(in, line);
        }

        in.close();
}

template <typename graph>
int graph_io::writeGraphWeighted(graph &G, const std::string &filename)
{
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() << " " << G.number_of_edges() / 2 << " 11" << std::endl;

        forall_nodes(G, node)
        {
                f << G.getNodeWeight(node);
                forall_out_edges(G, e, node)
                {
                        f << " " << (G.getEdgeTarget(e) + 1) << " " << G.getEdgeWeight(e);
                }
                endfor
                        f
                    << "\n";
        }
        endfor

            f.close();
        return 0;
}

template <typename graph>
int graph_io::writeGraph(graph &G, const std::string &filename)
{
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() << " " << G.number_of_edges() / 2 << std::endl;

        forall_nodes(G, node)
        {
                forall_out_edges(G, e, node)
                {
                        f << (G.getEdgeTarget(e) + 1) << " ";
                }
                endfor
                        f
                    << "\n";
        }
        endfor

            f.close();
        return 0;
}

template <typename graph>
int graph_io::readPartition(graph &G, const std::string &filename)
{
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in)
        {
                std::cerr << "Error opening file" << filename << std::endl;
                return 1;
        }

        PartitionID max = 0;
        forall_nodes(G, node)
        {
                // fetch current line
                std::getline(in, line);
                if (line[0] == '%')
                { // Comment
                        node--;
                        continue;
                }

                // in this line we find the block of Node node
                G.setPartitionIndex(node, (PartitionID)atol(line.c_str()));

                if (G.getPartitionIndex(node) > max)
                        max = G.getPartitionIndex(node);
        }
        endfor

            G.set_partition_count(max + 1);
        in.close();

        return 0;
}

template <typename graph>
int graph_io::readGraphWeighted(graph &G, const std::string &filename)
{
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in)
        {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        long nmbNodes;
        long nmbEdges;

        std::getline(in, line);
        // skip comments
        while (line[0] == '%')
        {
                std::getline(in, line);
        }

        int ew = 0;
        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        if (2 * nmbEdges > std::numeric_limits<int>::max() || nmbNodes > std::numeric_limits<int>::max())
        {
                std::cerr << "The graph is too large. Currently only 32bit supported!" << std::endl;
                exit(0);
        }

        bool read_ew = false;
        bool read_nw = false;

        if (ew == 1)
        {
                read_ew = true;
        }
        else if (ew == 11)
        {
                read_ew = true;
                read_nw = true;
        }
        else if (ew == 10)
        {
                read_nw = true;
        }
        nmbEdges *= 2; // since we have forward and backward edges

        NodeID node_counter = 0;
        EdgeID edge_counter = 0;
        long long total_nodeweight = 0;

        G.start_construction(nmbNodes, nmbEdges);

        while (std::getline(in, line))
        {

                if (line[0] == '%')
                { // a comment in the file
                        continue;
                }

                NodeID node = G.new_node();
                node_counter++;
                G.setPartitionIndex(node, 0);

                std::stringstream ss(line);

                NodeWeight weight = 1;
                if (read_nw)
                {
                        ss >> weight;
                        total_nodeweight += weight;
                        if (total_nodeweight > (long long)std::numeric_limits<NodeWeight>::max())
                        {
                                std::cerr << "The sum of the node weights is too large (it exceeds the node weight type)." << std::endl;
                                std::cerr << "Currently not supported. Please scale your node weights." << std::endl;
                                exit(0);
                        }
                }
                G.setNodeWeight(node, weight);

                NodeID target;
                while (ss >> target)
                {
                        // check for self-loops
                        if (target - 1 == node)
                        {
                                std::cerr << "The graph file contains self-loops. This is not supported. Please remove them from the file." << std::endl;
                        }

                        EdgeWeight edge_weight = 1;
                        if (read_ew)
                        {
                                ss >> edge_weight;
                        }
                        edge_counter++;
                        EdgeID e = G.new_edge(node, target - 1);

                        G.setEdgeWeight(e, edge_weight);
                }

                if (in.eof())
                {
                        break;
                }
        }

        if (edge_counter != (EdgeID)nmbEdges)
        {
                std::cerr << "number of specified edges mismatch" << std::endl;
                std::cerr << edge_counter << " " << nmbEdges << std::endl;
                exit(0);
        }

        if (node_counter != (NodeID)nmbNodes)
        {
                std::cerr << "number of specified nodes mismatch" << std::endl;
                std::cerr << node_counter << " " << nmbNodes << std::endl;
                exit(0);
        }

        G.finish_construction();
        return 0;
}

template <typename graph>
void graph_io::writePartition(graph &G, const std::string &filename)
{
        std::ofstream f(filename.c_str());
        std::cout << "writing partition to " << filename << " ... " << std::endl;

        forall_nodes(G, node)
        {
                f << G.getPartitionIndex(node) << "\n";
        }
        endfor

            f.close();
}
*/
#endif /*GRAPHIO_H_*/
