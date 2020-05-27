/**
 * bipartite_vertex_cover.h
 * Purpose: Create a minimal vertex cover for a given graph.
 *
 *****************************************************************************/

#ifndef _BIPARTITE_VERTEX_COVER_H_
#define _BIPARTITE_VERTEX_COVER_H_

#include <vector>

#include "mis_config.h"
#include "definitions.h"
#include "data_structure/graph_access.h"

class bipartite_vertex_cover {
    public:
        /**
         * Default Constructor.
         */
        bipartite_vertex_cover();

        /**
         * Default Destructor.
         */
        virtual ~bipartite_vertex_cover();

        /**
         * Creates a minimal vertex cover for the given graph.
         * A maximum matching is created using the Hopcroft-Karp algorithm.
         * In parallel a vertex cover is calculated using Koenig's theorem.
         *
         * @param G Graph representation.
         * @param vertex_cover The resulting vertex cover.
         * @param lhs Nodes on the left side of the bipartite graph.
         * @param rhs Nodes on the right side of the bipartite graph.
         */
        void hopcroft_cover(graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, std::vector<NodeID> & vertex_cover);

    private:
        /**
         * Hopcroft-Karp algorithm using Koenig's theorem to calculate a vertex cover.
         *
         * @param G Graph representation.
         * @param edge_matching Maximal matching.
         * @param T Array used for the calculation of the vertex cover.
         */
        int hopcroft(graph_access & G, Matching & edge_matching, std::vector<bool> & T);

        /**
         * DFS used in the Hopcroft-Karp algorithm.
         *
         * @param G Graph representation.
         * @param edge_matching Maximal matching.
         * @param dist Distance array used for determining augmented paths.
         * @param T Array used for the calculation of the vertex cover.
         * @param v Node for which an augmented path should be found.
         * @return True if the DFS was successful.
         */
        bool hopcroft_dfs(graph_access & G, Matching & edge_matching, std::vector<int> & dist, std::vector<bool> & T, NodeID v);

        /**
         * BFS used in the Hopcroft-Karp algorithm.
         *
         * @param G Graph representation.
         * @param edge_matching Maximal matching.
         * @param dist Distance array used for determining augmented paths.
         * @return True if the BFS was successful.
         */
        bool hopcroft_bfs(graph_access & G, Matching & edge_matching, std::vector<int> & dist);

        /**
         * Koenig's theorem for finding a minimal vertex cover.
         *
         * @param G Graph representation.
         * @param edge_matching Maximal matching.
         * @param T Array used for the calculation of the vertex cover.
         */
        void hopcroft_koenig(graph_access & G, Matching & edge_matching, std::vector<bool> & T);

        /**
         * Creates a bipartition for the given graph.
         *
         * @param G Graph representation.
         * @param bipartite The resulting bipartite graph.
         */
        void create_bipartite(graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, graph_access & bipartite, std::vector<NodeID> & bipartite_mapping);
};

#endif
