/**
 * bipartite_weighted_vertex_cover.h
 * Purpose: Create a minimal vertex cover for a given graph.
 *
 *****************************************************************************/

#ifndef _BIPARTITE_WEIGHTED_VERTEX_COVER_H_
#define _BIPARTITE_WEIGHTED_VERTEX_COVER_H_

#include <vector>

#include "mmwis_config.h"
#include "definitions.h"
#include "data_structure/flow_graph.h"
#include "mmwis_graph_access.h"

class bipartite_weighted_vertex_cover {
    public:
        /**
         * Default Constructor.
         */
        bipartite_weighted_vertex_cover();

        /**
         * Default Destructor.
         */
        virtual ~bipartite_weighted_vertex_cover();

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
        void max_flow_cover(mmwis::graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, std::vector<NodeID> & vertex_cover);

    private:
        /**
         * Creates a bipartition for the given graph.
         *
         * @param G Graph representation.
         * @param bipartite The resulting bipartite graph.
         */
        void create_bipartite(mmwis::graph_access & G, std::vector<NodeID> & lhs, std::vector<NodeID> & rhs, mmwis::graph_access & bipartite, std::vector<NodeID> & bipartite_mapping);

        /**
         * Creates a flow network for bipartition.
         *
         * @param bipartite The bipartite graph.
         */
        void create_flow_network(mmwis::graph_access & bipartite, flow_graph & flow_network, NodeID & source, NodeID & sink);
};

#endif
