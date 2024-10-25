/******************************************************************************
 * graph_extractor.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_EXTRACTOR_PDUTVIEF
#define GRAPH_EXTRACTOR_PDUTVIEF

#include "data_structure/graph_access.h"
#include "definitions.h"

class graph_extractor {
        public:
                graph_extractor();
                virtual ~graph_extractor();

                template <typename graph>
                void extract_block(graph & G, 
                                   graph & extracted_block, 
                                   PartitionID block, 
                                   std::vector<NodeID> & mapping);

                template <typename graph>
                void extract_two_blocks(graph & G, 
                                        graph & extracted_block_lhs, 
                                        graph & extracted_block_rhs, 
                                        std::vector<NodeID> & mapping_lhs, 
                                        std::vector<NodeID> & mapping_rhs,
                                        NodeWeight & partition_weight_lhs,
                                        NodeWeight & partition_weight_rhs); 

                template <typename graph>
               void extract_two_blocks_connected(graph & G, 
                                                 std::vector<NodeID> lhs_nodes,
                                                 std::vector<NodeID> rhs_nodes,
                                                 PartitionID lhs, 
                                                 PartitionID rhs,
                                                 graph & pair,
                                                 std::vector<NodeID> & mapping) ;


};

#include "graph_extractor.tpp"
#endif /* end of include guard: GRAPH_EXTRACTOR_PDUTVIEF */
