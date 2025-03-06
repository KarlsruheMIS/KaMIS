/******************************************************************************
 * kaffpa_interface.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/


#pragma once
#ifdef __cplusplus

extern "C"
{
#endif

const int FAST           = 0;
const int ECO            = 1;
const int STRONG         = 2;
const int FASTSOCIAL     = 3;
const int ECOSOCIAL      = 4;
const int STRONGSOCIAL   = 5;

const int MAPMODE_MULTISECTION = 0;
const int MAPMODE_BISECTION = 1;

// same data structures as in metis 
// edgecut and part are output parameters
// part has to be an array of n ints
void kaffpa(int* n, int* vwgt, int* xadj, 
                   int* adjcwgt, int* adjncy, int* nparts, 
                   double* imbalance, bool suppress_output, int seed, int mode, 
                   int* edgecut, int* part);

// balance constraint on nodes and edges
void kaffpa_balance_NE(int* n, int* vwgt, int* xadj, 
                int* adjcwgt, int* adjncy, int* nparts, 
                double* imbalance,  bool suppress_output, int seed, int mode,
                int* edgecut, int* part);

void node_separator(int* n, int* vwgt, int* xadj, 
                    int* adjcwgt, int* adjncy, int* nparts, 
                    double* imbalance,  bool suppress_output, int seed, int mode,
                    int* num_separator_vertices, int** separator, int* part); 

#ifdef __cplusplus
}
#endif

