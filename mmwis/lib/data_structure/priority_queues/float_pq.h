/**
 * bucket_array.h
 * Purpose: A bucket priority queue used for building initial maximum independent sets.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _FLOAT_ARRAY_H_
#define _FLOAT_ARRAY_H_

#include "definitions.h"
struct NodeID_value_comparison 
{
    bool operator()(std::pair<NodeID, float> & a, std::pair<NodeID, float> & b)
    {
        return a.second<b.second;
    }
};
#endif

