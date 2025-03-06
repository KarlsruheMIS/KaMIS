/**
 * bucket_array.cpp
 * Purpose: A bucket priority queue used for building initial maximum independent sets.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "bucket_array.h"

#include "random_functions.h"

bucket_array::bucket_array(unsigned int number_of_nodes) {
    size = number_of_nodes;
    init();
}

bucket_array::~bucket_array() {
    delete [] position;
    delete [] array;
    delete [] value;
    delete [] first;
}

void bucket_array::init() {
    position = new unsigned int[size]; 
    array = new NodeID[size];
    value = new int[size];
    first = new unsigned int[size + 1];

    for (unsigned int i = 0; i < size; ++i) {
        position[i] = i;
        array[i] = i;
        value[i] = 0;
        first[i] = size;
    }
    first[0] = 0;
    first[size] = size;
    non_negative = 0;
}

void bucket_array::increment(NodeID node, int value) {
    for (; value > 0; --value) {
        increment(node);
    }
}

void bucket_array::increment(NodeID node) {
    int node_value = value[node];
    int last = first[node_value+1] - 1; 
    
    swap (last, position[node]);

    first[node_value+1]--; 
    value[node]++;
}

void bucket_array::decrement(NodeID node) {
    int node_value = value[node];

    swap(first[node_value], position[node]);
    first[node_value]++;
    value[node]--;

    if (value[node] < 0) non_negative++;
}

bool bucket_array::contains(NodeID node) {
    return (value[node] >= 0); 
}

void bucket_array::remove(NodeID node) {
    while (value[node] >= 0) {
        decrement(node);
    }
}

NodeID bucket_array::pickSmallest() {
    if (non_negative >= size) return -1;
    int block = value[array[non_negative]];
    unsigned int pick = random_functions::nextInt(first[block], first[block + 1] - 1);

    return array[pick];
}

void bucket_array::swap(int first_pos, int second_pos) {
    NodeID first = array[first_pos];
    NodeID second = array[second_pos];

    array[first_pos] = second;
    array[second_pos] = first;
    position[first] = second_pos;
    position[second] = first_pos;
}

bool bucket_array::is_empty(int block) {
    return (first[block+1] == first[block]); 
}

