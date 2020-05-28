/**
 * candidate_list.cpp
 * Purpose: Represents the list used to store the candidates for the local search.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#include "candidate_list.h"

#include "macros_assertions.h"
#include "random_functions.h"

candidate_list::candidate_list() { count = 0; }

candidate_list::~candidate_list() {}

void candidate_list::init(unsigned int size) {
        nodes.clear();
        position.clear();
        nodes.resize(size);
        position.resize(size);
        for (unsigned int i = 0; i < size; ++i) {
                position[i] = -1;
                nodes[i] = 0;
        }
        count = 0;
}

void candidate_list::random_permute() {
        if (count <= 1) return;

        for (unsigned int i = 0; i < count; ++i) {
                unsigned int j = random_functions::nextInt(i, count - 1);
                if (i != j) swap(i, j);
        }
}

void candidate_list::swap(unsigned int i, unsigned int j) {
        NodeID node_i = nodes[i];
        NodeID node_j = nodes[j];
        nodes[i] = node_j;
        nodes[j] = node_i;
        position[node_j] = i;
        position[node_i] = j;
}

NodeID candidate_list::remove_random() {
        int random = random_functions::nextInt(0, count - 1);
        return remove_position(random);
}

NodeID candidate_list::pick_random() {
        int rand = random_functions::nextInt(0, count - 1);
        return pick(rand);
}

NodeID candidate_list::pick(unsigned int pos) { return nodes[pos]; }

void candidate_list::insert(NodeID node) {
        position[node] = count;
        nodes[count] = node;
        count++;
}

bool candidate_list::contains(NodeID node) {
        bool in = false;
        if ((unsigned)position[node] < count && position[node] != -1) in = true;
        return in;
}

bool candidate_list::is_empty() { return (count == 0); }

void candidate_list::reset() {
        while (!is_empty()) remove_position(0);
        ASSERT_EQ(count, 0);
}

void candidate_list::remove(NodeID node) { remove_position(position[node]); }

NodeID candidate_list::remove_position(unsigned int pos) {
        ASSERT_LT(pos, count);

        NodeID node = nodes[pos];
        if (pos < count - 1) {
                NodeID last = nodes[count - 1];
                position[last] = pos;
                nodes[pos] = last;
        }
        position[node] = -1;
        count--;

        return node;
}

void candidate_list::print() {
        printf("Candidates\n");
        for (unsigned int i = 0; i < count; ++i) printf("Node: %d\n", nodes[i]);
}
