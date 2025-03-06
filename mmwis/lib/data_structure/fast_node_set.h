/******************************************************************************
 * fast_node_set.h
 *
 *****************************************************************************/

#pragma once


#include <vector>
#include "definitions.h"

class fast_node_set {

	std::vector<NodeID> used;
	NodeID uid;

  public:
	fast_node_set(int const n) : used(n, 0), uid(1)
    {}

    int size() {
        return used.size();
    }

    int capacity() {
        return used.capacity();
    }

	void clear() {
		uid++;
		if (uid < 0) {
			for (unsigned int i = 0; i < used.size(); i++) used[i] = 0;
			uid = 1;
		}
	}

	bool add(int i) {
		assert(i < used.size() && "fast_node_set::add: index out of bounds");
		bool const res(used[i] != uid);
		used[i] = uid;
		return res;
	}

	void remove(int i) {
		assert(i < used.size() && "fast_node_set::remove: index out of bounds");
		used[i] = uid - 1;
	}

	bool get(int i) const {
		assert(i < used.size() && "fast_node_set::get: index out of bounds");
		return (used[i] == uid);
	}

    void print() {
        std::cout << "Fast node set: ";
        for (NodeID n = 0; n< used.size(); n++) std::cout << used[n] << " ";
        std::cout << std::endl;
    }

};

