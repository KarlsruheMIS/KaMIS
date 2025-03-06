//
// Created by alex on 08.01.20.
//

#ifndef _struction_COMPONENTS_VERTEX_MARKER_H
#define _struction_COMPONENTS_VERTEX_MARKER_H

#include <definitions.h>
#include <data_structure/sized_vector.h>
#include "fast_set.h"

class struction_vertex_marker {
private:

public:

    sized_vector<NodeID> current;
    sized_vector<NodeID> next;
    fast_set added_vertices;

    struction_vertex_marker() : struction_vertex_marker(10) {}
	struction_vertex_marker(size_t size) : current(size), next(size), added_vertices(size) {};

	void add(NodeID vertex) {
		if (!added_vertices.get(vertex)) {
			next.push_back(vertex);
			added_vertices.add(vertex);
		}
	}

	void get_next() {
		//if (next.size() != 0) {
			current.swap(next);
			clear_next();
		//}
	}

	void clear_next() {
		next.clear();
		added_vertices.clear();
	}

	void fill_current_ascending(size_t n) {
		current.clear();
		for (size_t i = 0; i < n; i++) {
			current.push_back(i);
		}
	}

	NodeID current_vertex(size_t index) {
		return current[index];
	}

	size_t current_size() {
		return current.size();
	}

    void resize(size_t size) {
        added_vertices.resize(size);
        current.resize(size);
        next.resize(size);
    }

    sized_vector<NodeID> &current_vec() {
	    return current;
	}
};

#endif //_struction_COMPONENTS_VERTEX_MARKER_H
