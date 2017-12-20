/******************************************************************************
 * fast_set.h
 *
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
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

#ifndef FAST_SET_H
#define FAST_SET_H

#include <vector>

class fast_set {
	
	std::vector<int> used;
	int uid;

public:
	fast_set(int const n) : used(n, 0), uid(1)
    { }
	
	void clear() {
		uid++;
		if (uid < 0) {
			for (unsigned int i = 0; i < used.size(); i++) used[i] = 0;
			uid = 1;
		}
	}
	
	bool add(int i) {
		bool const res(used[i] != uid);
		used[i] = uid;
		return res;
	}
	
	void remove(int i) {
		used[i] = uid - 1;
	}
	
	bool get(int i) {
		return (used[i] == uid);
	}
};

#endif // FAST_SET_H
