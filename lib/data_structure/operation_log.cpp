/**
 * operation_log.cpp
 * Purpose: Singleton deque to store the operations performed during ILS.
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

#include "operation_log.h"

operation_log::operation_log() {
    
}

operation_log::~operation_log() {

}

void operation_log::init(unsigned int size) {
    active = false;
    op_log.resize(3*size);
}

void operation_log::report_insert(NodeID x) {
    if (!active) return;
    if (is_full()) printf("Log overflow\n");
    op_log.push_front(x);
}

void operation_log::report_remove(NodeID x) {
    if (!active) return;
    if (is_full()) printf("Log overflow\n");
    op_log.push_front(-x);
}

int operation_log::peek(unsigned int pos) {
    return op_log[pos];
}

unsigned int operation_log::get_size() {
    return op_log.size();
}

int operation_log::unwind() {
    if (is_empty()) printf("Log empty\n");
    int first = op_log.front();
    op_log.pop_front();
    return first;
}

bool operation_log::is_empty() {
    return op_log.empty();
}

bool operation_log::is_full() {
    return (op_log.max_size() == op_log.size());
}

void operation_log::activate() {
    active = true;
}

void operation_log::deactivate() {
    active = false;
}

void operation_log::reset() {
    op_log.clear(); 
}

void operation_log::print() {
    printf("Log size:\t%d", get_size());
    for (unsigned int i = 0; i < get_size(); ++i) {
        int x = peek(i);
        if (x < 0) printf("\t%d", x);
        else printf("\t+%d", x);
    }
    printf("\n");
}

