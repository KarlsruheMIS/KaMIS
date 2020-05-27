/**
 * operation_log.cpp
 * Purpose: Singleton deque to store the operations performed during ILS.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
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

