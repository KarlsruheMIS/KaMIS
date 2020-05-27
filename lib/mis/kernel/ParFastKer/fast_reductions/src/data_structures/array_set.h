/**
 * array_set.h
 * Purpose: Dynamic array representation for sets with a given range of elements.
 *          Lookups can be performed in constant time using a dedicated array.
 *
 * The original code was kindly provided by Darren Strash.
 *
 *****************************************************************************/

#ifndef _ARRAY_SET_H
#define _ARRAY_SET_H

#include <set>
#include <vector>
#include <iostream>
#include <cassert>
#include <utility>

class array_set {
    public:
        /**
         * Constructor.
         * Initializes the array with the given size.
         *
         * @param size Size of the array.
         */
        array_set(unsigned int size) : lookup(size, -1), elements(size, -1), first(0), last(-1) { }

        /**
         * Default Constructor.
         */
        array_set() : lookup(), elements(), first(0), last(-1) { }

        /**
         * Default Destructor.
         */
        ~array_set() { }

        /**
         * Resizes the array with a given size.
         *
         * @param size Size of the array.
         */
        void resize(unsigned int size) {
            lookup.resize(size, -1);
            elements.resize(size, -1);
        }

        /**
         * Initializes the array with a given size.
         *
         * @param size Size of the array.
         */
        void init(unsigned int size) {
            lookup.resize(size, -1);
            elements.resize(size, -1);
            first = 0;
            last = -1;
        }

        /**
         * Initializes the array with the size of a given adjacency array.
         * Inserts all neighbors for the given node.
         *
         * @param adj Adjacency array of the array.
         * @param node Node for which neighbors will be added.
         */
        void init_from_adj(std::vector<std::vector<int>> & adj, int node) {
            resize(adj.size());
            for (int neighbor : adj[node]) {
                insert(neighbor);
            }
        }

        /**
         * Check if the set contains a given node.
         * Lookup can be performed in constant time.
         *
         * @param x Node to search for.
         * @return True, if the node exists in the array.
         *         False, otherwise.
         */
        bool contains(int const x) const {
            if (x < 0 || x >= static_cast<int>(lookup.size())) return false;
            int const location_x = lookup[x];
            return location_x >= first && location_x <= last;
        }

        /**
         * Insert a given node to the set.
         * If the node is already present, do nothing.
         *
         * @param x Node to insert.
         */
        void insert(int const x) {
            if (contains(x)) return;
            if (elements[last + 1] != -1 && !contains(elements[last + 1])) 
                lookup[elements[last + 1]] = -1; 
            last++;
            lookup[x] = last;
            elements[last] = x;
        }

        /**
         * Remove a given node from the set.
         * If the node is not present, do nothing.
         *
         * @param x Node to remove.
         */
        void remove(int x) {
            if (!contains(x)) return;
            int location_x(lookup[x]);
            elements[location_x] = elements[last];
            lookup[elements[location_x]] = location_x;
            lookup[x] = last;
            elements[last] = x;
            last--;
        }

        /**
         * Move a node from the current set to another one.
         *
         * @param x Node to move.
         * @param other Other array set.
         */
        void move_to(int x, array_set &other) {
            remove(x);
            other.insert(x);
        }

        /**
         * Returns the size of the set.
         *
         * @return Size of the set.
         */
        unsigned int size() const { return last - first + 1; }

        /**
         * Check if the set is empty.
         *
         * @return True, if the set is empty.
         *         False, otherwise.
         */
        bool empty() const { return last < first; }

        /**
         * Returns an iterator to the beginning of the set.
         *
         * @return Iterator to the beginning of the set.
         */
        std::vector<int>::iterator begin() { return elements.begin() + first; }

        /**
         * Returns an iterator to the end of the set.
         *
         * @return Iterator to the end of the set.
         */
        std::vector<int>::iterator end()   { return elements.begin() + last + 1; }

        /**
         * Returns a _const_ iterator to the beginning of the set.
         *
         * @return Iterator to the beginning of the set.
         */
        std::vector<int>::const_iterator begin() const { return elements.begin() + first; }

        /**
         * Returns a _const_ iterator to the end of the set.
         *
         * @return Iterator to the end of the set.
         */
        std::vector<int>::const_iterator end() const { return elements.begin() + last + 1; }

        /**
         * Returns the element at a given index.
         *
         * @param index Index of the element.
         * @return Element at the index.
         */
        int at(unsigned int index) {
            return elements[index];
        }

        /**
         * Returns the element at a given index.
         *
         * @param index Index of the element.
         * @return Element at the index.
         */
        int operator[](unsigned int index) { return at(index); }

        /**
         * Clear the set.
         * Simply resets indiciators for the beginning and end.
         */
        void clear() {
            first = 0;
            last = -1;
        }

        /**
         * Checks if two sets are the same.
         * Two sets are equal if they contain the same elements.
         *
         * @return True, if the sets are the same.
         *         False, otherwise.
         */
        bool operator==(array_set const &that) const {
            if (size() != that.size()) return false;
            for (int const value : *this) {
                if (!that.contains(value)) return false;
            }
            return true;
        }

        /**
         * Checks if two sets are _not_ the same.
         * Two sets are equal if they contain the same elements.
         *
         * @return True, if the sets not are the same.
         *         False, otherwise.
         */
        bool operator!=(array_set const &that) const {
            return !(*this == that);
        }

    private:
        // Lookup array.
        std::vector<int> lookup;
        // Array containing the actual elements.
        std::vector<int> elements;
        // Index of the first element.
        int first;
        // Index of the last element.
        int last;
};

#endif 
