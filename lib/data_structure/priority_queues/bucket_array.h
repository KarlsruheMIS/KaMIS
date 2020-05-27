/**
 * bucket_array.h
 * Purpose: A bucket priority queue used for building initial maximum independent sets.
 *
 * The original code from Andrade et al. was kindly provided by Renato Werneck.
 *
 *****************************************************************************/

#ifndef _BUCKET_ARRAY_H_
#define _BUCKET_ARRAY_H_

#include "definitions.h"

class bucket_array {
    public:
        /**
         * Constructor.
         * Initialize the bucket with nodes from 0..n-1.
         *
         * @param number_of_nodes Number of nodes.
         */
        bucket_array(unsigned int number_of_nodes);

        /**
         * Default Destructor.
         */
        virtual ~bucket_array();

        /**
         * Increment the value of a node by a given value.
         *
         * @param node Node whose value should be incremented.
         * @param value Value for incrementation.
         */
        void increment(NodeID node, int value);

        /**
         * Increment the value of a node by one.
         *
         * @param node Node whose value should be incremented.
         */
        void increment(NodeID node);

        /**
         * Decrement the value of a node by one.
         *
         * @param node Node whose value should be decremented.
         */
        void decrement(NodeID node);

        /**
         * Check if the bucket queue contains the given node.
         *
         * @param node Node to be looked for.
         * @return True if the queue contains the node.
         */
        bool contains(NodeID node);

        /**
         * Remove a node from the bucket queue.
         *
         * @param node Node to be removed.
         */
        void remove(NodeID node);

        /**
         * Remove a node from the bucket with the smallest value.
         * The node is chosen at random.
         *
         * @return Node from the lowest bucket.
         */
        NodeID pickSmallest();

    private:
        // Array containing the position of each node.
        unsigned int *position;
        // The array representing the buckets in increasing order.
        NodeID *array;
        // Array containing the value of each node.
        int *value;
        // Array containing the positition of each buckets start.
        unsigned int *first;
        // Total number of elements
        unsigned int size;
        // Removed elements. Starting point for non-removed elements.
        unsigned int non_negative;

        /**
         * Initializes the bucket with default values.
         */
        void init();

        /** 
         * Swap two nodes within the queue.
         *
         * @param first_pos Position of the first node.
         * @param second_pos Position of the second node.
         */
        void swap(int first_pos, int second_pos);

        /**
         * Check if a given bucket is empty:
         *
         * @param block Bucket to be checked.
         * @return True if the bucket is empty.
         */
        bool is_empty(int block);
};

#endif

