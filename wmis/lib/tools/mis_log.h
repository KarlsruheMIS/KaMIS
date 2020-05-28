/******************************************************************************
 * mis_log.h
 *
 *****************************************************************************/

#ifndef _MIS_LOG_H_
#define _MIS_LOG_H_

#include <sstream>

#include "timer.h"
#include "mis_config.h"
#include "data_structure/graph_access.h"

class mis_log {
    public:
        /**
         * Get the singleton logger instance.
         * 
         * @return Instance of the logger.
         */
        static mis_log *instance() {
            static mis_log inst;
            return &inst;
        };

        /**
         * Set the config.
         *
         * @param config Config for the evolutionary algorithm.
         */
        void set_config(MISConfig & config);

        /**
         * Set the graph.
         *
         * @param G Graph representation.
         */
        void set_graph(graph_access & G);

        /**
         * Write the log to a file.
         */
        void write_log();

        /**
         * Add a newline to the log.
         */
        void print_newline();

        /**
         * Print information about the graph.
         */
        void print_graph();

        /**
         * Print the current config.
         */
        void print_config();

        /**
         * Print information about the reduction step.
         * Includes number of extracted nodes and resulting kernel size
         *
         * @param mis_config Config for the logger.
         * @param extracted_nodes Number of removed nodes.
         * @param kernel_size Number of remaining nodes.
         */
        void print_reduction(MISConfig & mis_config, unsigned int extracted_nodes, unsigned int kernel_size);

        /**
         * Print banner for the initialization phase.
         */
        void print_init_title();

        /**
         * Print the final results.
         */
        void print_results();

        /**
         * Print a title.
         */
        void print_title();
        
        /**
         * Restart the timer for the total time including IO, etc.
         */
        void restart_total_timer();

        /**
         * Update the size of the best solution.
         *
         * @param mis_config Config for the logger.
         * @param size Candidate to replace the best solution size.
         */
        void set_best_size(MISConfig & mis_config, unsigned int size);

        /**
         * Reset the size of the best solution.
         */
        void reset_best_size();

    private:
        // General information
        timer total_timer;
        std::stringstream filebuffer_string;
        MISConfig log_config;

        // Graph informations
        std::string graph_name;
        unsigned int number_of_nodes;
        unsigned int number_of_edges;
        unsigned int arc_scans;
        double avg_degree;
        double density;

        // Reduction information
        unsigned int current_is_size;
        unsigned int optimum_size;

        // Results information
        double total_time_taken;
        double time_taken_best;

        /**
         * Default Constructor.
         */
        mis_log();

        /**
         * Default Destructor.
         */
        virtual ~mis_log();
};

#endif
