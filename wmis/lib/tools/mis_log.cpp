/******************************************************************************
 * mis_log.cpp
 *****************************************************************************/

#include "mis_log.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

mis_log::mis_log() {
    number_of_nodes = 0;
    number_of_edges = 0;
    avg_degree = 0.0;
    density = 0.0;
    arc_scans = 0;
    total_time_taken = 0.0;
    time_taken_best = 0.0;
    optimum_size = 0;
}

mis_log::~mis_log() {
    
}

void mis_log::set_config(MISConfig & config) {
    log_config = config; 
}

void mis_log::set_graph(graph_access & G) {
    number_of_nodes = G.number_of_nodes();
    number_of_edges = G.number_of_edges();
    avg_degree = (double) number_of_edges / number_of_nodes;
    density = (double) (2 * number_of_edges) / (number_of_nodes * (number_of_nodes - 1));
}

void mis_log::write_log() {
    std::stringstream filename_stream;
    filename_stream << "./logs/log_"<<  log_config.graph_filename <<   
                       "_seed_" <<  log_config.seed;
    std::ofstream f(filename_stream.str());
    f << filebuffer_string.str();
    f.close();
}

void mis_log::print_newline() {
    filebuffer_string << std::endl; 

    if (log_config.console_log) {
        std::cout << std::endl; 
    }
}

void mis_log::print_title() {
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "\t\tWeightedMIS"                                                          << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;

    std::cout << "=========================================="                           << std::endl;
    std::cout << "\t\tWeightedMIS"                                                          << std::endl;
    std::cout << "=========================================="                           << std::endl;
}

void mis_log::print_graph() {
    filebuffer_string << "\t\tGraph"        << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "IO time:\t\t\t\t"         << total_timer.elapsed()                    << std::endl;
    filebuffer_string << "Filename:\t\t\t\t"        << log_config.graph_filename                << std::endl;
    filebuffer_string << "|-Nodes:\t\t\t\t"         << number_of_nodes                          << std::endl;
    filebuffer_string << "|-Edges:\t\t\t\t"         << number_of_edges                          << std::endl;
    filebuffer_string << std::endl;

    std::cout << "\t\tGraph"        << std::endl;
    std::cout << "=========================================="                           << std::endl;
    std::cout << "IO time:\t\t\t"           << total_timer.elapsed()                    << std::endl;
    std::cout << "Filename:\t\t\t"          << log_config.graph_filename                << std::endl;
    std::cout << "|-Nodes:\t\t\t"           << number_of_nodes                          << std::endl;
    std::cout << "|-Edges:\t\t\t"           << number_of_edges                          << std::endl;
    std::cout << std::endl;
}

void mis_log::print_config() {
    filebuffer_string << "\t\tConfiguration"        << std::endl;
    filebuffer_string << "=========================================="                            << std::endl;
    filebuffer_string << "Time limit:\t\t\t"         << log_config.time_limit                    << std::endl; 
    filebuffer_string << "Seed:\t\t\t\t"             << log_config.seed                          << std::endl; 
    filebuffer_string << "ILS iterations:\t\t\t"     << log_config.ils_iterations                << std::endl;
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << std::endl;
    
    std::cout << "\t\tConfiguration"        << std::endl;
    std::cout << "=========================================="                            << std::endl;
    std::cout << "Time limit:\t\t\t"         << log_config.time_limit                    << std::endl; 
    std::cout << "Seed:\t\t\t\t"             << log_config.seed                          << std::endl; 
    std::cout << "ILS iterations:\t\t\t"     << log_config.ils_iterations                << std::endl;
    std::cout << "---"                       << std::endl;
    std::cout << std::endl;
}

void mis_log::print_reduction(MISConfig & mis_config, unsigned int extracted_nodes, unsigned int kernel_size) {
    filebuffer_string << "\t\tReduction" << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "Offset:\t\t\t\t" << extracted_nodes << std::endl;
    filebuffer_string << "Kernel size:\t\t\t" << kernel_size << "\n" << std::endl;
    if (log_config.console_log) {
        std::cout << "=========================================="                           << std::endl;
        std::cout << "\t\tReduction" << std::endl;
        std::cout << "=========================================="                           << std::endl;
        std::cout << "Offset:\t\t\t\t" << extracted_nodes << std::endl;
        std::cout << "Kernel size:\t\t\t" << kernel_size << "\n" << std::endl;
    }
}

void mis_log::print_results() {
    filebuffer_string << std::endl;
    filebuffer_string << "\t\tStatistics"                                                          << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "Total time:\t\t\t\t"       << total_timer.elapsed()                   << std::endl;
    std::cout << std::endl;
    std::cout << "\t\tBest"                                                             << std::endl;
    std::cout << "=========================================="                           << std::endl;
    std::cout << "Size:\t\t\t\t"              << optimum_size                                          << std::endl;
    std::cout << "Time found:\t\t\t"          << time_taken_best                        << std::endl;
    std::cout << std::endl;
}

void mis_log::print_init_title() {
    filebuffer_string << "\t\tInitialization"        << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;

    if (log_config.console_log) {
        std::cout << "\t\tInitialization"        << std::endl;
        std::cout << "=========================================="                           << std::endl;
    }
}

void mis_log::restart_total_timer() {
    total_timer.restart();
}

void mis_log::reset_best_size() {
    optimum_size = 0;
}

void mis_log::set_best_size(MISConfig & mis_config, unsigned int size) {
    if (size > optimum_size) {
        optimum_size = size;
        time_taken_best = total_timer.elapsed();
    }
}

