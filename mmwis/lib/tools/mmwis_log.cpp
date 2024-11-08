/******************************************************************************
 * mmwis_log.cpp
 *****************************************************************************/

#include "mmwis_log.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace mmwis;

mmwis_log::mmwis_log() {
    number_of_nodes = 0;
    number_of_edges = 0;
    avg_degree = 0.0;
    density = 0.0;
    arc_scans = 0;
    number_of_rounds = 0;
    number_of_repetitions = 0;
    best_solution_size = 0;
    best_solution_weight = 0;
    result_operator = 0;
    total_time_taken = 0.0;
    separator_selected = 0;
    cover_selected = 0;
    multiway_selected = 0;
    separator_improv = 0;
    cover_improv = 0;
    multiway_improv = 0;
    total_separator_time = 0.0;
    total_cover_time = 0.0;
    total_multiway_time = 0.0;
    avg_separator_time = 0.0;
    avg_cover_time = 0.0;
    avg_multiway_time = 0.0;
    total_separator_improv = 0;
    total_cover_improv = 0;
    total_multiway_improv = 0;
    avg_separator_improv = 0.0;
    avg_cover_improv = 0.0;
    avg_multiway_improv = 0.0;
    time_taken_best = 0.0;
    repetition_best = 0;
    round_best = 0;
    time_for_building_pool = 0.0;
    time_since_building_pool = 0.0;
    current_is_size = 0;
    total_time_taken = 0.0;
    time_taken_best = 0.0;
    optimum_size = 0;
    optimum_weight = 0;
}

mmwis_log::~mmwis_log() {
    
}


void mmwis_log::print_graph() {
    filebuffer_string << "\t\tGraph"        << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "IO time:\t\t\t\t"         << overall_total_timer.elapsed()            << std::endl;
    filebuffer_string << "Filename:\t\t\t\t"        << log_config.graph_filename                << std::endl;
    filebuffer_string << "|-Nodes:\t\t\t\t"         << number_of_nodes                          << std::endl;
    filebuffer_string << "|-Edges:\t\t\t\t"         << number_of_edges                          << std::endl;
    filebuffer_string << std::endl;

    std::cout << "\t\tGraph"        << std::endl;
    std::cout << "==========================================="                          << std::endl;
    std::cout << "IO time:\t\t\t"           << overall_total_timer.elapsed()            << std::endl;
    std::cout << "Filename:\t\t\t"          << log_config.graph_filename                << std::endl;
    std::cout << "|-Nodes:\t\t\t"           << number_of_nodes                          << std::endl;
    std::cout << "|-Edges:\t\t\t"           << number_of_edges                          << std::endl;
    std::cout << std::endl;
}

void mmwis_log::print_config() {
    filebuffer_string << "\t\tConfiguration"        << std::endl;
    filebuffer_string << "=========================================="                            << std::endl;
    filebuffer_string << "Configuration:\t\t"        << log_config.config_name                   << std::endl; 
    filebuffer_string << "Vertex Selection:\t\t"     << log_config.vertex_selection              << std::endl; 
    filebuffer_string << "threshold to solve exact:\t"<< log_config.V_solve_exact                << std::endl; 
    filebuffer_string << "time to solve exact:\t"    << log_config.time_solve_exact              << std::endl; 
    filebuffer_string << "fraction:\t\t"             << log_config.fraction                      << std::endl; 
    filebuffer_string << "Population size:\t\t"      << log_config.population_size               << std::endl; 
    filebuffer_string << "Repetitions:\t\t\t"        << log_config.repetitions                   << std::endl; 
    filebuffer_string << "Time limit:\t\t\t"         << log_config.time_limit                    << std::endl; 
    filebuffer_string << "Evo time limit:\t\t"       << log_config.evo_time_limit                    << std::endl; 
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "KaHIP mode:\t\t\t"         << log_config.kahip_mode                    << std::endl;
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "Seed:\t\t\t\t"             << log_config.seed                          << std::endl; 
    filebuffer_string << "Base imbalance:\t\t\t"     << log_config.imbalance                     << std::endl;
    filebuffer_string << "Random imbalance:\t\t"     << log_config.randomize_imbalance           << std::endl;
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "Tournament size:\t\t"      << log_config.tournament_size               << std::endl; 
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "Mutation rate:\t\t\t"      << log_config.flip_coin                     << std::endl; 
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "Use Multiway-VC:\t\t"      << log_config.use_multiway_vc               << std::endl; 
    filebuffer_string << "Multiway blocks:\t\t"      << log_config.multiway_blocks               << std::endl; 
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "Apply all reductions:\t\t" << log_config.all_reductions                << std::endl; 
    filebuffer_string << "Reductions threshold:\t\t" << log_config.reduction_threshold           << std::endl;
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "Insert threshold:\t\t"     << log_config.insert_threshold              << std::endl;
    filebuffer_string << "Pool threshold:\t\t\t"     << log_config.pool_threshold                << std::endl;
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "ILS iterations:\t\t\t"     << log_config.ils_iterations                << std::endl;
    filebuffer_string << "---"                       << std::endl;
    filebuffer_string << "No. of partitions:\t\t"    << log_config.number_of_partitions          << std::endl;
    filebuffer_string << "No. of separators:\t\t"    << log_config.number_of_separators          << std::endl; 
    filebuffer_string << "No. of k-partitions:\t\t"  << log_config.number_of_k_partitions        << std::endl;
    filebuffer_string << "No. of k-separators:\t\t"  << log_config.number_of_k_separators        << std::endl; 
    filebuffer_string << std::endl;
    
    std::cout << "\t\tConfiguration"        << std::endl;
    std::cout << "==========================================="                           << std::endl;
    std::cout << "Configuration:\t\t\t"        << log_config.config_name                   << std::endl; 
    std::cout << "Vertex Selection:\t\t"     << log_config.vertex_selection              << std::endl; 
    std::cout << "Threshold to solve exact:\t"<< log_config.V_solve_exact                << std::endl; 
    std::cout << "Time to solve exact:\t\t"    << log_config.time_solve_exact              << std::endl; 
    std::cout << "Fraction:\t\t\t"           << log_config.fraction                      << std::endl; 
    std::cout << "Population size:\t\t"      << log_config.population_size               << std::endl; 
    std::cout << "Repetitions:\t\t\t"        << log_config.repetitions                   << std::endl; 
    std::cout << "Time limit:\t\t\t"         << log_config.time_limit                    << std::endl; 
    std::cout << "Evo time limit:\t\t\t"     << log_config.evo_time_limit                << std::endl; 
    std::cout << "---"                       << std::endl;
    std::cout << "KaHIP mode:\t\t\t"         << log_config.kahip_mode                    << std::endl;
    std::cout << "---"                       << std::endl;
    std::cout << "Seed:\t\t\t\t"             << log_config.seed                          << std::endl; 
    std::cout << "Base imbalance:\t\t\t"     << log_config.imbalance                     << std::endl;
    std::cout << "Random imbalance:\t\t"     << log_config.randomize_imbalance           << std::endl;
    std::cout << "---"                       << std::endl;
    std::cout << "Tournament size:\t\t"      << log_config.tournament_size               << std::endl; 
    std::cout << "---"                       << std::endl;
    std::cout << "Mutation rate:\t\t\t"      << log_config.flip_coin                     << std::endl; 
    std::cout << "---"                       << std::endl;
    std::cout << "Use Multiway-VC:\t\t"      << log_config.use_multiway_vc               << std::endl; 
    std::cout << "Multiway blocks:\t\t"      << log_config.multiway_blocks               << std::endl; 
    std::cout << "---"                       << std::endl;
    std::cout << "Apply all reductions:\t\t" << log_config.all_reductions                << std::endl; 
    std::cout << "Reductions threshold:\t\t" << log_config.reduction_threshold           << std::endl;
    std::cout << "---"                       << std::endl;
    std::cout << "Insert threshold:\t\t"     << log_config.insert_threshold              << std::endl;
    std::cout << "Pool threshold:\t\t\t"     << log_config.pool_threshold                << std::endl;
    std::cout << "---"                       << std::endl;
    std::cout << "ILS iterations:\t\t\t"     << log_config.ils_iterations                << std::endl;
    std::cout << "---"                       << std::endl;
    std::cout << "No. of partitions:\t\t"    << log_config.number_of_partitions          << std::endl;
    std::cout << "No. of separators:\t\t"    << log_config.number_of_separators          << std::endl; 
    std::cout << "No. of k-partitions:\t\t"  << log_config.number_of_k_partitions        << std::endl;
    std::cout << "No. of k-separators:\t\t"  << log_config.number_of_k_separators        << std::endl; 

    std::cout << std::endl;
}



void mmwis_log::print_repetition(mmwis::MISConfig & mis_config) {
    std::string skip = (evo_operator == "Initial")? "\t\t" : "\t";
    std::string skip_operator;
    if (evo_operator == "Multiway" || 
        evo_operator == "Initial" || 
        evo_operator == "Collect") skip_operator = "\t\t\t";
    else if (evo_operator == "Vertex cover" || 
             evo_operator == "Node separator" || 
             evo_operator == "Global collect") skip_operator = "\t\t";
    else skip_operator = "\t";
    if (mis_config.print_repetition) {
        filebuffer_string               << number_of_repetitions << "\t\t"
                                        << evo_operator << skip
                                        << result_operator << "\t"
                                        << best_solution_weight << "\t"
                                        << total_timer.elapsed()                          << std::endl;
    }

    if (mis_config.print_repetition) {
        std::cout               << number_of_repetitions << "\t\t"
                                << evo_operator << skip
                                << result_operator << "\t"
                                << best_solution_weight << "\t"
                                << total_timer.elapsed()                                  << std::endl;
    }
}

void mmwis_log::print_results() {
    compute_avg();
    filebuffer_string << std::endl;
    filebuffer_string << "\t\tStatistics"                                                                               << std::endl;
    filebuffer_string << "==========================================="                                                  << std::endl;
    filebuffer_string << "Overall total time: (including IO etc)\t\t\t\t"    << overall_total_timer.elapsed()           << std::endl;
    filebuffer_string << "Total time:\t\t\t\t"                               << total_timer.elapsed()                   << std::endl;
    filebuffer_string << "Evolution time:\t\t\t"                             << evo_timer.elapsed()                     << std::endl;
    filebuffer_string << "Number of rounds:\t\t"                             << number_of_rounds                        << std::endl;
    filebuffer_string << "Number of repetitions:\t"                          << number_of_repetitions                   << std::endl;
    filebuffer_string << std::endl;
    filebuffer_string << "Operator\t\tTimes selected\tImprovements\tAvg. time\tAvg. improv"                             << std::endl;
    filebuffer_string << "--------------------------------------------------------------------------------"
              << std::endl;
    filebuffer_string                   << "Vertex cover\t" 
                                        << cover_selected << "\t\t\t\t" 
                                        << cover_improv << "\t\t\t\t" 
                                        << avg_cover_time << "\t\t"                                  
                                        << avg_cover_improv                                     << std::endl;
    filebuffer_string                   << "Node separator\t" 
                                        << separator_selected << "\t\t\t\t" 
                                        << separator_improv << "\t\t\t\t" 
                                        << avg_separator_time << "\t\t"
                                        << avg_separator_improv                                 << std::endl;
    filebuffer_string                   << "Multiway\t\t" 
                                        << multiway_selected << "\t\t\t\t" 
                                        << multiway_improv << "\t\t\t\t" 
                                        << avg_multiway_time << "\t\t"                            
                                        << avg_multiway_improv                                  << std::endl;
    filebuffer_string << std::endl;
    filebuffer_string << "\t\tBest"                                                             << std::endl;
    filebuffer_string << "==========================================="                          << std::endl;
    filebuffer_string << "Weight:\t\t\t\t\t"          << optimum_weight                         << std::endl;
    filebuffer_string << "final:\t\t\t\t\t"           << optimum_weight                         << std::endl;
    filebuffer_string << "Time found:\t\t\t\t"        << time_taken_best                        << std::endl;
    filebuffer_string << "Operator:\t\t\t\t"          << operator_best                          << std::endl;
    filebuffer_string << "Repetition number:\t\t"     << repetition_best                        << std::endl;
    filebuffer_string << "Round number:\t\t\t"        << round_best                             << std::endl;
    filebuffer_string << std::endl;

    std::cout << std::endl;
    std::cout << "\t\tStatistics"                                                                                       << std::endl;
    std::cout << "==========================================="                                                          << std::endl;
    filebuffer_string << "Overall total time: (including IO etc)\t\t\t\t"    << overall_total_timer.elapsed()           << std::endl;
    std::cout << "Total time:\t\t\t"                                         << total_timer.elapsed()                   << std::endl;
    std::cout << "Evolution time:\t\t\t"                                     << evo_timer.elapsed()                     << std::endl;
    std::cout << "Number of rounds:\t\t"                                     << number_of_rounds                        << std::endl;
    std::cout << "Number of repetitions:\t\t"                                << number_of_repetitions                   << std::endl;
    std::cout << std::endl;
    std::cout << "Operator\tTimes selected\tImprovements\tAvg. time\tAvg. improv"                                       << std::endl;
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
    std::cout                   << "Vertex cover\t" 
                                << cover_selected << "\t\t" 
                                << cover_improv << "\t\t" 
                                << avg_cover_time << "\t\t"                                  
                                << avg_cover_improv                                     << std::endl;
    std::cout                   << "Node separator\t" 
                                << separator_selected << "\t\t" 
                                << separator_improv << "\t\t" 
                                << avg_separator_time << "\t\t"
                                << avg_separator_improv                                 << std::endl;
    std::cout                   << "Multiway\t" 
                                << multiway_selected << "\t\t" 
                                << multiway_improv << "\t\t" 
                                << avg_multiway_time << "\t\t"                            
                                << avg_multiway_improv                                  << std::endl;
    std::cout << std::endl;
    std::cout << "\t\tBest"                                                             << std::endl;
    std::cout << "==========================================="                          << std::endl;
    std::cout << "Weight:\t\t\t\t"            << optimum_weight                         << std::endl;
    std::cout << "Time found:\t\t\t"          << time_taken_best                        << std::endl;
    std::cout << "Operator:\t\t\t"            << operator_best                          << std::endl;
    std::cout << "Repetition number:\t\t"     << repetition_best                        << std::endl;
    std::cout << "Round number:\t\t\t"        << round_best                             << std::endl;
    std::cout << std::endl;
}

void mmwis_log::print_pool_title() {
    filebuffer_string << "\t\tPartitioning"                                             << std::endl;
    filebuffer_string << "==========================================="                  << std::endl;

    if (log_config.console_log) {
        std::cout << "\t\tPartitioning"                                                 << std::endl;
        std::cout << "==========================================="                      << std::endl;
    }
}

void mmwis_log::print_evolution_title() {
    filebuffer_string << "\t\tEvolution"                                                << std::endl;
    filebuffer_string << "==========================================="                  << std::endl;

    if (log_config.console_log) {
        std::cout << "\t\tEvolution"                                                    << std::endl; 
        std::cout << "==========================================="                      << std::endl;
    }
}

void mmwis_log::print_init_title() {
    filebuffer_string << "\t\tInitialization"                                           << std::endl;
    filebuffer_string << "==========================================="                  << std::endl;

    if (log_config.console_log) {
        std::cout << "\t\tInitialization"                                               << std::endl;
        std::cout << "==========================================="                      << std::endl;
    }
}

void mmwis_log::print_separator() {
    filebuffer_string << "Created separator pool:\t" << pool_timer.elapsed()            << std::endl;
    filebuffer_string << std::endl;
    if (log_config.console_log) {
        std::cout << "Created separator pool:\t" << pool_timer.elapsed()                << std::endl;
        std::cout << std::endl;
    }
}

void mmwis_log::restart_overall_total_timer() {
    overall_total_timer.restart();
}

void mmwis_log::restart_total_timer() {
    total_timer.restart();
}

double mmwis_log::get_total_timer() {
    return total_timer.elapsed();
}


void mmwis_log::restart_evo_timer() {
    evo_timer.restart();
}

double mmwis_log::get_evo_timer() {
    return evo_timer.elapsed();
}

void mmwis_log::restart_operator_timer() {
    operator_timer.restart();
}

void mmwis_log::restart_building_pool_timer() {
    pool_timer.restart();
}


double mmwis_log::get_pool_building_time() {
    double elapsed = pool_timer.elapsed();
    non_pool_timer.restart();
    return elapsed;
}

double mmwis_log::get_after_pool_time() {
    return non_pool_timer.elapsed();
}

void mmwis_log::inc_rounds() {
    number_of_rounds++;
}

void mmwis_log::inc_repetitions() {
    number_of_repetitions++; 
}

void mmwis_log::set_operator(std::string operator_name) {
    evo_operator = operator_name;
    if (evo_operator == "Node separator") {
        separator_selected++;
    } else if (evo_operator == "Vertex cover") {
        cover_selected++;
    } else if (evo_operator == "Multiway") {
        multiway_selected++;
    }
}

void mmwis_log::set_result_operator(unsigned int result) {
    result_operator = result;    
    if (evo_operator == "Node separator") {
        total_separator_time += operator_timer.elapsed();
        if (result_operator > best_solution_weight) {
            total_separator_improv += result_operator - best_solution_weight;
            separator_improv++;
        }
    } else if (evo_operator == "Vertex cover") {
        total_cover_time += operator_timer.elapsed();
        if (result_operator > best_solution_weight) {
            total_cover_improv += result_operator - best_solution_weight;
            cover_improv++;
        }
    } else if (evo_operator == "Multiway") {
        total_multiway_time += operator_timer.elapsed();
        if (result_operator > best_solution_weight) {
            total_multiway_improv += result_operator - best_solution_weight;
            multiway_improv++;
        }
    }
}

void mmwis_log::compute_avg() {
    avg_cover_time = (cover_selected > 0)? total_cover_time / cover_selected : 0.0;
    avg_separator_time = (separator_selected > 0)? total_separator_time / separator_selected : 0.0;
    avg_multiway_time = (multiway_selected > 0)? total_multiway_time / multiway_selected : 0.0;

    avg_cover_improv = (cover_improv > 0)? total_cover_improv / cover_improv : 0.0;
    avg_separator_improv = (separator_improv > 0)? total_separator_improv / separator_improv : 0.0;
    avg_multiway_improv = (multiway_improv > 0)? total_multiway_improv / multiway_improv : 0.0;
}

void mmwis_log::set_avg_solution_size(double avg_size) {
    avg_solution_size = avg_size;
}

void mmwis_log::reset_best_size() {
    optimum_size = 0;
}

void mmwis_log::set_best_weight(mmwis::MISConfig & mis_config, NodeWeight  weight) {
    best_solution_weight = weight;
    if (best_solution_weight > optimum_weight) {
        print_repetition(mis_config);
        optimum_weight = best_solution_weight;
        time_taken_best = total_timer.elapsed();
        operator_best = evo_operator;
        repetition_best = number_of_repetitions;
        round_best = number_of_rounds;
    }
}

void mmwis_log::set_best_size(mmwis::MISConfig & mis_config, unsigned int size) {
    best_solution_size = size;
    if (best_solution_size > optimum_size) {
        print_repetition(mis_config);
        optimum_size = best_solution_size;
        time_taken_best = total_timer.elapsed();
        operator_best = evo_operator;
        repetition_best = number_of_repetitions;
        round_best = number_of_rounds;
    }
}

void mmwis_log::set_config(mmwis::MISConfig & config) {
    log_config = config; 
}

void mmwis_log::write_log() {
    std::stringstream filename_stream;
    filename_stream << "./logs/log_"<<  log_config.graph_filename <<   
                       "_seed_" <<  log_config.seed;
    std::ofstream f(filename_stream.str());
    f << filebuffer_string.str();
    f.close();
}

void mmwis_log::print_newline() {
    filebuffer_string << std::endl; 

    if (log_config.console_log) {
        std::cout << std::endl; 
    }
}

void mmwis_log::print_memetic_title() {
    filebuffer_string << "==========================================="                   << std::endl;
    filebuffer_string << "\t\tMemeticMWIS"                                               << std::endl;
    filebuffer_string << "==========================================="                   << std::endl;

    std::cout << "==========================================="                           << std::endl;
    std::cout << "\t\tMemeticMWIS"                                                       << std::endl;
    std::cout << "==========================================="                           << std::endl;
}

void mmwis_log::print_title() {
    filebuffer_string << "==========================================="                   << std::endl;
    filebuffer_string << "\t\tWeightedMIS"                                               << std::endl;
    filebuffer_string << "==========================================="                   << std::endl;

    std::cout << "==========================================="                           << std::endl;
    std::cout << "\t\tWeightedMIS"                                                       << std::endl;
    std::cout << "==========================================="                           << std::endl;
}



void mmwis_log::print_reduction(mmwis::MISConfig & mis_config, unsigned int extracted_nodes, unsigned int kernel_size) {
    filebuffer_string << "\t\tReduction"                                                         << std::endl;
    filebuffer_string << "==========================================="                           << std::endl;
    filebuffer_string << "Offset:\t\t\t\t"            << extracted_nodes                         << std::endl;
    filebuffer_string << "Kernel size:\t\t\t"         << kernel_size                             << std::endl;
    filebuffer_string << "Time:\t\t\t\t"              << total_timer.elapsed()           << "\n" << std::endl;
        std::cout << "==========================================="                               << std::endl;
        std::cout << "\t\tReduction"                                                             << std::endl;
        std::cout << "==========================================="                               << std::endl;
        std::cout << "Offset:\t\t\t\t"                << extracted_nodes                         << std::endl;
        std::cout << "Kernel size:\t\t\t"             << kernel_size                             << std::endl;
        std::cout << "Time:\t\t\t\t"                  << total_timer.elapsed()           << "\n" << std::endl;
}



void mmwis_log::print_round(mmwis::MISConfig & mis_config) {
    filebuffer_string << std::endl;
    filebuffer_string << "Round:\t\t\t\t"            << number_of_rounds                         << std::endl;
    filebuffer_string << "Best solution:\t\t\t"      << best_solution_weight                     << std::endl; 
    filebuffer_string << "Evo Time:\t\t\t"           << evo_timer.elapsed()                      << std::endl;
    filebuffer_string << "Total Time:\t\t\t\t"       << total_timer.elapsed()                    << std::endl;
    filebuffer_string << std::endl;

    std::cout << std::endl;
    std::cout << "Round:\t\t\t\t"                    << number_of_rounds                         << std::endl;
    std::cout << "Best solution:\t\t\t"              << best_solution_weight                     << std::endl; 
    std::cout << "Evo Time:\t\t\t"                   << evo_timer.elapsed()                      << std::endl;
    std::cout << "Time:\t\t\t\t"                     << total_timer.elapsed()                    << std::endl;
    std::cout << std::endl;
}
