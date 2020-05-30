/**
 * online_mis.cpp
 * Purpose: Main program for the online reduction algorithm.
 *
 *****************************************************************************/

#include <argtable3.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "data_structure/graph_access.h"
#include "data_structure/mis_permutation.h"
#include "graph_io.h"
#include "greedy_mis.h"
#include "ils/online_ils.h"
#include "ils/local_search.h"
#include "mis_config.h"
#include "mis_log.h"
#include "parse_parameters_omis.h"
#include "timer.h"

int main(int argn, char **argv) {
  mis_log::instance()->restart_total_timer();
  mis_log::instance()->print_title();

  MISConfig mis_config;
  std::string graph_filepath;

  // Parse the command line parameters;
  int ret_code = parse_parameters(argn, argv, mis_config, graph_filepath);
  if (ret_code) {
    return 0;
  }
  mis_config.graph_filename =
      graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
  mis_log::instance()->set_config(mis_config);

  // Read the graph
  graph_access G;
  graph_io::readGraphWeighted(G, graph_filepath);
  mis_log::instance()->set_graph(G);

  // Graph copy for unfolding
  graph_access fG;
  graph_io::readGraphWeighted(fG, graph_filepath);

  // Print setup information
  mis_log::instance()->print_graph();
  //mis_log::instance()->print_config();

  // Perform the online reduction algorithm
  online_ils online;
  online.perform_ils(mis_config, G, mis_config.ils_iterations);

  // Gather best solution
  std::vector<NodeID> independent_set(G.number_of_nodes(), false);
  forall_nodes(G, node) {
    independent_set[node] = G.getPartitionIndex(node);
    fG.setPartitionIndex(node, G.getPartitionIndex(node));
  } endfor

  mis_log::instance()->print_results_online();
  //if (mis_config.print_log) mis_log::instance()->write_log();
  if (mis_config.write_graph)
    graph_io::writeIndependentSet(G, mis_config.output_filename);

  //std::cout << "insert folds ..." << std::endl;
  //std::stack<std::tuple<NodeID, NodeID, NodeID>> &folds = online.get_folds();
  //while (!folds.empty()) {
    //std::tuple<NodeID, NodeID, NodeID> fold = folds.top();
    //folds.pop();
    //if (fG.getPartitionIndex(std::get<0>(fold)) == 1) {
      //fG.setPartitionIndex(std::get<0>(fold), 0);
      //fG.setPartitionIndex(std::get<1>(fold), 1);
      //fG.setPartitionIndex(std::get<2>(fold), 1);
    //} else {
      //fG.setPartitionIndex(std::get<0>(fold), 1);
      //fG.setPartitionIndex(std::get<1>(fold), 0);
      //fG.setPartitionIndex(std::get<2>(fold), 0);
    //}
  //}

  //std::cout << "checking solution ..." << std::endl;
  //int counter = 0;
  //forall_nodes(fG, node) {
    //if (fG.getPartitionIndex(node)) {
      //counter++;
      //forall_out_edges(fG, e, node) {
        //NodeID target = fG.getEdgeTarget(e);
        //if (fG.getPartitionIndex(target)) {
          //std::cout << "not an independent set!" << std::endl;
          //exit(0);
        //}
      //} endfor
    //}
  //} endfor 
  //std::cout << "done " << counter << std::endl;

  return 0;
}
