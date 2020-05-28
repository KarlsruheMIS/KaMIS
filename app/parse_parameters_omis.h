/**
 * parse_parameters.h
 * Purpose: Parse command line parameters.
 *
 *****************************************************************************/

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <omp.h>

#include "configuration_mis.h"

/**
 * Parse the given parameters and apply them to the config.
 *
 * @param argn Number of parameters.
 * @param argv Values of the parameters.
 * @param mis_config Config to store the values in.
 * @param graph_filename String to store the filename of the graph.
 *
 * @return -1 if there was an error. 0 otherwise.
 */
int parse_parameters(int argn, char **argv, MISConfig &mis_config,
                     std::string &graph_filename) {
  const char *progname = argv[0];

  // Setup the argtable structs
  struct arg_lit *help = arg_lit0(NULL, "help", "Print help.");
  struct arg_int *user_seed =
      arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
  struct arg_str *user_conf =
      arg_str0(NULL, "config", NULL, "Configuration to use. ([standard]).");

  struct arg_str *filename =
      arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
  struct arg_str *output = arg_str0(NULL, "output", NULL,
                                    "Path to store resulting independent set.");
  struct arg_dbl *time_limit =
      arg_dbl0(NULL, "time_limit", NULL, "Time limit in s. Default 1000s.");
  struct arg_lit *console_log =
      arg_lit0(NULL, "console_log", "Stream the log into the console");
  struct arg_lit *disable_checks =
      arg_lit0(NULL, "disable_checks", "Disable sortedness check during I/O.");
  struct arg_lit *adaptive_greedy =
      arg_lit0(NULL, "adaptive_greedy", "Use adaptive greedy solution");

  struct arg_end *end = arg_end(100);

  // Setup the argtable
  void *argtable[] = {help,      filename,        output,      user_seed,
                      time_limit,      console_log, disable_checks,
                      adaptive_greedy, end};

  // Choose standard configuration
  configuration_mis cfg;
  cfg.standard(mis_config);

  // Parse the arguments
  int nerrors = arg_parse(argn, argv, argtable);

  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, "  %-40s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 1;
  }

  if (nerrors > 0) {
    arg_print_errors(stderr, end, progname);
    printf("Try '%s --help' for more information.\n", progname);
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 1;
  }

  if (user_conf->count > 0) {
    if (strcmp(user_conf->sval[0], "standard") == 0) cfg.standard(mis_config);
  }

  if (filename->count > 0) {
    graph_filename = filename->sval[0];
  }

  if (user_seed->count > 0) {
    mis_config.seed = user_seed->ival[0];
  }

  if (time_limit->count > 0) {
    mis_config.time_limit = time_limit->dval[0];
  }

  mis_config.console_log = console_log->count > 0;
  mis_config.print_log = !mis_config.console_log;

  mis_config.check_sorted = disable_checks->count > 0;

  mis_config.start_greedy_adaptive = adaptive_greedy->count > 0;

  if (output->count > 0) {
    mis_config.output_filename = output->sval[0];
    mis_config.write_graph = true;
  } else {
    mis_config.write_graph = false;
  }

  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return 0;
}

#endif
