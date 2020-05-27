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
int parse_parameters(int argn, char **argv,
                     MISConfig & mis_config,
                     std::string & graph_filename) {
    const char *progname = argv[0];

    // Setup the argtable structs
    struct arg_lit *help                = arg_lit0(NULL, "help", "Print help.");
    struct arg_int *user_seed           = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
    struct arg_str *user_conf           = arg_str0(NULL, "config", NULL, "Configuration to use. ([standard|social|full_standard|full_social]). Standard/social use different modes of the graph partitioning tool. Full configurations use more time consuming parameters.");
    
    struct arg_int *kahip_mode          = arg_int0(NULL, "kahip_mode", NULL, "KaHIP mode to use.");

    struct arg_str *kernelization_mode  = arg_str0(NULL, "kernelization", NULL, "Kernelization to use. ([FastKer|full]). Full is slower but produces smaller kernels (default: FastKer)");

    // struct arg_int *repetitions         = arg_int0(NULL, "repetitions", NULL, "Number of repetitions per round.");

    // struct arg_lit *use_multiway_ns     = arg_lit0(NULL, "use_multiway_ns", "Use the multiway combine operator with node separators.");
    // struct arg_lit *use_multiway_vc     = arg_lit0(NULL, "use_multiway_vc", "Use the multiway combine operator with vertex covers.");
    // struct arg_lit *use_hopcroft        = arg_lit0(NULL, "use_hopcroft", "Use Hopcroft-Karp to fix vertex cover candidates.");

    struct arg_str *filename            = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
    struct arg_str *output              = arg_str0(NULL, "output", NULL, "Path to store resulting independent set.");
    struct arg_dbl *time_limit          = arg_dbl0(NULL, "time_limit", NULL, "Time limit in s. Default 1000s.");
    struct arg_lit *console_log         = arg_lit0(NULL, "console_log", "Stream the log into the console");
    struct arg_lit *disable_checks      = arg_lit0(NULL, "disable_checks", "Disable sortedness check during I/O.");

    // Reduction
    // struct arg_dbl *best_degree_frac    = arg_dbl0(NULL, "best_degree_frac", NULL, "Fraction of degree nodes to remove before reduction.");
    struct arg_int *red_thres           = arg_int0(NULL, "red_thres", NULL, "Number of unsuccessful iterations before reduction.");

    struct arg_end *end                 = arg_end(100);

    // Setup the argtable
    void *argtable[] = {
            help, 
            filename, 
            output,
            user_seed, 
            user_conf, 
            kahip_mode,
	    kernelization_mode,
            // use_hopcroft,
            // use_multiway_ns,
            // use_multiway_vc,
            // repetitions, 
            red_thres,
            time_limit, 
            console_log,
            // best_degree_frac,
            disable_checks,
            end
    };    

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
        else if (strcmp(user_conf->sval[0], "social") == 0) cfg.social(mis_config);
        else if (strcmp(user_conf->sval[0], "full_standard") == 0) cfg.full_standard(mis_config);
        else if (strcmp(user_conf->sval[0], "full_social") == 0) cfg.full_social(mis_config);
    }

    mis_config.fullKernelization = false;
    if (kernelization_mode->count > 0) {
        if (strcmp(kernelization_mode->sval[0], "FastKer") == 0) mis_config.fullKernelization = false;
        else if (strcmp(kernelization_mode->sval[0], "full") == 0) mis_config.fullKernelization = true;
    }

    if (filename->count > 0) {
        graph_filename = filename->sval[0];
    }   

    if (kahip_mode->count > 0) {
        mis_config.kahip_mode = kahip_mode->ival[0];
    }

    if (user_seed->count > 0) {
        mis_config.seed = user_seed->ival[0];
    }

    // if (use_multiway_ns->count > 0) {
    //     mis_config.use_multiway_vc = false;
    // }

    // if (use_multiway_vc->count > 0) {
    //     mis_config.use_multiway_vc = true;
    // }

    // if (use_hopcroft->count > 0) {
    //     mis_config.use_hopcroft = true;
    // }

    // if (repetitions->count > 0) {
    //     mis_config.repetitions = repetitions->ival[0];
    // }

    if (red_thres->count > 0) {
        mis_config.reduction_threshold = red_thres->ival[0];
    }

    if (time_limit->count > 0) {
        mis_config.time_limit = time_limit->dval[0];
    }

    //if (best_degree_frac->count > 0) {
        //mis_config.remove_fraction = best_degree_frac->dval[0];
    //}

    if (console_log->count > 0) {
        mis_config.console_log = true;
        mis_config.print_log = false;
    } else {
        mis_config.print_log = true;
    }

    if (disable_checks->count > 0) {
        mis_config.check_sorted = false;
    }

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
