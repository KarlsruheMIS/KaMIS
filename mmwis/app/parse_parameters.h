/**
 * parse_parameters.h
 * Purpose: Parse command line parameters.
 *
 *****************************************************************************/

#pragma once
#include <omp.h>

#include "configuration_mis.h"
#include "argtable3.h"
#include "string.h"

namespace mmwis {

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
    struct arg_str *user_conf           = arg_str0(NULL, "config", NULL, "Configuration to use. ([mmwis, mmwiss]). mmwiss is additionally using the struction algorithm in the initial population calculation.");
    struct arg_int *kahip_mode          = arg_int0(NULL, "kahip_mode", NULL, "preconfiguration KaHIP mode to use.");
    struct arg_dbl *imbalance           = arg_dbl0(NULL, "imbalance", NULL, "Set imbalance for partition and separator computation between. Default: 0.1");
    struct arg_str *filename            = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
    struct arg_str *output              = arg_str0(NULL, "output", NULL, "Path to store resulting independent set.");
    struct arg_str *kernel              = arg_str0(NULL, "kernel", NULL, "Path to store resulting kernel.");
    struct arg_dbl *time_limit          = arg_dbl0(NULL, "time_limit", NULL, "Total time limit in s. Default 1000s.");
    struct arg_dbl *evo_time_limit      = arg_dbl0(NULL, "evo_time_limit", NULL, "Time limit for one evolutionary round in s. Default 1000s.");
    struct arg_dbl *ils_time_limit      = arg_dbl0(NULL, "ils_time_limit", NULL, "Time limit for one ils round in s. Default 1000s.");
    struct arg_dbl *V_solve_exact       = arg_dbl0(NULL, "V_solve_exact", NULL, "Threshold of number of nodes, when the exact solver is started. Default V=0.");
    struct arg_dbl *time_solve_exact    = arg_dbl0(NULL, "time_solve_exact", NULL, "Time to take to solve kernel smaller than Threshold V_solve_exact.");
    struct arg_lit *console_log         = arg_lit0(NULL, "console_log", "Stream the log into the console");
    struct arg_lit *disable_checks      = arg_lit0(NULL, "disable_checks", "Disable sortedness check during I/O.");
    struct arg_str *weight_source       = arg_str0(NULL, "weight_source", NULL, "Choose how the weights are assigned. Can be either: file (default), hybrid, uniform, geometric, unit.");
    struct arg_str *vertex_selection    = arg_str0(NULL, "vertex_selection", NULL, "Choose vertex selection strategie to be applied. Can be either: degree, weight, hybrid, weight_over_degree, solution_participation (default).");
    struct arg_str *reduction_style     = arg_str0(NULL, "reduction_style", NULL, "Choose the type of reductions appropriate for the input graph. Can be either: initial (default), weight, time, time_and_weight.");
	struct arg_int *heavy_set           = arg_int0(NULL, "heavy_set",NULL, " Set size constraint of neighborhood in heavy_set reduction. Set to 0, then don't perform heavy set reduction.");
    struct arg_lit *disable_fold1       = arg_lit0(NULL, "disable_fold1", "Disable fold1 reduction.");
    struct arg_lit *disable_neighborhood       = arg_lit0(NULL, "disable_neighborhood", "Disable neighborhood reduction.");
    struct arg_lit *disable_twin       = arg_lit0(NULL, "disable_twin", "Disable twin reduction.");
    struct arg_lit *disable_clique       = arg_lit0(NULL, "disable_clique", "Disable clique reduction.");
    struct arg_lit *disable_triangle       = arg_lit0(NULL, "disable_triangle", "Disable triangle reduction.");
    struct arg_lit *disable_v_shape_min       = arg_lit0(NULL, "disable_v_shape_min", "Disable v_shape_min reduction.");
    struct arg_lit *disable_v_shape_mid       = arg_lit0(NULL, "disable_v_shape_mid", "Disable v_shape_max reduction.");
    struct arg_lit *disable_v_shape_max       = arg_lit0(NULL, "disable_v_shape_max", "Disable v_shape_max reduction.");
    struct arg_lit *disable_basic_se = arg_lit0(NULL, "disable_basic_se", "Disable basic_se reduction.");
    struct arg_lit *disable_extended_se = arg_lit0(NULL, "disable_extended_se", "Disable extended reduction.");
    struct arg_lit *disable_generalized_fold = arg_lit0(NULL, "disable_generalized_fold", "Disable generalized_fold reduction.");
    struct arg_lit *disable_heavy_set = arg_lit0(NULL, "disable_heavy_set", "Disable heavy_set reduction.");
    struct arg_lit *disable_critical_set = arg_lit0(NULL, "disable_critical_set", "Disable critical set reduction.");
    struct arg_lit *disable_clique_neighborhood = arg_lit0(NULL, "disable_clique_neighborhood", "Disable clique_neighborhood reduction.");
    struct arg_lit *disable_generalized_neighborhood = arg_lit0(NULL, "disable_generalized_neighborhood", "Disable generalized neighborhood reduction.");

    // for struction
    struct arg_int *use_struction_initial_sol = arg_int0(NULL, "use_struction_initial_sol",NULL, "use struction for given initial solution additionally.");
    struct arg_dbl *fraction                  = arg_dbl0(NULL, "fraction",NULL, "percentage fraction in heuristicReduce instead of a single node. Value has to be between 0 and 100. default:100");
    struct arg_end *end                     = arg_end(100);

    // Setup the argtable
    void *argtable[] = {
            help,
            filename,
            output,
            kernel,
            user_seed,
            user_conf, 
            kahip_mode,
            imbalance,
            // use_max_flow,
            // use_multiway_ns,
            // use_multiway_vc,
            // repetitions, 
            /* red_thres, */
            /* local_search_thres, */
            disable_basic_se,
            disable_extended_se,
            disable_generalized_fold,
            disable_fold1,
            disable_neighborhood,
            disable_twin,
            disable_clique,
            disable_triangle,
            disable_v_shape_min,
            disable_v_shape_mid,
            disable_v_shape_max,
            disable_heavy_set,
            disable_critical_set,
            disable_clique_neighborhood,
            disable_generalized_neighborhood,
            time_limit,
            evo_time_limit,
            ils_time_limit,
            V_solve_exact,
            time_solve_exact,
            console_log,
            disable_checks,
			// random_freenodes,
			// disable_reduction,
            weight_source,
            vertex_selection,
            reduction_style,
            use_struction_initial_sol,
            fraction,
            end
    };

    // Choose standard configuration
    configuration_mis cfg;
    cfg.mmwis(mis_config);

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
        if (strcmp(user_conf->sval[0], "mmwiss") == 0) cfg.mmwiss(mis_config);
    }

    if (filename->count > 0) {
        graph_filename = filename->sval[0];
    }

    if (imbalance->count > 0) {
        mis_config.imbalance = imbalance->dval[0];
    }

    if (kahip_mode->count > 0) {
        mis_config.kahip_mode = kahip_mode->ival[0];
    }

    if (user_seed->count > 0) {
        mis_config.seed = user_seed->ival[0];
    }

    if (disable_fold1->count > 0) {
        mis_config.disable_fold1 = true;
    }

    if (disable_neighborhood->count > 0) {
        mis_config.disable_neighborhood = true;
    }

    if (disable_twin->count > 0) {
        mis_config.disable_twin = true;
    }
    
    if (disable_clique->count > 0) {
        mis_config.disable_clique = true;
    }

    if (disable_triangle->count > 0) {
        mis_config.disable_triangle = true;
    }

    if (disable_v_shape_min->count > 0) {
        mis_config.disable_v_shape_min = true;
    }

    if (disable_v_shape_mid->count > 0) {
        mis_config.disable_v_shape_mid = true;
    }

    if (disable_v_shape_max->count > 0) {
        mis_config.disable_v_shape_max = true;
    }

    if (disable_basic_se->count > 0) {
        mis_config.disable_basic_se = true;
    }

    if (disable_extended_se->count > 0) {
        mis_config.disable_extended_se = true;
    }

    if (disable_generalized_fold->count > 0) {
        mis_config.disable_generalized_fold = true;
    }

    if (disable_heavy_set->count > 0) {
        mis_config.disable_heavy_set = true;
    }

    if (disable_critical_set->count > 0) {
        mis_config.disable_critical_set = true;
    }

    if (disable_clique_neighborhood->count > 0) {
        mis_config.disable_clique_neighborhood = true;
    }

    if (disable_generalized_neighborhood->count > 0) {
        mis_config.disable_generalized_neighborhood = true;
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

    /* if (red_thres->count > 0) { */
    /*     mis_config.reduction_threshold = red_thres->ival[0]; */
    /* } */

    /* if (local_search_thres->count > 0) { */
    /*     mis_config.local_search_threshold = local_search_thres->ival[0]; */
    /* } */

    if (time_limit->count > 0) {
        mis_config.time_limit = time_limit->dval[0];
    }

    if (ils_time_limit->count > 0) {
        mis_config.ils_time_limit = ils_time_limit->dval[0];
    }

    if (evo_time_limit->count > 0) {
        mis_config.evo_time_limit = evo_time_limit->dval[0];
    } else {
        mis_config.evo_time_limit = 0.1 * mis_config.time_limit;
    }

    if (V_solve_exact->count > 0) {
        mis_config.V_solve_exact= V_solve_exact->dval[0];
    } else {
        mis_config.V_solve_exact= 0;
    }

    if (time_solve_exact->count > 0) {
        mis_config.time_solve_exact= time_solve_exact->dval[0];
    } else {
        mis_config.time_solve_exact= mis_config.time_limit;
    }


    if (console_log->count > 0) {
        mis_config.console_log = true;
        mis_config.print_log = false;
    } else {
        mis_config.print_log = true;
    }

    if (disable_checks->count > 0) {
        mis_config.check_sorted = false;
    }

/* 	if (random_freenodes->count > 0) { */
  /* mis_config.sort_freenodes = false; */
/* 	} */

/* 	if (disable_reduction->count > 0) { */
  /* mis_config.perform_reductions = false; */
/* 	} */


    if (heavy_set->count > 0) {
        mis_config.heavy_set = heavy_set->ival[0];
    }

	if (weight_source->count > 0) {
		mis_config.setWeightSource(weight_source->sval[0]);
	}

	if (vertex_selection->count > 0) {
		mis_config.setVertexSelection(vertex_selection->sval[0]);
	}

	if (reduction_style->count > 0) {
		mis_config.setReductionStyle(reduction_style->sval[0]);
	} else {
		mis_config.setReductionStyle("initial");
    }


    if (output->count > 0) {
        mis_config.output_filename = output->sval[0];
        mis_config.write_graph = true;
    } else {
        mis_config.write_graph = false;
    }

    if (kernel->count > 0) {
        mis_config.kernel_filename = kernel->sval[0];
        mis_config.write_kernel = true;
    } else {
        mis_config.write_kernel = false;
    }

    if (fraction->count > 0) {
        mis_config.fraction = 0.01 * fraction->dval[0];
    } else {
        mis_config.fraction = 100;
    }

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}}

