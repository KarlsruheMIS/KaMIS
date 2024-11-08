/**
 * parse_parameters.h
 * Purpose: Parse command line parameters.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _PARSE_PARAMETERS_H_
#define _PARSE_PARAMETERS_H_

#include <omp.h>

#include "configuration_struction.h"

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
                     ::mmwis::MISConfig & mis_config,
                     std::string & graph_filename) {
    const char *progname = argv[0];

    // Setup the argtable structs
    struct arg_lit *help                    = arg_lit0(NULL, "help", "Print help.");
    struct arg_int *user_seed               = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
    struct arg_str *filename                = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file.");
    struct arg_str *output                  = arg_str0(NULL, "output", NULL, "Path to store resulting independent set.");
    struct arg_str *kernel                  = arg_str0(NULL, "kernel", NULL, "Path to store resulting kernel.");
    struct arg_dbl *time_limit              = arg_dbl0(NULL, "time_limit", NULL, "Time limit in s. Default 1000s.");
    struct arg_lit *console_log             = arg_lit0(NULL, "console_log", "Stream the log into the console");
    struct arg_lit *disable_checks          = arg_lit0(NULL, "disable_checks", "Disable sortedness check during I/O.");
	struct arg_lit *random_freenodes        = arg_lit0(NULL, "random_freenodes", "Randomly picks free nodes to maximize to IS instead of sorting them by weight.");
	struct arg_lit *disable_reduction       = arg_lit0(NULL, "disable_reduction", "Don't perforn any reductions.");
    struct arg_str *weight_source           = arg_str0(NULL, "weight_source", NULL, "Choose how the weights are assigned. Can be either: file (default), hybrid, uniform, geometric.");
    struct arg_str *reduction_style         = arg_str0(NULL, "reduction_style", NULL, "Choose the type of reductions appropriate for the input graph. Can be either: normal/sparse (default), dense/osm.");
    struct arg_lit *cyclicStrong            = arg_lit0(NULL, "cyclicStrong", "Switch to cyclicStrong configuration (default is cyclicFast)");

    struct arg_int *set_limit               = arg_int0(NULL, "set_limit", NULL, "Choose maximum number of new vertices allowed to be created during struction application");
    struct arg_int *struction_degree        = arg_int0(NULL, "struction_degree", NULL, "Choose maximum degree of vertex to perform struction on it.");
    struct arg_str *struction_type          = arg_str0(NULL, "struction_type", NULL, "Choose used struction type. Can be either: original, modified, extended (default), extended_reduced or none");
    struct arg_str *key_type                = arg_str0(NULL, "key_type", NULL, "Choose used vertex selection strategy. Can be either: random, degree, increase (default), approximate_increase");
    struct arg_dbl *key_reinsert_factor     = arg_dbl0(NULL, "key_reinsert_factor", NULL, "Choose reinsert factor beta for approximate increase selection. default 2");
    struct arg_dbl *global_blow_up_percent  = arg_dbl0(NULL, "global_blow_up_factor", NULL, "Choose global blow up factor alpha.");
    struct arg_dbl *phase_blow_up_factor    = arg_dbl0(NULL, "phase_blow_up_factor", NULL, "Choose percentual phase blow up factor gamma.");
    struct arg_int *phase_blow_ups          = arg_int0(NULL, "phase_blow_ups", NULL, "Choose fixed blow ups per phase Z.");
    struct arg_int *max_unimproving_phases  = arg_int0(NULL, "max_unimproving_phases", NULL, "Choose maximum unimproving blow up phases X.");
    struct arg_str *backtrack_style         = arg_str0(NULL, "backtrack_style", NULL, "Choose backtrack strategy. Can be either: immediate_tie_breaking, immediate_exclude (default), end_min_kernel or no_backtrack");
    /* struct arg_lit *disable_neighbor_clique = arg_lit0(NULL, "disable_neighborhood_clique", "Disable neighborhood clique reduction."); */
    struct arg_lit *disable_generalized_fold= arg_lit0(NULL, "disable_generalized_fold", "Disable generalized fold reduction.");
    struct arg_lit *disable_critical_set    = arg_lit0(NULL, "disable_critical_set", "Disable critical set reduction.");
    struct arg_lit *disable_clique          = arg_lit0(NULL, "disable_clique", "Disable clique reduction.");
    struct arg_lit *disable_blow_up         = arg_lit0(NULL, "disable_blow_up", "Disable cyclic blow up algorithm.");
    struct arg_lit *plain_struction         = arg_lit0(NULL, "plain_struction", "Only use struction to reduce graph.");
    struct arg_lit *reduce_and_peel         = arg_lit0(NULL, "reduce_and_peel", "Use reduce-and-peel as initial solution for local search");
    struct arg_lit *ils         = arg_lit0(NULL, "ils", "Use ils as local search");

    struct arg_end *end                 = arg_end(100);

    // Setup the argtable
    void *argtable[] = {
            help,
            filename,
            output,
            kernel,
            user_seed,
            time_limit,
            console_log,
            cyclicStrong,
            disable_checks,
            random_freenodes,
            disable_reduction,
            weight_source,
            reduction_style,
            struction_degree,
            set_limit,
            struction_type,
            key_type,
            key_reinsert_factor,
            global_blow_up_percent,
            phase_blow_up_factor,
            phase_blow_ups,
            max_unimproving_phases,
            backtrack_style,
            reduce_and_peel,
            ils,
            plain_struction,
            disable_generalized_fold,
            /* disable_neighbor_clique, */
            disable_critical_set,
            disable_clique,
            disable_blow_up,
            end
    };

    // Choose standard configuration
    configuration_struction cfg;
    cfg.cyclicFast(mis_config);

    if (cyclicStrong->count > 0 ) {
        cfg.cyclicStrong(mis_config);
    }

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

    if (filename->count > 0) {
        graph_filename = filename->sval[0];
    }

    if (user_seed->count > 0) {
        mis_config.seed = user_seed->ival[0];
    }

    if (time_limit->count > 0) {
        mis_config.time_limit = time_limit->dval[0];
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

	if (random_freenodes->count > 0) {
		mis_config.sort_freenodes = false;
	}

	if (disable_reduction->count > 0) {
		mis_config.perform_reductions = false;
	}

	if (weight_source->count > 0) {
		mis_config.setWeightSource(weight_source->sval[0]);
	}

	if (reduction_style->count > 0) {
		mis_config.setReductionStyle(reduction_style->sval[0]);
	}

    if (output->count > 0) {
        mis_config.output_filename = output->sval[0];
        mis_config.write_graph = true;
    } else {
        mis_config.write_graph = false;
    }

    if (kernel->count > 0) {
        mis_config.kernel_filename = kernel->sval[0];
        mis_config.write_kernel= true;
    } else {
        mis_config.write_kernel= false;
    }



    if (struction_degree->count > 0) {
        mis_config.struction_degree = struction_degree->ival[0];
    }
    if (set_limit->count > 0) {
        mis_config.set_limit = set_limit->ival[0];
    }

    if (struction_type->count > 0) {
        mis_config.setStructionType(struction_type->sval[0]);
    }
    if (key_type->count > 0) {
        mis_config.setKeyType(key_type->sval[0]);
    }
    if (key_reinsert_factor->count > 0) {
        mis_config.key_reinsert_factor = key_reinsert_factor->dval[0];
    }
    if (global_blow_up_percent->count > 0) {
        mis_config.global_blow_up_factor = global_blow_up_percent->dval[0];
    }
    if (phase_blow_up_factor->count > 0) {
        mis_config.phase_blow_up_factor = phase_blow_up_factor->dval[0];
    }
    if (phase_blow_ups->count > 0) {
        mis_config.phase_blow_ups = phase_blow_ups->ival[0];
    }
    if (max_unimproving_phases->count > 0) {
        mis_config.max_unimproving_phases = max_unimproving_phases->ival[0];
    }
    if (backtrack_style->count > 0) {
        mis_config.setBacktrackType(backtrack_style->sval[0]);
    }
    if (reduce_and_peel->count > 0) {
        mis_config.reduce_and_peel = true;
    }
    if (ils->count > 0) {
        mis_config.perform_hils = false;
    }

    if (plain_struction->count) {
        mis_config.plain_struction = true;
    }
    if (disable_generalized_fold->count > 0) {
        mis_config.disable_generalized_fold = true;
    }
    /* if (disable_neighbor_clique->count > 0) { */
    /*     mis_config.neighborhood_clique_disabled = true; */
    /* } */
    if (disable_critical_set->count > 0) {
        mis_config.disable_critical_set = true;
    }
    if (disable_clique->count > 0) {
        mis_config.disable_clique = true;
    }
    if (disable_blow_up->count > 0) {
        mis_config.disable_blow_up = true;
    }
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

    return 0;
}

#endif
