//
// Created by alex on 22.03.20.
//

#ifndef COMPONENTS_HILS_H
#define COMPONENTS_HILS_H

#include "definitions.h"
#include "graph_access.h"
#include "mmwis_config.h"
#include "random_functions.h"
#include "timer.h"

#include "Solution.h"

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <vector>

class hils {
private:
    timer t;
    const mmwis::MISConfig &config;

public:
    hils(const mmwis::MISConfig &config) : config(config) {}

	template <typename graph>
    void make_maximal(graph & G) {
        srand(config.seed);
        random_functions::setSeed(config.seed);

        Solution s(&G);
        while (!s.isMaximal()) {
            s.addRandomVertex();
        }
        NodeID i = 0;
        for (int state : s.solution()) {
            G.setPartitionIndex(i++, state);
        }
        assert(s.integrityCheck());
    }
    
	template <typename graph>
    void add_candidates(graph &G, std::vector<NodeID>& cand_list) {
        Solution s(&G);
        for (NodeID cand : cand_list) {
            s.force_candidate(cand);
        }
    }

	template <typename graph>
    void direct_improvement(graph &G, NodeID node) {
        double p[4] = {2, 4, 4, 1}; 
        Solution s(&G);

        while (!s.isMaximal()) {
                s.addRandomVertex();
        }

        bool improvement = s.candOmegaImprovement(node);
        if (!improvement) improvement = s.candTwoImprovement(node);

        NodeID i = 0;
        for (int state : s.solution()) {
            G.setPartitionIndex(i++, state);
        }
        assert(s.integrityCheck());
    }

	template <typename graph>
    void direct_improvement(graph&G, std::vector<NodeID>& candidates) {
        double p[4] = {2, 4, 4, 1}; // intensification/exploration parameters
        Solution s(&G);

        make_maximal(G);
        bool improvement = false;
        for (NodeID c : candidates) {
            while (!s.isMaximal()) {
                s.addRandomVertex();
            }
            improvement = s.candOmegaImprovement(c);
            if (!improvement) improvement = s.candTwoImprovement(c);
        }

        NodeID node = 0;
        for (int state : s.solution()) {
            G.setPartitionIndex(node++, state);
        }
        assert(s.integrityCheck());
    }


	template <typename graph>
    void perform_ils(graph &G, unsigned int max_iterations, double time_limit,  NodeWeight weight_offset = 0) {
        double p[4] = {2, 4, 4, 1}; // intensification/exploration parameters
        Solution s(&G);


        // Init RNG
        srand(config.seed);
        random_functions::setSeed(config.seed);

        forall_nodes(G, n)
            if (G.getPartitionIndex(n))
                s.addVertex(n);
        endfor
        while (!s.isMaximal()) {
            s.addRandomVertex();
        }
        // found initial solution
        do {
            while (!s.isMaximal()) {
                s.addRandomVertex();
            }
            if (t.elapsed() > time_limit) break;
        } while (s.omegaImprovement() || s.twoImprovement() );

        Solution best_s(s);
        double best_t = 0;

        // run ILS iterations

        int k = 1;
        int local_best = s.weight();
        int iter;
        for (iter = 0; t.elapsed() < time_limit && iter < max_iterations; iter++) {
            Solution next_s(s);

            // shake
            next_s.force(p[0]);

            assert(next_s.integrityCheck());

            do {
                while (!next_s.isMaximal()) {
                    next_s.addRandomVertex();
                }
            } while (next_s.omegaImprovement() || next_s.twoImprovement());

            assert(best_s.integrityCheck());

            if (next_s.weight() > s.weight()) {
                k = 1;
                s = next_s;

                if (local_best < next_s.weight()) {
                    k -= s.size() / p[1];
                    local_best = next_s.weight();
                }

                if (best_s.weight() < s.weight()) {
                    best_s = s;
                    k -= s.size() * p[2];
                    best_t = t.elapsed();
                }
            } else if (k <= s.size() / p[1]) {
                k++;
            } else {
                local_best = s.weight();
                s.force(p[3]);
                k = 1;
            }
        }

        NodeID i = 0;
        for (int state : best_s.solution()) {
            G.setPartitionIndex(i++, state);
        }
        assert(best_s.integrityCheck());
        std::cout << "hils weight: " << best_s.weight()<< std::endl;
        std::cout << "hils time: " << best_t<< std::endl;
    }
};

#endif //COMPONENTS_HILS_H
