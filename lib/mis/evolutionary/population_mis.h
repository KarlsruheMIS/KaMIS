/**
 * population_mis.h
 * Purpose: Represents the population used in the evolutionary framework.
 *          This includes combine, mutate and selection operators.
 *
 *****************************************************************************/

#ifndef _POPULATION_MIS_H_
#define _POPULATION_MIS_H_

#include "definitions.h"
#include "mis_config.h"
#include "kaHIP_interface.h"
#include "data_structure/graph_access.h"

/**
 * Representation of a single individuum.
 * An individuum contains a solution, its corresponding size and an id.
 */
struct individuum_mis {
    NodeID *solution;
    unsigned int solution_size;
    unsigned int id;
};

class population_mis {
    public:
        /**
         * Default Constructor.
         */
        population_mis();

        /**
         * Default Destructor.
         */
        virtual ~population_mis();

        /** 
         * Initialize the population.
         * 
         * @param config Config used for the population.
         * @param G Graph representation.
         */
        void init(MISConfig & config, graph_access & G);

        /**
         * Create a single individuum.
         * An individuum is created by producing a greedy initial solution.
         * The greedy algorithm is chosen at random.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Individuum to create/refine.
         */
        void create_individuum(MISConfig & config, graph_access & G, individuum_mis & ind);

        /**
         * Returns a certain individuum.
         *
         * @param id The id of the individuum.
         * @param ind The specific individuum.
         */
        void get_individuum(unsigned int id, individuum_mis & ind);

        /**
         * Returns a random individuum.
         *
         * @param ind The random individuum.
         */
        void get_random_individuum(individuum_mis & ind);

        /**
         * Returns a random individuum using tournament selection.
         *
         * @param config Config used for the population.
         * @param ind The random individuum.
         */
        void get_one_individuum_tournament(MISConfig & config, individuum_mis & ind);

        /**
         * Return two random individuals.
         *
         * @param first First random individuum.
         * @param second Second random individuum.
         */
        void get_two_random_individuals(individuum_mis & first, individuum_mis & second);

        /**
         * Return a number of random individuals.
         *
         * @param k The number of random individuals.
         * @param parents The resulting random individuals.
         */
        void get_random_individuals(unsigned int k, std::vector<individuum_mis> & parents);

        /**
         * Return two random individuals using tournament selection.
         *
         * @param config Config used for the population.
         * @param first First random individuum.
         * @param second Second random individuum.
         */
        void get_two_individuals_tournament(MISConfig & config, individuum_mis & first, individuum_mis & second);

        /**
         * Returns the individuum with the best solution.
         *
         * @param ind Individuum with the best solution.
         */
        void get_best_individuum(individuum_mis & ind);

        /**
         * Mutate a random individuum.
         * Mutation means an additional application of the ILS or GLP.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Mutated individuum.
         */
        void mutate_random(MISConfig & config, graph_access & G, individuum_mis & ind);

        /**
         * Mutate a given individuum.
         * Mutation means an additional application of the ILS or GLP.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Mutated individuum.
         */
        void mutate(MISConfig & config, graph_access & G, individuum_mis & ind);
        
        /**
         * Insert a individuum in the population.
         * The individuum is only inserted if its solution is
         * greater than an existing one or a certain threshold
         * is passed.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Individuum to insert.
         * @return True if the insertion was successful;
         */
        bool insert(MISConfig & config, graph_access & G, individuum_mis & ind);

        /**
         * Get the most similar individuum to a given one.
         * Similarity is measured using the hamming distance.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Individuum for which the most similar one should be found.
         * @param replacement The most similar individuum.
         * @return Whether or not the replacement is valid.
         */
        bool get_most_similar_replacement(MISConfig & config, graph_access & G, individuum_mis & ind, individuum_mis & replacement);

        /**
         * Replace an individuum by another.
         *
         * @param remove Individuum to remove.
         * @param insert Individuum to insert.
         */
        void replace(individuum_mis & remove, individuum_mis & insert);

        /**
         * Create a solution array for the current configuration of the graph.
         *
         * @param G Graph representation.
         * @param solution Node array containing the solution.
         * @return The size of the solution.
         */
        unsigned int create_solution(graph_access & G, NodeID *solution);

        /**
         * Takes the solution of the individuum and applies it to the graph.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Individuum whose solution should be applied.
         * @param secondary Whether or not the solution should be applied to the second partition index.
         */
        void set_mis_for_individuum(MISConfig & config, graph_access & G, individuum_mis & ind, bool secondary = false);

        /**
         * Resets the population by removing all individuals.
         */
        void extinction();

        /**
         * Prints the current state of the population.
         * Prints all individuals and their solution size.
         */
        void print(MISConfig & config);

        /**
         * Sets the maximum number of individuals.
         *
         * @param size Limit of individuals.
         */
        void set_population_size(unsigned int size);

        /**
         * Returns the average solution size.
         *
         * @return Average solution size.
         */
        double get_avg_solution_size();

        /**
         * Checks whether or not the population is full.
         * 
         * @return True if the population is full.
         */
        bool is_full();

        /**
         * Check if the individuum is a valid mis.
         *
         * @param config Config used for the population.
         * island.
         * @param G Graph representation.
         * @param ind Individuum to be checked.
         */
        bool is_mis(MISConfig & config, graph_access & G, individuum_mis & ind);

        /**
         * Check if the individuum is a valid vertex cover.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param ind Individuum to be checked.
         */
        bool is_vertex_cover(MISConfig & config, graph_access & G, individuum_mis & ind);

        /**
         * Get a fraction of the independent set nodes from the best individual.
         * Nodes are picked by increasing order of degree.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         * @param nodes Resulting set of nodes.
         */
        void get_best_individual_nodes(MISConfig & config, graph_access & G, std::vector<NodeID> & nodes);

        /**
         * Set the fraction of nodes that will be picked from the best individual.
         *
         * @param fraction Fraction of nodes.
         */
        void set_extract_fraction(double fraction);

        /**
         * Removes all individuals from the population and performs reinitialization.
         *
         * @param config Config used for the population.
         * @param G Graph representation.
         */
        void reset(MISConfig & config, graph_access & G);

    private:
        // Vector containing the individuals.
        std::vector<individuum_mis> internal_population;
        // Maximum size of the population.
        unsigned int population_size;
        // Number of insertions without a new solution.
        unsigned int insert_no_change;

        // Number of nodes that will be removed from the best individual.
        double remove_fraction;
};

#endif

