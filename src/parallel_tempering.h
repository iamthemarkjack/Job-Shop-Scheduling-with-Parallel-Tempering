#ifndef PARALLEL_TEMPERING_H
#define PARALLEL_TEMPERING_H

#include <vector>
#include <functional>
#include <random>
#include <omp.h>
#include "tsp_problem.h"
#include "solution.h"

class ParallelTempering {
public:
    // Constructor
    ParallelTempering(const TSPProblem& problem, 
                      int numReplicas = 8,
                      double minTemp = 0.1,
                      double maxTemp = 100.0,
                      int swapInterval = 100,
                      int maxIterations = 10000);

    // Run the parallel tempering algorithm
    Solution solve();

    // Getters
    std::vector<Solution> getReplicas() const { return replicas; }
    std::vector<double> getTemperatures() const { return temperatures; }
    
    // Setters
    void setNumThreads(int threads) { numThreads = threads; }

private:
    // TSP problem instance
    const TSPProblem& problem;
    
    // Replicas (solutions at different temperatures)
    std::vector<Solution> replicas;
    
    // Temperature ladder
    std::vector<double> temperatures;
    
    // Algorithm parameters
    int numReplicas;
    int numThreads;
    double minTemp;
    double maxTemp;
    int swapInterval;
    int maxIterations;
    
    // Random generators (one per thread)
    std::vector<std::mt19937> rngPool;
    
    // Initialize temperature ladder (geometric progression)
    void initializeTemperatures();
    
    // Initialize replicas with random solutions
    void initializeReplicas();
    
    // Perform a Metropolis step at a given temperature
    void metropolisStep(Solution& solution, double temperature, std::mt19937& rng);
    
    // Attempt to swap configurations between adjacent temperatures
    void attemptSwaps();
    
    // Calculate the acceptance probability for replica exchange
    double calculateSwapAcceptance(const Solution& replica1, const Solution& replica2, 
                                  double temp1, double temp2);
};

#endif // PARALLEL_TEMPERING_H