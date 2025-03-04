#ifndef ADAPTIVE_PARALLEL_TEMPERING_H
#define ADAPTIVE_PARALLEL_TEMPERING_H

#include "tsp_problem.h"
#include "solution.h"

#include <vector>
#include <functional>
#include <random>
#include <omp.h>

class AdaptiveParallelTempering {
public:
    // Constructor with parameters for adaptive tempering
    AdaptiveParallelTempering(const TSPProblem& problem, 
                      int initialReplicas = 16,     // Start with more replicas than needed
                      double minTemp = 1.0,
                      double maxTemp = 10.0,
                      int swapInterval = 100,
                      int maxIterations = 10000,
                      double targetSAP = 0.23,       // Target swap acceptance probability
                      double initialAlpha = 1.0,    // Initial damping factor
                      double t0 = 1000);           // Time constant for damping decay

    // Run the adaptive parallel tempering algorithm
    Solution solve(const std::string& historyFilename = "");

    // Getters
    std::vector<Solution> getReplicas() const { return replicas; }
    std::vector<double> getTemperatures() const { return temperatures; }
    std::vector<double> getSAPs() const { return swapAcceptanceProbabilities; }
    
    // Setters
    void setNumThreads(int threads) { numThreads = threads; }
    
    // Method to save history to file
    void saveHistoryToFile(const std::string& filename) const;

private:
    // TSP problem instance
    const TSPProblem& problem;
    
    // Replicas (solutions at different temperatures)
    std::vector<Solution> replicas;
    
    // Temperature ladder
    std::vector<double> temperatures;
    
    // Statistics for adaptive algorithm
    std::vector<int> acceptedSwaps;    // Number of accepted swaps for each temperature pair
    std::vector<int> attemptedSwaps;   // Number of attempted swaps for each temperature pair
    std::vector<double> swapAcceptanceProbabilities; // SAPs between adjacent temperatures

    // For tracking temperatures, SAP and Energy 
    std::vector<std::vector<double>> temperatureHistory;
    std::vector<std::vector<double>> sapHistory;
    std::vector<std::vector<double>> energyHistory;
    std::vector<int> iterationHistory;
    
    // Algorithm parameters
    int numReplicas;
    int numThreads;
    double minTemp;
    double maxTemp;
    int swapInterval;     // N_swap
    int maxIterations;
    
    // Adaptive parameters
    double targetSAP;     // Target swap acceptance probability
    double alpha;         // Current damping factor
    double alpha0;        // Initial damping factor
    double t0;            // Time constant for damping decay
    int totalSwapAttempts; // t in equation (5)
    double thresholdDrop;  // Threshold for dropping replicas
    double thresholdInsert; // Threshold for inserting replicas
    
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
    
    // Adjust temperatures to achieve uniform SAP
    void adjustTemperatures();
    
    // Manage replica insertion/removal
    void manageReplicas();
    
    // Update damping factoor
    void updateDampingFactor();
    
    // Update SAP statistics
    void updateSAPs();

    // Method to record current state
    void recordCurrentState(int iteration);
};

#endif // ADAPTIVE_PARALLEL_TEMPERING_H