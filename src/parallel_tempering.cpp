#include "parallel_tempering.h"
#include <cmath>
#include <algorithm>
#include <iostream>

ParallelTempering::ParallelTempering(const TSPProblem& problem, 
                                     int numReplicas,
                                     double minTemp,
                                     double maxTemp,
                                     int swapInterval,
                                     int maxIterations) 
    : problem(problem),
      numReplicas(numReplicas),
      minTemp(minTemp),
      maxTemp(maxTemp),
      swapInterval(swapInterval),
      maxIterations(maxIterations) {
    
    // Set default number of threads to number of replicas if not overridden later
    numThreads = numReplicas;
    
    // Initialize RNG pool (one per thread)
    unsigned seed = std::random_device{}();
    for (int i = 0; i < numReplicas; ++i) {
        rngPool.emplace_back(seed + i);  // Different seed for each thread
    }
    
    // Initialize temperatures and replicas
    initializeTemperatures();
    initializeReplicas();
}

void ParallelTempering::initializeTemperatures() {
    temperatures.resize(numReplicas);
    
    // Create a geometric progression of temperatures
    double factor = std::pow(maxTemp / minTemp, 1.0 / (numReplicas - 1));
    
    for (int i = 0; i < numReplicas; ++i) {
        temperatures[i] = minTemp * std::pow(factor, i);
    }
}

void ParallelTempering::initializeReplicas() {
    replicas.resize(numReplicas);
    
    // Create random initial solutions for each replica
    for (int i = 0; i < numReplicas; ++i) {
        replicas[i] = Solution(problem);
        replicas[i].randomize(rngPool[i]);
    }
}

Solution ParallelTempering::solve() {
    // Set the number of OpenMP threads
    omp_set_num_threads(numThreads);
    
    // Main loop
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Parallel Metropolis steps for each replica
        #pragma omp parallel for
        for (int i = 0; i < numReplicas; ++i) {
            int threadId = omp_get_thread_num();
            std::mt19937& rng = rngPool[threadId];
            
            // Perform multiple Metropolis steps before swap attempt
            for (int step = 0; step < swapInterval; ++step) {
                metropolisStep(replicas[i], temperatures[i], rng);
            }
        }
        
        // Attempt to swap replicas between adjacent temperatures
        attemptSwaps();
        
        // Optional: Add progress reporting
        if (iter % (maxIterations / 10) == 0) {
            std::cout << "Iteration " << iter << "/" << maxIterations 
                      << ", Best solution: " << replicas[0].getCost() << std::endl;
        }
    }
    
    // Find the best solution across all replicas
    Solution bestSolution = replicas[0];
    for (int i = 1; i < numReplicas; ++i) {
        if (replicas[i].getCost() < bestSolution.getCost()) {
            bestSolution = replicas[i];
        }
    }
    
    return bestSolution;
}

void ParallelTempering::metropolisStep(Solution& solution, double temperature, std::mt19937& rng) {
    // Create a candidate solution by applying a random move
    Solution candidate = solution;
    
    // Apply a random move (e.g., 2-opt, insertion, etc.)
    candidate.applyRandomMove(rng);
    
    // Calculate the cost difference
    double currentCost = solution.getCost();
    double candidateCost = candidate.getCost();
    double delta = candidateCost - currentCost;
    
    // Accept or reject the move based on Metropolis criterion
    if (delta <= 0) {
        // Always accept if the new solution is better
        solution = candidate;
    } else {
        // Accept with probability exp(-delta/T) if the new solution is worse
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        if (dist(rng) < std::exp(-delta / temperature)) {
            solution = candidate;
        }
    }
}

void ParallelTempering::attemptSwaps() {
    // Random number generator for swap decisions
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    // Try to swap adjacent replicas
    for (int i = 0; i < numReplicas - 1; i += 2) {
        double acceptance = calculateSwapAcceptance(replicas[i], replicas[i+1], 
                                                  temperatures[i], temperatures[i+1]);
        
        if (dist(rng) < acceptance) {
            std::swap(replicas[i], replicas[i+1]);
        }
    }
    
    // Try to swap odd-indexed pairs
    for (int i = 1; i < numReplicas - 1; i += 2) {
        double acceptance = calculateSwapAcceptance(replicas[i], replicas[i+1], 
                                                  temperatures[i], temperatures[i+1]);
        
        if (dist(rng) < acceptance) {
            std::swap(replicas[i], replicas[i+1]);
        }
    }
}

double ParallelTempering::calculateSwapAcceptance(const Solution& replica1, const Solution& replica2, 
                                                double temp1, double temp2) {
    // Calculate the energy (cost) difference
    double energy1 = replica1.getCost();
    double energy2 = replica2.getCost();
    
    // Calculate the Metropolis criterion for replica exchange
    double deltaE = (1.0 / temp1 - 1.0 / temp2) * (energy1 - energy2);
    
    // Return the acceptance probability
    return std::min(1.0, std::exp(deltaE));
}