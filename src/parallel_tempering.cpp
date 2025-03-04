#include "parallel_tempering.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

// Implementation of AdaptiveParallelTempering

AdaptiveParallelTempering::AdaptiveParallelTempering(const TSPProblem& problem, 
                                              int initialReplicas,
                                              double minTemp,
                                              double maxTemp,
                                              int swapInterval,
                                              int maxIterations,
                                              double targetSAP,
                                              double initialAlpha,
                                              double t0)
    : problem(problem), 
      numReplicas(initialReplicas),
      minTemp(minTemp),
      maxTemp(maxTemp),
      swapInterval(swapInterval),
      maxIterations(maxIterations),
      targetSAP(targetSAP),
      alpha(initialAlpha),
      alpha0(initialAlpha),
      t0(t0),
      totalSwapAttempts(0),
      thresholdDrop(0.05),
      thresholdInsert(0.05) {
    
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
    
    // Initialize swap statistics
    acceptedSwaps.resize(numReplicas - 1, 0);
    attemptedSwaps.resize(numReplicas - 1, 0);
    swapAcceptanceProbabilities.resize(numReplicas - 1, 0.0);
}

void AdaptiveParallelTempering::initializeTemperatures() {
    temperatures.resize(numReplicas);
    
    // Geometric progression from minTemp to maxTemp
    double ratio = std::pow(maxTemp / minTemp, 1.0 / (numReplicas - 1));
    
    for (int i = 0; i < numReplicas; ++i) {
        temperatures[i] = minTemp * std::pow(ratio, i);
    }
}

void AdaptiveParallelTempering::initializeReplicas() {
    replicas.resize(numReplicas);
    
    // Create random initial solutions for each replica
    for (int i = 0; i < numReplicas; ++i) {
        replicas[i] = Solution(problem);
        replicas[i].randomize(rngPool[i]);
    }
}

void AdaptiveParallelTempering::metropolisStep(Solution& solution, double temperature, std::mt19937& rng) {
    // Create a candidate solution by applying a random move
    Solution candidate = solution;
    
    // Apply a random move
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

double AdaptiveParallelTempering::calculateSwapAcceptance(const Solution& replica1, const Solution& replica2, 
    double temp1, double temp2) {
    // Calculate the energy (cost) difference
    double energy1 = replica1.getCost();
    double energy2 = replica2.getCost();

    // Calculate the Metropolis criterion for replica exchange
    double deltaE = (1.0 / temp1 - 1.0 / temp2) * (energy1 - energy2);

    // Return the acceptance probability
    return std::min(1.0, std::exp(deltaE));
}

void AdaptiveParallelTempering::attemptSwaps() {
    // Random number generator for swap decisions
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    // Swap only between adjacent temperatures
    for (int i = 0; i < numReplicas - 1; ++i) {
        totalSwapAttempts++;
        attemptedSwaps[i]++;
        
        double acceptProb = calculateSwapAcceptance(
            replicas[i], replicas[i+1], 
            temperatures[i], temperatures[i+1]
        );
        
        if (dist(rng) < acceptProb) {
            // Accept swap
            std::swap(replicas[i], replicas[i+1]);
            acceptedSwaps[i]++;
        }
    }
}

void AdaptiveParallelTempering::updateSAPs() {
    // Update swap acceptance probabilities
    for (int i = 0; i < numReplicas - 1; ++i) {
        if (attemptedSwaps[i] > 0) {
            swapAcceptanceProbabilities[i] = static_cast<double>(acceptedSwaps[i]) / attemptedSwaps[i];
        }
    }
}

void AdaptiveParallelTempering::updateDampingFactor() {
    // Update damping factor according to equation (5)
    alpha = alpha0 * (t0 / (t0 + totalSwapAttempts));
}

void AdaptiveParallelTempering::adjustTemperatures() {
    // Implement Algorithm 1 from the paper
    std::vector<double> newTemperatures = temperatures;
    
    for (int i = 0; i < numReplicas - 1; ++i) {
        // Calculate gap between actual SAP and target SAP
        double gap = swapAcceptanceProbabilities[i] - targetSAP;
        
        if (gap > 0) {
            // SAP is too high, increase temperature difference
            double deltaT = (temperatures[i+1] - temperatures[i]) * alpha;
            newTemperatures[i+1] += deltaT;
        } else if (gap < 0) {
            // SAP is too low, decrease temperature difference
            double deltaT = (temperatures[i+1] - temperatures[i]) * alpha;
            newTemperatures[i+1] -= deltaT;
        }
    }
    
    // Update temperatures
    temperatures = newTemperatures;
}

void AdaptiveParallelTempering::manageReplicas() {
    // Check if we need to drop or insert replicas
    
    // If numReplicas < 2, we can't perform any more operations
    if (numReplicas < 2) return;
    
    // Get SAP between the two highest temperatures
    double highestSAP = swapAcceptanceProbabilities[numReplicas - 2];
    
    if (highestSAP > targetSAP + thresholdDrop) {
        // Drop the highest temperature replica
        temperatures.pop_back();
        replicas.pop_back();
        acceptedSwaps.pop_back();
        attemptedSwaps.pop_back();
        swapAcceptanceProbabilities.pop_back();
        numReplicas--;
        
        std::cout << "Dropped replica. New count: " << numReplicas << std::endl;
    } 
    else if (highestSAP < targetSAP - thresholdInsert) {
        // Insert new replica between two highest temperatures
        double newTemp = 0.5 * (temperatures[numReplicas-1] + temperatures[numReplicas-2]);
        temperatures.push_back(newTemp);
        
        // Create new replica by cloning the second highest one
        replicas.push_back(replicas[numReplicas-2]);
        
        // Initialize statistics for the new pair
        acceptedSwaps.push_back(0);
        attemptedSwaps.push_back(0);
        swapAcceptanceProbabilities.push_back(0.0);
        
        numReplicas++;
        
        std::cout << "Inserted replica. New count: " << numReplicas << std::endl;
    }
}

Solution AdaptiveParallelTempering::solve(const std::string& historyFilename) {
    // Set the number of OpenMP threads
    omp_set_num_threads(numThreads);

    // Main algorithm loop
    int iteration = 0;
    
    while (iteration < maxIterations) {
        // Run MCMC steps for each replica in parallel
        #pragma omp parallel for num_threads(numThreads)
        for (int i = 0; i < numReplicas; ++i) {
            int threadId = omp_get_thread_num();
            
            // For each iteration, perform multiple Metropolis steps
            // to ensure good mixing at that temperature
            for (int step = 0; step < 100; step++) {
                metropolisStep(replicas[i], temperatures[i], rngPool[threadId]);
            }
        }
        
        // Attempt swaps between adjacent temperatures
        #pragma omp critical
        {
            attemptSwaps();
        }
        
        // Every N_swap iterations, adjust temperatures and manage replicas
        if (iteration % swapInterval == 0 && iteration > 0) {
            updateSAPs();
            updateDampingFactor();
            adjustTemperatures();
            manageReplicas();
        }

        // Record the current state
        recordCurrentState(iteration);

        // Add progress reporting
        if (iteration % (maxIterations / 10) == 0) {
            std::cout << "Iteration " << iteration << "/" << maxIterations 
                      << ", Best solution: " << replicas[0].getCost() << std::endl;
        }
        
        for (int i = 0; i < numReplicas - 1; ++i) {
            std::cout << swapAcceptanceProbabilities[i];
        }

        iteration++;
    }

    // Save history to file if filename is provided
    if (!historyFilename.empty()) {
        saveHistoryToFile(historyFilename);
    }
    
    // Return the best solution (should be at lowest temperature)
    Solution bestSolution = replicas[0];
    double bestCost = bestSolution.getCost();
    
    for (int i = 1; i < numReplicas; ++i) {
        double cost = replicas[i].getCost();
        if (cost < bestCost) {
            bestCost = cost;
            bestSolution = replicas[i];
        }
    }
    
    return bestSolution;
}

void AdaptiveParallelTempering::recordCurrentState(int iteration) {
    // Record the current temperatures
    temperatureHistory.push_back(temperatures);
    
    // Record current energies
    std::vector<double> currentEnergies;
    for (int i = 0; i < numReplicas ; ++i){
        currentEnergies.push_back(replicas[i].getCost());
    }
    energyHistory.push_back(currentEnergies);

    // Record the current SAPs
    if (iteration % swapInterval == 0 && iteration > 0) {
        sapHistory.push_back(swapAcceptanceProbabilities);
        iterationHistory.push_back(iteration);
    }
}

void AdaptiveParallelTempering::saveHistoryToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    
    // Write header
    file << "Iteration,ReplicaIndex,Temperature,Energy,SAP" << std::endl;
    
    // Write temperature data for all iterations
    for (size_t i = 0; i < temperatureHistory.size(); ++i) {
        int iteration = i; // Assuming we record every iteration
        
        for (size_t j = 0; j < temperatureHistory[i].size(); ++j) {
            file << iteration << "," << j << "," << temperatureHistory[i][j] << ",";

            // Write energy data
            if (j < energyHistory[i].size()) {
                file << energyHistory[i][j] << ",";
            } else {
                file << "NA,"; // Shouldn't happen unless replicas were added/removed
            }
            
            // Find the corresponding SAP data
            bool sapFound = false;
            for (size_t k = 0; k < iterationHistory.size(); ++k) {
                if (iterationHistory[k] == iteration) {
                    // SAP is only defined between replicas, so we have one less SAP than temperatures
                    if (j < sapHistory[k].size()) {
                        file << sapHistory[k][j];
                    } else {
                        file << "NA"; // For the last replica, there's no SAP
                    }
                    sapFound = true;
                    break;
                }
            }
            
            if (!sapFound) {
                file << "NA"; // No SAP recorded for this iteration
            }
            
            file << std::endl;
        }
    }
    
    file.close();
    std::cout << "History data saved to " << filename << std::endl;
}