#include "parallel_tempering.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>

AdaptiveParallelTempering::AdaptiveParallelTempering(const JSSPProblem& problem,
                                                    int numReplicas,
                                                    double minTemp,
                                                    double rho,
                                                    int maxIterations,
                                                    double target_SAP,
                                                    double C,
                                                    double eta,
                                                    std::string hisPrefix)
    : problem(problem),
      numReplicas(numReplicas),
      minTemp(minTemp),
      rho(rho),
      maxIterations(maxIterations),
      target_SAP(target_SAP),
      C(C),
      eta(eta),
      hisPrefix(hisPrefix) {

    unsigned seed = std::random_device{}();
    for (int i = 0; i < numReplicas; ++i){
        rngPool.emplace_back(seed + i);
    }

    replicas.resize(numReplicas);

    rhos.resize(numReplicas - 1, rho);
    betas.resize(numReplicas, 1 / minTemp);

    attemptedSwaps.resize(numReplicas - 1, 0);
    acceptedSwaps.resize(numReplicas - 1, 0);
    ratios.resize(numReplicas - 1, 0.0);

    betaHistory.resize(numReplicas);
    energyHistory.resize(numReplicas);
    SAPHistory.resize(numReplicas - 1);

    initialize();
}

void AdaptiveParallelTempering::initialize(){
    for(int i = 1; i < numReplicas; ++i){
        betas[i] = betas[i - 1] / 10.0;
    }
    
    for (int i = 0; i < numReplicas; ++i){
        replicas[i] = JSSPSolution(problem);
        replicas[i].randomize(rngPool[i]);
    }
}

void AdaptiveParallelTempering::metropolisStep(JSSPSolution& solution, double beta, std::mt19937_64& rng){
    JSSPSolution candidate = solution;
    candidate.generateNeighbor(rng);

    double currentCost = solution.getMakespan();
    double candidateCost = candidate.getMakespan();
    double delta = candidateCost - currentCost;

    std::uniform_real_distribution<double> dist(0.0, 1.0);
    if (dist(rng) < std::min(1.0, std::exp(-delta * beta))){
        solution = candidate;
    }
}

double AdaptiveParallelTempering::calculateSAP(const JSSPSolution& replica1, const JSSPSolution& replica2, double beta1, double beta2){
    double E1 = replica1.getMakespan();
    double E2 = replica2.getMakespan();

    double exponent = (beta1 - beta2) * (E1 - E2);

    return std::min(1.0 , std::exp(exponent));
}

void AdaptiveParallelTempering::attemptSwap(){
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<int> int_dist(0, numReplicas - 2);
    std::uniform_real_distribution<double> real_dist(0.0, 1.0);

    int k = int_dist(rng);

    double acceptProb = calculateSAP(replicas[k], replicas[k + 1], betas[k], betas[k + 1]);
    attemptedSwaps[k]++;

    if (real_dist(rng) < acceptProb){
        std::swap(replicas[k], replicas[k + 1]);
        acceptedSwaps[k]++;
    }

    for (int i = 0; i < numReplicas - 1; ++i){
        double val = attemptedSwaps[i] > 0 ? static_cast<double>(acceptedSwaps[i]) / attemptedSwaps[i] : 0.0;
        ratios[i] = val;
    }
}

void AdaptiveParallelTempering::updateBetas(int iteration){
    for (int i = 0; i < numReplicas - 1; ++i){
        double SAP = calculateSAP(replicas[i], replicas[i + 1], betas[i], betas[i + 1]);
        rhos[i] += gamma(iteration) * (target_SAP - SAP);
    }

    for (int i = 0; i < numReplicas - 1; ++i){
        betas[i + 1] = betas[i] * (1 / (1 + std::exp( -rhos[i])));

    }
}

double AdaptiveParallelTempering::gamma(int t){
    return C * std::pow(1 + t, -eta);
}

JSSPSolution AdaptiveParallelTempering::solve() {
    int iteration = 0;
    
    while (iteration < maxIterations) {
        for (int i = 0; i < numReplicas; ++i) {
            metropolisStep(replicas[i], betas[i], rngPool[i]);
        }
        
        attemptSwap();
        updateBetas(iteration);

        iteration += 1;

        for (int i = 0; i < numReplicas; ++i) {
            betaHistory[i].push_back(betas[i]);
            energyHistory[i].push_back(replicas[i].getMakespan());
            if (i < numReplicas - 1) SAPHistory[i].push_back(ratios[i]);
        }

        if (iteration % (maxIterations / 10) == 0) {
            std::cout << "Iteration " << iteration << "/" << maxIterations 
                      << ", best solution: " << replicas[0].getMakespan() << std::endl;
        }
    }

    std::cout << "Final Rho : ";
    for (double rho_val : rhos){
        std::cout <<  rho_val << " ";
    }
    std::cout << std::endl;

    saveHistory(hisPrefix);
        
    JSSPSolution bestSolution = replicas[0];
    double bestCost = bestSolution.getMakespan();
    for (int i = 1; i < numReplicas; ++i) {
        double cost = replicas[i].getMakespan();
        if (cost < bestCost) {
            bestCost = cost;
            bestSolution = replicas[i];
        }
    }
    return bestSolution;
}

void AdaptiveParallelTempering::saveHistory(const std::string& prefix) const {
    // Save beta history
    std::ofstream betaFile("histories/" + prefix + "_beta_history.csv");
    betaFile << "Iteration";
    for (int i = 0; i < numReplicas; ++i) {
        betaFile << ",Replica" << i;
    }
    betaFile << std::endl;
    
    for (size_t iter = 0; iter < betaHistory[0].size(); ++iter) {
        betaFile << iter;
        for (int i = 0; i < numReplicas; ++i) {
            betaFile << "," << betaHistory[i][iter];
        }
        betaFile << std::endl;
    }
    betaFile.close();
    
    // Save energy history
    std::ofstream energyFile("histories/" + prefix + "_energy_history.csv");
    energyFile << "Iteration";
    for (int i = 0; i < numReplicas; ++i) {
        energyFile << ",Replica" << i;
    }
    energyFile << std::endl;
    
    for (size_t iter = 0; iter < energyHistory[0].size(); ++iter) {
        energyFile << iter;
        for (int i = 0; i < numReplicas; ++i) {
            energyFile << "," << energyHistory[i][iter];
        }
        energyFile << std::endl;
    }
    energyFile.close();
    
    // Save SAP history
    std::ofstream sapFile("histories/" + prefix + "_sap_history.csv");
    sapFile << "Iteration";
    for (int i = 0; i < numReplicas - 1; ++i) {
        sapFile << ",Pair" << i << "_" << i+1;
    }
    sapFile << std::endl;
    
    for (size_t iter = 0; iter < SAPHistory[0].size(); ++iter) {
        sapFile << iter;
        for (int i = 0; i < numReplicas - 1; ++i) {
            sapFile << "," << SAPHistory[i][iter];
        }
        sapFile << std::endl;
    }
    sapFile.close();
}