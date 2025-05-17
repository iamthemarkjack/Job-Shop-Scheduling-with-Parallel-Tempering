#ifndef ADAPTIVE_PARALLEL_TEMPERING_H
#define ADAPTIVE_PARALLEL_TEMPERING_H

#include "jssp_problem.h"
#include "solution.h"

#include <vector>
#include <functional>
#include <random>

class AdaptiveParallelTempering{
public:
    // Constructor
    AdaptiveParallelTempering(const JSSPProblem& problem,
                            int numReplicas = 8,
                            double minTemp = 1.0,
                            double rho  = 1.0,
                            double maxTimeSeconds = 60,
                            int N_sweep = 100,
                            double target_SAP = 0.234,
                            double C = 1,
                            double eta = 2/3,
                            std::string hisPrefix = "pt_run");
    
    JSSPSolution solve();

private:
        const JSSPProblem& problem;
        int numReplicas;
        double minTemp;
        double rho;
        double maxTimeSeconds;
        int N_sweep;

        double target_SAP;
        double C;
        double eta;
        std::string hisPrefix;

        double bestMakespan;
        JSSPSolution bestSolution;

        std::vector<JSSPSolution> replicas;
        std::vector<JSSPSolution> bestReplicas;
        std::vector<double> rhos;
        std::vector<double> betas;

        std::vector<int> acceptedSwaps;
        std::vector<int> attemptedSwaps;
        std::vector<double> ratios;

        std::vector<std::vector<double>> betaHistory;
        std::vector<std::vector<double>> energyHistory;
        std::vector<std::vector<double>> SAPHistory;

        // Random generators (one per thread)
        std::vector<std::mt19937_64> rngPool;

        void initialize();

        void metropolisStep(JSSPSolution& solution, double beta, std::mt19937_64& rng);

        void attemptSwap();

        void updateBetas(int iteration);

        void updateBestSolutions();

        double gamma(int iteration); // Decay function (Polynomial Decay)

        double calculateSAP(const JSSPSolution& replica1, const JSSPSolution& replica2, double beta1, double beta);

        void saveHistory(const std::string& prefix) const;
};

#endif // ADAPTIVE_PARALLEL_TEMPERING_H