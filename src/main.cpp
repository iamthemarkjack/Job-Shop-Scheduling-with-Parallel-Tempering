#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>
#include "parallel_tempering.h"
#include "tsp_problem.h"
#include "solution.h"

void printHelp() {
    std::cout << "TSP Solver using Parallel Tempering\n";
    std::cout << "Usage:\n";
    std::cout << "  --file <filename>       Load TSP instance from file\n";
    std::cout << "  --random <num_cities>   Generate random TSP instance\n";
    std::cout << "  --replicas <num>        Number of temperature replicas (default: 8)\n";
    std::cout << "  --threads <num>         Number of OpenMP threads (default: same as replicas)\n";
    std::cout << "  --min-temp <value>      Minimum temperature (default: 0.1)\n";
    std::cout << "  --max-temp <value>      Maximum temperature (default: 100.0)\n";
    std::cout << "  --iterations <num>      Maximum iterations (default: 10000)\n";
    std::cout << "  --swap-interval <num>   Swap attempt interval (default: 100)\n";
    std::cout << "  --seed <num>            Random seed (default: from system clock)\n";
    std::cout << "  --output <filename>     Save best tour to file\n";
    std::cout << "  --help                  Display this help message\n";
}

void saveTourToFile(const std::string& filename, const Solution& solution) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }
    
    const auto& tour = solution.getTour();
    
    outFile << "Tour length: " << std::fixed << std::setprecision(2) << solution.getCost() << "\n";
    outFile << "Tour: ";
    for (size_t i = 0; i < tour.size(); ++i) {
        outFile << tour[i];
        if (i < tour.size() - 1) {
            outFile << " ";
        }
    }
    outFile << "\n";
    
    std::cout << "Tour saved to " << filename << "\n";
}

int main(int argc, char* argv[]) {
    // Default parameters
    std::string inputFile;
    int randomSize = 0;
    int numReplicas = 8;
    int numThreads = 0;  // 0 means use same as replicas
    double minTemp = 0.1;
    double maxTemp = 100.0;
    int maxIterations = 10000;
    int swapInterval = 100;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::string outputFile;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help") {
            printHelp();
            return 0;
        } else if (arg == "--file" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (arg == "--random" && i + 1 < argc) {
            randomSize = std::stoi(argv[++i]);
        } else if (arg == "--replicas" && i + 1 < argc) {
            numReplicas = std::stoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            numThreads = std::stoi(argv[++i]);
        } else if (arg == "--min-temp" && i + 1 < argc) {
            minTemp = std::stod(argv[++i]);
        } else if (arg == "--max-temp" && i + 1 < argc) {
            maxTemp = std::stod(argv[++i]);
        } else if (arg == "--iterations" && i + 1 < argc) {
            maxIterations = std::stoi(argv[++i]);
        } else if (arg == "--swap-interval" && i + 1 < argc) {
            swapInterval = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = std::stoul(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            outputFile = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            printHelp();
            return 1;
        }
    }
    
    // Check that either a file is provided or a random instance size
    if (inputFile.empty() && randomSize <= 0) {
        std::cerr << "Error: Must specify either --file or --random\n";
        printHelp();
        return 1;
    }
    
    // Set threads to replicas if not specified
    if (numThreads <= 0) {
        numThreads = numReplicas;
    }
    
    // Create the TSP problem
    TSPProblem problem;
    
    try {
        if (!inputFile.empty()) {
            std::cout << "Loading TSP instance from file: " << inputFile << "\n";
            if (!problem.loadFromFile(inputFile)) {
                std::cerr << "Error: Failed to load TSP instance from file.\n";
                return 1;
            }
            std::cout << "Loaded " << problem.getNumCities() << " cities.\n";
        } else {
            std::cout << "Generating random TSP instance with " << randomSize << " cities.\n";
            problem.generateRandom(randomSize);
        }
        
        // Create and configure the solver
        std::cout << "Initializing parallel tempering solver with " << numReplicas << " replicas.\n";
        std::cout << "Using " << numThreads << " OpenMP threads.\n";
        std::cout << "Temperature range: " << minTemp << " to " << maxTemp << "\n";
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        ParallelTempering solver(problem, numReplicas, minTemp, maxTemp, swapInterval, maxIterations);
        solver.setNumThreads(numThreads);
        
        // Run the solver
        std::cout << "Starting parallel tempering algorithm...\n";
        Solution bestSolution = solver.solve();
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        
        // Print results
        std::cout << "Optimization completed in " << duration / 1000.0 << " seconds.\n";
        std::cout << "Best tour length: " << std::fixed << std::setprecision(2) << bestSolution.getCost() << "\n";
        
        // Save results if requested
        if (!outputFile.empty()) {
            saveTourToFile(outputFile, bestSolution);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}