#include "jssp_problem.h"
#include "solution.h"
#include <iostream>
#include <random>
#include <iomanip>

void printJobOrdering(const JSSPSolution& solution) {
    const std::vector<std::vector<int>>& machineJobOrder = solution.getMachineJobOrder();
    
    for (size_t m = 0; m < machineJobOrder.size(); ++m) {
        std::cout << "Machine " << m << ": ";
        for (size_t j = 0; j < machineJobOrder[m].size(); ++j) {
            std::cout << machineJobOrder[m][j];
            if (j < machineJobOrder[m].size() - 1) {
                std::cout << " -> ";
            }
        }
        std::cout << std::endl;
    }
}

int main() {
    JSSPProblem problem;
    if (!problem.loadFromFiles("../data/instances/processing_times.txt", "../data/instances/machines.txt")) {
        std::cerr << "Failed to load JSSP instance." << std::endl;
        return 1;
    }
    
    // Initialize the solution
    JSSPSolution solution(problem);
    
    // Print initial job ordering on machines
    std::cout << "Initial job ordering on machines:" << std::endl;
    printJobOrdering(solution);
    
    // Print initial makespan
    std::cout << "Initial Makespan: " << solution.getMakespan() << std::endl;
    
    // Set up random number generator
    std::random_device rd;
    std::mt19937_64 rng(rd());
    
    // Generate and evaluate 10 neighbors
    for (int i = 0; i < 100; ++i) {
        JSSPSolution neighbor = solution.generateNeighbor(rng);
        std::cout << "Neighbor " << i + 1 << " Makespan: " << neighbor.getMakespan() << std::endl;
        std::cout << "Job ordering:" << std::endl;
        printJobOrdering(neighbor);
        std::cout << std::endl;
    }
    
    return 0;
}