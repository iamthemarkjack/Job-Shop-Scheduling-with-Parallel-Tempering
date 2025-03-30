#include "jssp_problem.h"
#include "solution.h"
#include <iostream>
#include <random>
#include <iomanip>

void printJobOrdering(const JSSPSolution& solution) {
    const std::vector<std::vector<Operation>>& machineJobOrders = solution.getMachineJobOrders();
    
    for (size_t m = 0; m < machineJobOrders.size(); ++m) {
        std::cout << "Machine " << m << ": ";
        for (size_t j = 0; j < machineJobOrders[m].size(); ++j) {
            std::cout << machineJobOrders[m][j].job;
            if (j < machineJobOrders[m].size() - 1) {
                std::cout << " -> ";
            }
        }
        std::cout << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <processing_times_file> <machines_file>" << std::endl;
        return 1;
    }

    std::string processing_times_file = argv[1];
    std::string machines_file = argv[2];

    JSSPProblem problem;
    if (!problem.loadFromFiles(processing_times_file, machines_file)) {
        std::cerr << "Failed to load JSSP instance." << std::endl;
        return 1;
    }

    JSSPSolution solution(problem);
    
    std::cout << "Initial job ordering on machines:" << std::endl;
    printJobOrdering(solution);
    
    std::cout << "Initial Makespan: " << solution.getMakespan() << std::endl;
    std::cout << std::endl;
    
    std::random_device rd;
    std::mt19937_64 rng(rd());
    
    for (int i = 0; i < 10; ++i) {
        JSSPSolution neighbor = solution.generateNeighbor(rng);
        std::cout << "Neighbor " << i + 1 << " Makespan: " << neighbor.getMakespan() << std::endl;
        std::cout << "Job ordering:" << std::endl;
        printJobOrdering(neighbor);
        std::cout << std::endl;
    }
    
    return 0;
}