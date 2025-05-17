#include "jssp_problem.h"

#include <fstream>
#include <sstream>
#include <iostream>

JSSPProblem::JSSPProblem(int numJobs, int numMachines) 
    : numJobs(numJobs), numMachines(numMachines) {
    if (numJobs > 0) {
        jobs.resize(numJobs); // Size of jobs is numJobs x numMachines
    }
}

bool JSSPProblem::loadFromFiles(const std::string& processingTimeFile, const std::string& machineFile) {
    try {
        std::vector<std::vector<int>> processingTimes = readMatrix(processingTimeFile);
        std::vector<std::vector<int>> machines = readMatrix(machineFile);
        
        if (processingTimes.empty() || machines.empty()) {
            std::cerr << "Error: Empty input files" << std::endl;
            return false;
        }
        
        numJobs = processingTimes.size();
        numMachines = processingTimes[0].size();
        
        // Check if the number of jobs match between the files
        if (machines.size() != numJobs) {
            std::cerr << "Error: Inconsistent number of jobs between files" << std::endl;
            return false;
        }
        
        // Check if each jobs has the same number of operations equal to numMachines in both files
        for (int i = 0; i < numJobs; ++i) {
            if (processingTimes[i].size() != numMachines || machines[i].size() != numMachines) {
                std::cerr << "Error: Inconsistent number of machines in input files" << std::endl;
                return false;
            }
        }
        
        jobs.clear();
        jobs.resize(numJobs);
        
        // Load the Operations onto jobs
        for (int i = 0; i < numJobs; ++i) {
            for (int j = 0; j < numMachines; ++j) {
                int duration = processingTimes[i][j];
                int machine = machines[i][j];
                jobs[i].push_back({i, j, machine, duration});
            }
        }
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error loading JSSP problem: " << e.what() << std::endl;
        return false;
    }
}

std::vector<Operation> JSSPProblem::getMachineOperations(int machine) const {
    std::vector<Operation> machineOps;
    
    // Iteration over every Operations in jobs to fetch the Operations scheduled on machine
    for (int i = 0; i < numJobs; ++i) {
        for (const auto& op : jobs[i]) {
            if (op.machine == machine) {
                machineOps.push_back(op);
            }
        }
    }
    
    return machineOps;
}

std::vector<std::vector<int>> JSSPProblem::readMatrix(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<int>> matrix;
    std::string line;

    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        throw std::runtime_error("File not found");
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<int> row;
        int value;
        while (ss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}