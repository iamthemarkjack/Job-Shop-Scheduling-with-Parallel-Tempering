#include "jssp_problem.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>

JSSPProblem::JSSPProblem(int numJobs, int numMachines) 
    : numJobs(numJobs), numMachines(numMachines) {
    if (numJobs > 0) {
        jobs.resize(numJobs);
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
        
        if (machines.size() != numJobs) {
            std::cerr << "Error: Inconsistent number of jobs between files" << std::endl;
            return false;
        }
        
        for (int i = 0; i < numJobs; ++i) {
            if (processingTimes[i].size() != numMachines || machines[i].size() != numMachines) {
                std::cerr << "Error: Inconsistent number of machines in input files" << std::endl;
                return false;
            }
        }
        
        jobs.clear();
        jobs.resize(numJobs);
        
        for (int i = 0; i < numJobs; ++i) {
            for (int j = 0; j < numMachines; ++j) {
                int duration = processingTimes[i][j];
                int machine = machines[i][j];
                jobs[i].push_back({i, machine, duration, j});
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