#ifndef JSSP_PROBLEM_H
#define JSSP_PROBLEM_H

#include <vector>
#include <string>

// Structure which defines our Operation
struct Operation {
    int job;
    int jobIndex;
    int machine;
    int duration;

    bool operator==(const Operation& other) const {
        return job == other.job && jobIndex == other.jobIndex && machine == other.machine;
    }
};

class JSSPProblem {
public:
    JSSPProblem(int numJobs = 0, int numMachines = 0);
    
    bool loadFromFiles(const std::string& processingTimeFile, const std::string& machineFile);
    
    int getNumJobs() const { return numJobs; }
    int getNumMachines() const { return numMachines; }
    const std::vector<std::vector<Operation>>& getJobs() const { return jobs; } // Returns the n x m Operations
    
    std::vector<Operation> getMachineOperations(int machine) const; // Returns the operations associated with the machine

private:
    int numJobs;
    int numMachines;
    std::vector<std::vector<Operation>> jobs;
    
    std::vector<std::vector<int>> readMatrix(const std::string& filename); // Method to read matrix from a file
};

#endif // JSSP_PROBLEM_H