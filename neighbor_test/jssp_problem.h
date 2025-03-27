#ifndef JSSP_PROBLEM_H
#define JSSP_PROBLEM_H

#include <vector>
#include <string>
#include <unordered_map>
#include <set>

struct Operation {
    int job;
    int machine;
    int duration;
    int index; // Position within job
};

class JSSPProblem {
public:
    JSSPProblem(int numJobs = 0, int numMachines = 0);
    
    bool loadFromFiles(const std::string& processingTimeFile, const std::string& machineFile);
    
    int getNumJobs() const { return numJobs; }
    int getNumMachines() const { return numMachines; }
    const std::vector<std::vector<Operation>>& getJobs() const { return jobs; }
    
    std::vector<Operation> getMachineOperations(int machine) const;

private:
    int numJobs;
    int numMachines;
    std::vector<std::vector<Operation>> jobs;
    
    std::vector<std::vector<int>> readMatrix(const std::string& filename);
};

#endif // JSSP_PROBLEM_H