#ifndef JSSP_SOLUTION_H
#define JSSP_SOLUTION_H

#include "jssp_problem.h"
#include <vector>
#include <random>
#include <memory>
#include <map>

class JSSPSolution {
public:
    JSSPSolution();
    JSSPSolution(const JSSPProblem& problem);

    JSSPSolution(const JSSPSolution& other); // Copy Constructor to create deep copy
    
    void randomize(std::mt19937_64& rng); // Randomly initialize a solution
    
    JSSPSolution generateNeighbor(std::mt19937_64& rng) const;

    const std::vector<std::vector<Operation>>& getMachineJobOrders() const { return machineJobOrders; }
    
    int getMakespan() const { return makespan; }
            
private:
    struct Vertex {
        Operation op;
        std::vector<std::pair<int, int>> outEdges; // (target vertex ID, weight)
    };

    const JSSPProblem* problem; // Pointer to the JSSP Problem instance
    std::vector<std::vector<Operation>> machineJobOrders; // For each machine, the order of jobs
    int makespan;

    // Graph representation
    std::vector<Vertex> graph;
    int sourceVertex;
    int sinkVertex;
    std::map<std::pair<int, int>, int> operationVertices; // Maps (job, index) to vertex ID
    
    void buildGraph();
    void updateGraph();
    void updateMakespan();
    int calculateMakespan();
    
    bool hasCycle() const;
    void topoSortUtil(int v, std::vector<bool>& visited, std::vector<int>& order) const;
    std::vector<int> topoSort() const;
};

#endif // JSSP_SOLUTION_H