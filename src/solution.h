#ifndef JSSP_SOLUTION_H
#define JSSP_SOLUTION_H

#include "jssp_problem.h"

#include <vector>
#include <unordered_map>
#include <random>
#include <map>
#include <queue>
#include <functional>

class JSSPSolution {
public:    
    JSSPSolution();
    JSSPSolution(const JSSPProblem& problem);
    JSSPSolution(const JSSPSolution& other);
    void generateNeighbor(std::mt19937_64& rng);
    const std::vector<std::vector<Operation>>& getMachineJobOrders() const { return machineJobOrders; }
    int getMakespan() const { return makespan; }
    const std::vector<int>& getPi() const { return pi; }

private:
    struct Vertex {
        Operation op;
        std::vector<std::pair<int, int>> outEdges;
    };

    struct pair_hash {
        std::size_t operator()(const std::pair<int, int>& p) const {
            return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
        }
    };

    const JSSPProblem* problem;
    std::vector<std::vector<Operation>> machineJobOrders;
    int makespan;
    std::vector<Vertex> graph;
    int sourceVertex, sinkVertex;
    std::map<std::pair<int, int>, int> operationVertices;
    std::vector<int> pi;
    std::unordered_map<int, int> piDict;
    std::vector<int> longestPaths;
    std::vector<int> inTerminals;
    std::unordered_map<std::pair<int, int>, int, pair_hash> pq;
    std::unordered_map<int, std::unordered_map<int, int>> outEdges;
    std::unordered_map<int, std::unordered_map<int, int>> inEdges;

    void initGraph();
    void updateGraph();
    void updateMakespan();
    int calculateMakespan();
    bool hasCycle() const;
    void initialize();
    void buildPriorityQueue();
    void computeLongestPaths();
    void recomputeLongestPaths(int startVertex);
    void updateConnections(int vertex, int delta);
};

#endif