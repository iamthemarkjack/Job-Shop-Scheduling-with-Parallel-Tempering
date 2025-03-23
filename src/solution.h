#ifndef JSSP_SOLUTION_H
#define JSSP_SOLUTION_H

#include "jssp_problem.h"
#include <vector>
#include <random>
#include <memory>
#include <map>
#include <boost/graph/adjacency_list.hpp>

class JSSPSolution {
public:
    JSSPSolution();
    JSSPSolution(const JSSPProblem& problem);
    ~JSSPSolution();
    
    JSSPSolution(const JSSPSolution& other);
    JSSPSolution& operator=(const JSSPSolution& other);
    
    void randomize(std::mt19937_64& rng);
    
    JSSPSolution generateNeighbor(std::mt19937_64& rng) const;
    
    int getMakespan() const { return makespan; }
    
    const std::vector<std::vector<int>>& getMachineJobOrder() const { return machineJobOrder; }
        
private:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, 
                                  boost::property<boost::vertex_name_t, Operation>,
                                  boost::property<boost::edge_weight_t, int>> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;
    
    const JSSPProblem* problem;
    std::vector<std::vector<int>> machineJobOrder; // For each machine, the order of jobs
    int makespan;
    
    // Graph representation
    Graph* precedenceGraph;
    Vertex sourceVertex;
    Vertex sinkVertex;
    std::map<std::pair<int, int>, Vertex> operationVertices; // Maps (job, index) to vertex
    
    void buildGraph();
    void updateGraph();
    void updateMakespan();
    int calculateMakespan();
    
    bool hasCycle() const;
};

#endif // JSSP_SOLUTION_H