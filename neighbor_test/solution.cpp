#include "solution.h"

#include <algorithm>
#include <stdexcept>
#include <queue>

JSSPSolution::JSSPSolution() 
    : problem(nullptr), makespan(0), sourceVertex(-1), sinkVertex(-1) {
}

JSSPSolution::JSSPSolution(const JSSPProblem& problem) 
    : problem(&problem), makespan(0) {
    
    int numMachines = problem.getNumMachines();
    int numJobs = problem.getNumJobs();
    
    machineJobOrders.resize(numMachines);
    
    for (int m = 0; m < numMachines; ++m) {
        machineJobOrders[m] = problem.getMachineOperations(m);
    }
    
    buildGraph();
    updateMakespan();
}

JSSPSolution::JSSPSolution(const JSSPSolution& other)
    : problem(other.problem), 
      machineJobOrders(other.machineJobOrders), 
      makespan(other.makespan),
      graph(other.graph), 
      sourceVertex(other.sourceVertex), 
      sinkVertex(other.sinkVertex), 
      operationVertices(other.operationVertices) {
}

void JSSPSolution::randomize(std::mt19937_64& rng) {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a JSSP problem");
    }
    
    int numMachines = problem->getNumMachines();
    
    for (int m = 0; m < numMachines; ++m) {
        std::shuffle(machineJobOrders[m].begin(), machineJobOrders[m].end(), rng);
    }
    
    updateGraph();
    updateMakespan();
}

JSSPSolution JSSPSolution::generateNeighbor(std::mt19937_64& rng) const {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a JSSP problem");
    }
    
    JSSPSolution neighbor = *this;
    int numMachines = problem->getNumMachines();
    
    std::uniform_int_distribution<int> machineDist(0, numMachines - 1);

    const int MAX_ATTEMPTS = 100;
    int attempts = 0;
    
    bool validNeighborFound = false;
    
    while (!validNeighborFound && attempts < MAX_ATTEMPTS) {        
        attempts++;
        neighbor = *this;
        
        int machine = machineDist(rng);

        if (neighbor.machineJobOrders[machine].size() < 2) {
            continue;
        }
        
        std::uniform_int_distribution<int> posDist(0, neighbor.machineJobOrders[machine].size() - 2);
        int pos = posDist(rng);
        
        std::swap(neighbor.machineJobOrders[machine][pos], neighbor.machineJobOrders[machine][pos + 1]);
        
        neighbor.updateGraph();
        
        if (!neighbor.hasCycle()) {
            neighbor.updateMakespan();
            validNeighborFound = true;
        }
    }
    
    if (!validNeighborFound) {
        return *this;
    }
    
    return neighbor;
}

void JSSPSolution::buildGraph() {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a JSSP Problem");
    }

    graph.clear();
    operationVertices.clear();

    sourceVertex = 0;
    graph.push_back(Vertex()); // Source vertex

    int vertexCount = 1;

    // Add operation vertices
    for (int i = 0; i < problem->getNumJobs(); ++i){
        for (int j = 0; j < problem->getNumMachines(); ++j){
            Operation op = problem->getJobs()[i][j];
            graph.push_back(Vertex());
            graph.back().op = op;
            operationVertices[{i, j}] = vertexCount;
            vertexCount++;
        }
    }

    sinkVertex = vertexCount;
    graph.push_back(Vertex()); // Sink vertex

    // Add Conjunctive edges (source -> first op -> ... -> last op -> sink)
    for (int i = 0; i < problem->getNumJobs(); ++i){
        int firstOpVertex = operationVertices[{i, 0}];
        graph[sourceVertex].outEdges.push_back({firstOpVertex, 0}); // Weight 0 from source

        for (int j = 0; j < problem->getNumMachines() - 1; ++j){
            int currVertex = operationVertices[{i, j}];
            int nextVertex = operationVertices[{i, j + 1}];
            int duration = graph[currVertex].op.duration;

            graph[currVertex].outEdges.push_back({nextVertex, duration});
        }

        int lastOpVertex = operationVertices[{i, problem->getNumMachines() - 1}];
        int lastOpDuration = graph[lastOpVertex].op.duration;
        graph[lastOpVertex].outEdges.push_back({sinkVertex, lastOpDuration});
    }

    updateGraph();
}

void JSSPSolution::updateGraph(){
    if (problem == nullptr){
        return;
    }

    // Remove disjunctive edges keeping only conjunctive edges
    for (int v = 1; v < sinkVertex; ++v){
        auto& outEdges = graph[v].outEdges;
        auto it = outEdges.begin();
        while (it != outEdges.end()){
            int target = it->first;
            // Keep if it's sink or part of same job sequence
            if (target == sinkVertex || (target < sinkVertex && graph[v].op.job == graph[target].op.job)){
                ++it;
            } else {
                it = outEdges.erase(it);
            }
        }
    }

    // Adding the updated disjunctive edges
    for (int m = 0; m < problem->getNumMachines(); ++m){
        const auto& jobOrder = machineJobOrders[m];

        for (size_t i = 0; i < jobOrder.size() - 1; ++i){
            Operation operation1 = jobOrder[i];
            Operation operation2 = jobOrder[i + 1];

            int op1Vertex = operation1.job * problem->getNumMachines() + (operation1.machine + 1);
            int op2Vertex = operation2.job * problem->getNumMachines() + (operation2.machine + 1);

            graph[op1Vertex].outEdges.push_back({op2Vertex, operation1.duration});
        }
    }
}

void JSSPSolution::updateMakespan(){
    makespan = calculateMakespan();
}

int JSSPSolution::calculateMakespan(){
    if (problem == nullptr){
        throw std::runtime_error("Solution not properly initialized");
    }

    // Initialize longest path for all nodes
    std::vector<int> longestPath(graph.size(), 0);

    // Initialize in-degrees for all nodes
    std::vector<int> inDegree(graph.size(), 0);
    for (size_t u = 0; u < graph.size(); ++u){
        for (const auto& edge: graph[u].outEdges){
            int v = edge.first;
            inDegree[v]++;
        }
    }

    // Initialize the queue
    std::queue<int> q;
    for (size_t i = 0; i < inDegree.size(); ++i){
        if (inDegree[i] == 0){
            q.push(i);
        }
    }
    
    while (!q.empty()){
        int u = q.front();
        q.pop();
        
        for (const auto& edge: graph[u].outEdges){
            int v = edge.first;
            int weight = edge.second;

            longestPath[v] = std::max(longestPath[v], longestPath[u] + weight);
            inDegree[v]--;

            if (inDegree[v] == 0){
                q.push(v);
            }
        }
    }

    return longestPath[sinkVertex];
}

bool JSSPSolution::hasCycle() const {    
    std::vector<int> inDegree(graph.size(), 0);
    for (size_t u = 0; u < graph.size(); ++u) {
        for (const auto& edge : graph[u].outEdges) {
            int v = edge.first;
            inDegree[v]++;
        }
    }
    
    std::queue<int> q;
    for (size_t i = 0; i < inDegree.size(); ++i) {
        if (inDegree[i] == 0) {
            q.push(i);
        }
    }
    
    int count = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        count++;
        
        for (const auto& edge : graph[u].outEdges) {
            int v = edge.first;
            inDegree[v]--;
            if (inDegree[v] == 0) {
                q.push(v);
            }
        }
    }
    
    // If count is less than the number of vertices, there's a cycle
    return count < graph.size();
}