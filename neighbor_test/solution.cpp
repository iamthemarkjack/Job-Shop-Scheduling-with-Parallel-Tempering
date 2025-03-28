#include "solution.h"
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include <queue>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>

struct cycle_detector : public boost::dfs_visitor<> {
    cycle_detector(bool& has_cycle) : has_cycle(has_cycle) {}
    
    template <class Edge, class Graph>
    void back_edge(Edge, const Graph&) {
        has_cycle = true;
    }
    
private:
    bool& has_cycle;
};

JSSPSolution::JSSPSolution() 
    : problem(nullptr), makespan(0), precedenceGraph(nullptr) {
}

JSSPSolution::JSSPSolution(const JSSPProblem& problem) 
    : problem(&problem), makespan(0), precedenceGraph(new Graph()) {
    
    int numMachines = problem.getNumMachines();
    int numJobs = problem.getNumJobs();
    
    machineJobOrder.resize(numMachines);
    
    for (int m = 0; m < numMachines; ++m) {
        std::vector<int> jobsOnMachine;
        for (int j = 0; j < numJobs; ++j) {
            bool usesThisMachine = false;
            for (const auto& op : problem.getJobs()[j]) {
                if (op.machine == m) {
                    usesThisMachine = true;
                    break;
                }
            }
            if (usesThisMachine) {
                jobsOnMachine.push_back(j);
            }
        }
        machineJobOrder[m] = jobsOnMachine;
    }
    
    buildGraph();
    updateMakespan();
}

JSSPSolution::~JSSPSolution() {
    delete precedenceGraph;
}

JSSPSolution::JSSPSolution(const JSSPSolution& other)
    : problem(other.problem), 
      machineJobOrder(other.machineJobOrder), 
      makespan(other.makespan),
      precedenceGraph(nullptr) {
    
    if (problem) {
        precedenceGraph = new Graph();
        buildGraph();
    }
}

JSSPSolution& JSSPSolution::operator=(const JSSPSolution& other) {
    if (this != &other) {
        delete precedenceGraph;
        
        problem = other.problem;
        machineJobOrder = other.machineJobOrder;
        makespan = other.makespan;
        
        precedenceGraph = nullptr;
        
        if (problem) {
            precedenceGraph = new Graph();
            buildGraph();
        }
    }
    return *this;
}

void JSSPSolution::randomize(std::mt19937_64& rng) {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a JSSP problem");
    }
    
    int numMachines = problem->getNumMachines();
    
    for (int m = 0; m < numMachines; ++m) {
        std::shuffle(machineJobOrder[m].begin(), machineJobOrder[m].end(), rng);
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

    const int MAX_ATTEMPTS = 1000;
    int attempts = 0;
    
    bool validNeighborFound = false;
    
    while (!validNeighborFound && attempts < MAX_ATTEMPTS) {        
        attempts++;
        neighbor = *this;
        
        int machine = machineDist(rng);
        
        std::uniform_int_distribution<int> posDist(0, neighbor.machineJobOrder[machine].size() - 2);
        int pos = posDist(rng);
        
        std::swap(neighbor.machineJobOrder[machine][pos], neighbor.machineJobOrder[machine][pos + 1]);
        
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
        throw std::runtime_error("Solution not associated with a JSSP problem");
    }
    
    if (precedenceGraph) {
        delete precedenceGraph;
    }
    
    precedenceGraph = new Graph();
    operationVertices.clear();
    
    sourceVertex = boost::add_vertex(*precedenceGraph);
    sinkVertex = boost::add_vertex(*precedenceGraph);
    
    for (int i = 0; i < problem->getNumJobs(); ++i) {
        for (const auto& op : problem->getJobs()[i]) {
            Vertex v = boost::add_vertex(*precedenceGraph);
            boost::put(boost::vertex_name, *precedenceGraph, v, op);
            operationVertices[{op.job, op.index}] = v;
        }
    }
    
    for (int i = 0; i < problem->getNumJobs(); ++i) {
        Vertex firstOp = operationVertices[{i, 0}];
        Edge e = boost::add_edge(sourceVertex, firstOp, *precedenceGraph).first;
        boost::put(boost::edge_weight, *precedenceGraph, e, 0);
        
        for (int j = 0; j < problem->getNumMachines() - 1; ++j) {
            Vertex u = operationVertices[{i, j}];
            Vertex v = operationVertices[{i, j + 1}];
            Operation uOp = boost::get(boost::vertex_name, *precedenceGraph, u);
            
            Edge e = boost::add_edge(u, v, *precedenceGraph).first;
            boost::put(boost::edge_weight, *precedenceGraph, e, uOp.duration);
        }
        
        Vertex lastOp = operationVertices[{i, problem->getNumMachines() - 1}];
        Operation lastOpData = boost::get(boost::vertex_name, *precedenceGraph, lastOp);
        e = boost::add_edge(lastOp, sinkVertex, *precedenceGraph).first;
        boost::put(boost::edge_weight, *precedenceGraph, e, lastOpData.duration);
    }
    
    updateGraph();
}

void JSSPSolution::updateGraph() {
    if (!precedenceGraph || !problem) {
        return;
    }
    
    std::vector<Edge> edgesToRemove;
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(*precedenceGraph); ei != ei_end; ++ei) {
        Vertex u = boost::source(*ei, *precedenceGraph);
        Vertex v = boost::target(*ei, *precedenceGraph);
        
        if (u == sourceVertex || v == sinkVertex) {
            continue;
        }
        
        if (u != sourceVertex && v != sinkVertex) {
            Operation uOp = boost::get(boost::vertex_name, *precedenceGraph, u);
            Operation vOp = boost::get(boost::vertex_name, *precedenceGraph, v);
            
            if (uOp.job != vOp.job) {
                edgesToRemove.push_back(*ei);
            }
        }
    }
    
    for (const auto& edge : edgesToRemove) {
        boost::remove_edge(edge, *precedenceGraph);
    }
    
    for (int m = 0; m < problem->getNumMachines(); ++m) {
        const auto& jobOrder = machineJobOrder[m];
        
        for (size_t i = 0; i < jobOrder.size() - 1; ++i) {
            int job1 = jobOrder[i];
            int job2 = jobOrder[i + 1];
            
            Vertex op1Vertex;
            Vertex op2Vertex;
            bool foundOp1 = false;
            bool foundOp2 = false;
            
            for (int opIdx = 0; opIdx < problem->getJobs()[job1].size(); ++opIdx) {
                if (problem->getJobs()[job1][opIdx].machine == m) {
                    op1Vertex = operationVertices[{job1, opIdx}];
                    foundOp1 = true;
                    break;
                }
            }
            
            for (int opIdx = 0; opIdx < problem->getJobs()[job2].size(); ++opIdx) {
                if (problem->getJobs()[job2][opIdx].machine == m) {
                    op2Vertex = operationVertices[{job2, opIdx}];
                    foundOp2 = true;
                    break;
                }
            }
            
            if (foundOp1 && foundOp2) {
                Operation op1 = boost::get(boost::vertex_name, *precedenceGraph, op1Vertex);
                Edge e = boost::add_edge(op1Vertex, op2Vertex, *precedenceGraph).first;
                boost::put(boost::edge_weight, *precedenceGraph, e, op1.duration);
            }
        }
    }
}

void JSSPSolution::updateMakespan() {
    makespan = calculateMakespan();
}

int JSSPSolution::calculateMakespan() {
    if (!precedenceGraph || !problem) {
        throw std::runtime_error("Solution not properly initialized");
    }
    
    if (hasCycle()) {
        return std::numeric_limits<int>::max();
    }
    
    size_t numVertices = boost::num_vertices(*precedenceGraph);
    
    std::vector<int> earliestStart(numVertices, 0);
    std::vector<int> inDegree(numVertices, 0);
    
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(*precedenceGraph); ei != ei_end; ++ei) {
        Vertex v = boost::target(*ei, *precedenceGraph);
        inDegree[v]++;
    }
    
    std::queue<Vertex> q;
    
    q.push(sourceVertex);
    
    while (!q.empty()) {
        Vertex u = q.front();
        q.pop();
        
        boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(u, *precedenceGraph); ei != ei_end; ++ei) {
            Vertex v = boost::target(*ei, *precedenceGraph);
            int weight = boost::get(boost::edge_weight, *precedenceGraph, *ei);
            
            earliestStart[v] = std::max(earliestStart[v], earliestStart[u] + weight);
            inDegree[v]--;
            if (inDegree[v] == 0) {
                q.push(v);
            }
        }
    }
    
    return earliestStart[sinkVertex];
}

bool JSSPSolution::hasCycle() const {
    if (!precedenceGraph) {
        return false;
    }
    
    bool has_cycle = false;
    std::vector<boost::default_color_type> color_map(boost::num_vertices(*precedenceGraph));
    
    boost::depth_first_search(
        *precedenceGraph,
        boost::visitor(cycle_detector(has_cycle))
        .color_map(boost::make_iterator_property_map(
            color_map.begin(),
            boost::get(boost::vertex_index, *precedenceGraph)))
    );
    
    return has_cycle;
}