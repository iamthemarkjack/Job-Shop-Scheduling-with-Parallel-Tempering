#include "solution.h"

#include <algorithm>
#include <stdexcept>
#include <queue>
#include <deque>
#include <unordered_set>
#include <climits>
#include <iostream>
#include <chrono>

JSSPSolution::JSSPSolution() 
    : problem(nullptr), makespan(0), sourceVertex(-1), sinkVertex(-1) {
}

JSSPSolution::JSSPSolution(const JSSPProblem& problem) 
    : problem(&problem), makespan(0) {
    
    int numMachines = problem.getNumMachines();
    int numJobs = problem.getNumJobs();
    
    machineJobOrders.resize(numMachines);

    // Use greedy initialization
    initialize();
    
    // Compute longest paths and topological ordering
    computeLongestPaths();

    buildPriorityQueue();
    
    updateMakespan();
}

JSSPSolution::JSSPSolution(const JSSPSolution& other)
    : problem(other.problem), 
      machineJobOrders(other.machineJobOrders), 
      makespan(other.makespan),
      graph(other.graph), 
      sourceVertex(other.sourceVertex), 
      sinkVertex(other.sinkVertex), 
      operationVertices(other.operationVertices),
      pi(other.pi),
      piDict(other.piDict),
      longestPaths(other.longestPaths),
      inTerminals(other.inTerminals),
      outEdges(other.outEdges),
      inEdges(other.inEdges),
      pq(other.pq) {
}

void JSSPSolution::generateNeighbor(std::mt19937_64& rng) {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a JSSP problem");
    }

    bool found = false;
    std::pair<int,int> topPair;

    if (!pq.empty()){
        for (auto it = pq.begin(); it != pq.end(); ++it){
            if (it->second == 2){
                topPair = it->first;
                pq.erase(it);
                found = true;
                break;
            }
        }
    }
    
    if (found) {                    
            int a = topPair.first;
            int b = topPair.second;
            int aIdx = piDict[a];
            int bIdx = piDict[b];
            int c = pi[aIdx + 1];
            
            int aJob = graph[a].op.job;
            int bJob = graph[b].op.job;
            int cJob = graph[c].op.job;
            
            if ((cJob != aJob && cJob != bJob) || (cJob == aJob)) {
                // (a, c, b) to (b, a, c)
                pi[aIdx] = b;
                pi[aIdx + 1] = a;
                pi[bIdx] = c;
                piDict[a] = aIdx + 1;
                piDict[b] = aIdx;
                piDict[c] = bIdx;
                
                // Update priority queue and edge details
                pq[{b, a}] = 1;
                outEdges[b][a] = 1;
                inEdges[a][b] = 1;
                
                // Update connections
                updateConnections(a, 1);
                updateConnections(c, 1);
                updateConnections(b, -2);
            }
            else {
                // (a, c, b) to (c, b, a)
                pi[aIdx] = c;
                pi[aIdx + 1] = b;
                pi[bIdx] = a;
                piDict[a] = bIdx;
                piDict[b] = aIdx + 1;
                piDict[c] = aIdx;
                
                // Update priority queue and graph
                pq[{b, a}] = 1;
                outEdges[b][a] = 1;
                inEdges[a][b] = 1;
                
                // Update connections
                updateConnections(a, 2);
                updateConnections(c, -1);
                updateConnections(b, -1);
            }

            // Remove b from a's outgoing edges
            auto& aOutEdges = graph[a].outEdges;
            aOutEdges.erase(
                std::remove_if(aOutEdges.begin(), aOutEdges.end(),
                    [b](const std::pair<int, int>& edge) {
                        return edge.first == b;
                    }),
                aOutEdges.end()
            );

            // Add a to b's outgoing edges
            int aDuration = graph[a].op.duration;
            graph[b].outEdges.push_back({a, aDuration});

            // Update the outEdges and inEdges maps to match the graph changes
            if (outEdges[a].find(b) != outEdges[a].end()) {
                outEdges[a].erase(b);
            }
            if (inEdges[b].find(a) != inEdges[b].end()) {
                inEdges[b].erase(a);
            }
            
            recomputeLongestPaths(b);
        }
    
    if (!found) {        
        while (true) {            
            int numMachines = problem->getNumMachines();
            std::uniform_int_distribution<int> machineDist(0, numMachines - 1);
            int machine = machineDist(rng);
            
            if (machineJobOrders[machine].size() < 2) {
                continue;
            }
            
            std::uniform_int_distribution<int> posDist(0, machineJobOrders[machine].size() - 2);
            int pos = posDist(rng);
            
            // Swap operations
            std::swap(machineJobOrders[machine][pos], machineJobOrders[machine][pos + 1]);
            
            updateGraph();
            
            if (!hasCycle()) {
                // Recompute paths and priority queue
                computeLongestPaths();
                buildPriorityQueue();
                updateMakespan();
                break;
            }
            else {
                // Undo the swap if it creates a cycle
                std::swap(machineJobOrders[machine][pos], machineJobOrders[machine][pos + 1]);
            }
        }
    }
}

void JSSPSolution::initGraph() {
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
        
        // Track operations that go directly to sink
        inTerminals.push_back(lastOpVertex);
    }
    
    // Initialize longest paths vector
    longestPaths.resize(graph.size(), 0);
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

            int job1 = operation1.job;
            int job2 = operation2.job;
            int machine = operation1.machine; // Should be the same as operation2.machine
            
            int op1Vertex = operationVertices[{job1, operation1.jobIndex}];
            int op2Vertex = operationVertices[{job2, operation2.jobIndex}];

            graph[op1Vertex].outEdges.push_back({op2Vertex, operation1.duration});
        }
    }
}

void JSSPSolution::updateMakespan(){
    makespan = calculateMakespan();
}

int JSSPSolution::calculateMakespan(){
    return longestPaths[sinkVertex];
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

// void JSSPSolution::initialize() {
//     int numMachines = problem->getNumMachines();
//     int numJobs = problem->getNumJobs();
    
//     // Clear machine job orders
//     for (int m = 0; m < numMachines; ++m) {
//         machineJobOrders[m].clear();
//     }
    
//     // Track next operation index for each job
//     std::vector<int> jobProgress(numJobs, 0);
    
//     // Track when each machine and job becomes available
//     std::vector<int> machineTime(numMachines, 0);
//     std::vector<int> jobTime(numJobs, 0);
    
//     // Keep processing until all operations are scheduled
//     bool allDone = false;
//     while (!allDone) {
//         allDone = true;
        
//         // Find machine with earliest availability
//         int earliestMachine = 0;
//         int earliestTime = INT_MAX;
        
//         for (int m = 0; m < numMachines; ++m) {
//             if (machineTime[m] < earliestTime) {
//                 earliestMachine = m;
//                 earliestTime = machineTime[m];
//             }
//         }
        
//         // Find a job that needs this machine
//         bool foundJob = false;
//         int selectedJob = -1;
//         int earliestJobTime = INT_MAX;
        
//         for (int j = 0; j < numJobs; ++j) {
//             // Skip completed jobs
//             if (jobProgress[j] >= numMachines) {
//                 continue;
//             }
            
//             allDone = false; // At least one job isn't done
            
//             // Check if this job's next operation needs the earliest machine
//             Operation op = problem->getJobs()[j][jobProgress[j]];
//             if (op.machine == earliestMachine && jobTime[j] < earliestJobTime) {
//                 selectedJob = j;
//                 earliestJobTime = jobTime[j];
//                 foundJob = true;
//             }
//         }
        
//         if (foundJob) {
//             // Schedule this operation
//             Operation op = problem->getJobs()[selectedJob][jobProgress[selectedJob]];
//             machineJobOrders[earliestMachine].push_back(op);
            
//             // Update times
//             int startTime = std::max(machineTime[earliestMachine], jobTime[selectedJob]);
//             machineTime[earliestMachine] = startTime + op.duration;
//             jobTime[selectedJob] = startTime + op.duration;
            
//             // Move to next operation for this job
//             jobProgress[selectedJob]++;
//         } else if (!allDone) {
//             // No job needs this machine right now, advance the machine time
//             machineTime[earliestMachine] = INT_MAX; // Mark as temporarily unavailable
//         }
//     }
    
//     initGraph();
//     updateGraph();
// }

void JSSPSolution::initialize() {
    int numMachines = problem->getNumMachines();
    int numJobs = problem->getNumJobs();
    
    // Clear machine job orders
    for (int m = 0; m < numMachines; ++m) {
        machineJobOrders[m].clear();
    }
    
    // Create a list of operations with earliest start times
    std::vector<std::pair<int, int>> readyJobs; // (job, jobIndex)
    std::vector<int> jobProgress(numJobs, 0);
    
    // Initialize with first operation of each job
    for (int j = 0; j < numJobs; ++j) {
        readyJobs.push_back({j, 0});
    }
    
    // Track which machines are occupied and when they'll be free
    std::vector<int> machineReleaseTime(numMachines, 0);
    // Track when each job's last operation completes
    std::vector<int> jobReleaseTime(numJobs, 0);
    
    // Use priority dispatching to schedule all operations
    while (!readyJobs.empty()) {
        // Find operation with earliest possible start time
        int bestJob = -1;
        int bestJobIndex = -1;
        int earliestStart = INT_MAX;
        
        for (const auto& jp : readyJobs) {
            int j = jp.first;
            int idx = jp.second;
            
            Operation op = problem->getJobs()[j][idx];
            int machine = op.machine;
            
            // Operation can start when both job and machine are available
            int startTime = std::max(jobReleaseTime[j], machineReleaseTime[machine]);
            
            if (startTime < earliestStart) {
                earliestStart = startTime;
                bestJob = j;
                bestJobIndex = idx;
            }
        }
        
        // Schedule the selected operation
        Operation op = problem->getJobs()[bestJob][bestJobIndex];
        int machine = op.machine;
        
        // Add operation to machine's job order
        machineJobOrders[machine].push_back(op);
        
        // Update release times
        machineReleaseTime[machine] = earliestStart + op.duration;
        jobReleaseTime[bestJob] = earliestStart + op.duration;
        
        // Remove scheduled operation from ready list
        readyJobs.erase(
            std::remove_if(readyJobs.begin(), readyJobs.end(),
                [bestJob, bestJobIndex](const std::pair<int, int>& jp) {
                    return jp.first == bestJob && jp.second == bestJobIndex;
                }),
            readyJobs.end()
        );
        
        // Add next operation of the job if available
        if (bestJobIndex + 1 < numMachines) {
            readyJobs.push_back({bestJob, bestJobIndex + 1});
        }
    }
    
    initGraph();
    updateGraph();
}

void JSSPSolution::computeLongestPaths() {
    // Ensure graph is correctly built before computation
    if (graph.empty()) {
        throw std::runtime_error("Graph structure not initialized");
    }

    // Initialize data structures for topological sort
    std::vector<int> inDegree(graph.size(), 0);
    
    // Calculate in-degree for each vertex
    for (size_t u = 0; u < graph.size(); ++u) {
        for (const auto& edge : graph[u].outEdges) {
            int v = edge.first;
            inDegree[v]++;
        }
    }
    
    // Initialize longest paths to 0
    longestPaths.resize(graph.size(), 0);
    std::fill(longestPaths.begin(), longestPaths.end(), 0);
    
    // Reset topological ordering
    pi.clear();
    piDict.clear();
    
    // Queue for topological sort (start with vertices having no incoming edges)
    std::queue<int> q;
    for (size_t i = 0; i < inDegree.size(); ++i) {
        if (inDegree[i] == 0) {
            q.push(i);
        }
    }
    
    // Process vertices in topological order
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        // Add non-terminal vertices to pi for the topological ordering
        if (u != sourceVertex && u != sinkVertex) {
            pi.push_back(u);
            piDict[u] = pi.size() - 1;
        }
        
        // Process outgoing edges to update longest paths
        for (const auto& edge : graph[u].outEdges) {
            int v = edge.first;
            int weight = edge.second;
            
            // Update longest path to v if a longer path is found through u
            longestPaths[v] = std::max(longestPaths[v], longestPaths[u] + weight);
            
            // Decrement in-degree and enqueue if no incoming edges remain
            inDegree[v]--;
            if (inDegree[v] == 0) {
                q.push(v);
            }
        }
    }
}

void JSSPSolution::buildPriorityQueue() {
    // Clear previous data
    pq.clear();
    outEdges.clear();
    inEdges.clear();
    
    // Build priority queue based on machine job orders
    for (int m = 0; m < problem->getNumMachines(); ++m) {
        for (size_t i = 0; i < machineJobOrders[m].size() - 1; ++i) {
            Operation op1 = machineJobOrders[m][i];
            Operation op2 = machineJobOrders[m][i + 1];
            
            int job1 = op1.job;
            int job2 = op2.job;
            int idx1 = op1.jobIndex;
            int idx2 = op2.jobIndex;
            
            int u = operationVertices[{job1, idx1}];
            int v = operationVertices[{job2, idx2}];
            
            int uIdx = piDict[u];
            int vIdx = piDict[v];
            int distance = vIdx - uIdx;
            
            if (distance > 0) {
                pq[{u, v}] = distance;
                outEdges[u][v] = distance;
                inEdges[v][u] = distance;
            } else {
                throw std::runtime_error("Invalid topological ordering");
            }
        }
    }
}

void JSSPSolution::recomputeLongestPaths(int bVertex) {
    int startIdx = piDict[bVertex];
    
    // Recompute longest paths from b Vertex onwards
    for (size_t i = startIdx; i < pi.size(); ++i) {
        int v = pi[i];
        longestPaths[v] = 0; // Reset path length
        
        // Check all incoming edges from the graph structure
        for (int u = 0; u < graph.size(); ++u) {
            for (const auto& edge : graph[u].outEdges) {
                if (edge.first == v) {
                    int weight = edge.second;
                    longestPaths[v] = std::max(longestPaths[v], longestPaths[u] + weight);
                }
            }
        }
    }
    
    // Update makespan based on terminal nodes
    longestPaths[sinkVertex] = 0;
    for (int terminal : inTerminals) {
        int pathLength = 0;
        
        // Find the edge from terminal to sink
        for (const auto& edge : graph[terminal].outEdges) {
            if (edge.first == sinkVertex) {
                pathLength = longestPaths[terminal] + edge.second;
                break;
            }
        }
        
        longestPaths[sinkVertex] = std::max(longestPaths[sinkVertex], pathLength);
    }

    // Update the makespan
    makespan = longestPaths[sinkVertex];
}

void JSSPSolution::updateConnections(int vertex, int delta) {
    // Update incoming edges
    for (auto& pair : inEdges[vertex]) {
        int u = pair.first;
        int newDistance = pair.second + delta;
        
        pq[{u, vertex}] = newDistance;

        outEdges[u][vertex] += delta;
        inEdges[vertex][u] += delta;
    }
    
    // Update outgoing edges
    for (auto& pair : outEdges[vertex]) {
        int v = pair.first;
        int newDistance = pair.second - delta;
        
        pq[{vertex, v}] = newDistance;

        outEdges[vertex][v] -= delta;
        inEdges[v][vertex] -= delta;
    }
}