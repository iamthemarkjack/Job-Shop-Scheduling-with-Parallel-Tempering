#include "tsp_problem.h"
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <iostream>

TSPProblem::TSPProblem(int numCities) : numCities(numCities) {
    if (numCities > 0) {
        coordinates.resize(numCities, {0.0, 0.0});
        initializeDistanceMatrix();
    }
}

TSPProblem::TSPProblem(const std::vector<std::pair<double, double>>& coordinates) 
    : numCities(coordinates.size()), coordinates(coordinates) {
    initializeDistanceMatrix();
}

bool TSPProblem::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    // Clear current data
    coordinates.clear();
    
    // Simple parser for TSPLIB format
    std::string line;
    bool readingCoordinates = false;
    
    while (std::getline(file, line)) {
        if (line.find("NODE_COORD_SECTION") != std::string::npos) {
            readingCoordinates = true;
            continue;
        }
        
        if (line.find("EOF") != std::string::npos) {
            break;
        }
        
        if (readingCoordinates) {
            std::istringstream iss(line);
            int id;
            double x, y;
            
            if (iss >> id >> x >> y) {
                coordinates.push_back({x, y});
            }
        }
    }
    
    numCities = coordinates.size();
    
    if (numCities == 0) {
        std::cerr << "No valid city coordinates found in file" << std::endl;
        return false;
    }
    
    // Initialize distance matrix
    initializeDistanceMatrix();
    return true;
}

void TSPProblem::generateRandom(int numCities, double minCoord, double maxCoord) {
    if (numCities <= 0) {
        throw std::invalid_argument("Number of cities must be positive");
    }
    
    this->numCities = numCities;
    coordinates.clear();
    
    // Generate random coordinates
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minCoord, maxCoord);
    
    for (int i = 0; i < numCities; ++i) {
        coordinates.push_back({dist(gen), dist(gen)});
    }
    
    // Initialize distance matrix
    initializeDistanceMatrix();
}

double TSPProblem::getDistance(int city1, int city2) const {
    if (city1 < 0 || city1 >= numCities || city2 < 0 || city2 >= numCities) {
        throw std::out_of_range("City index out of range");
    }
    
    return distanceMatrix[city1][city2];
}

double TSPProblem::calculateTourLength(const std::vector<int>& tour) const {
    if (tour.size() != numCities) {
        throw std::invalid_argument("Tour must visit all cities exactly once");
    }
    
    double length = 0.0;
    
    for (size_t i = 0; i < tour.size(); ++i) {
        int current = tour[i];
        int next = tour[(i + 1) % tour.size()];
        
        length += getDistance(current, next);
    }
    
    return length;
}

void TSPProblem::initializeDistanceMatrix() {
    distanceMatrix.resize(numCities, std::vector<double>(numCities, 0.0));
    
    for (int i = 0; i < numCities; ++i) {
        for (int j = i + 1; j < numCities; ++j) {
            double dist = calculateEuclideanDistance(i, j);
            distanceMatrix[i][j] = dist;
            distanceMatrix[j][i] = dist;  // Symmetric matrix
        }
    }
}

double TSPProblem::calculateEuclideanDistance(int city1, int city2) const {
    const auto& coord1 = coordinates[city1];
    const auto& coord2 = coordinates[city2];
    
    double dx = coord1.first - coord2.first;
    double dy = coord1.second - coord2.second;
    
    return std::sqrt(dx * dx + dy * dy);
}