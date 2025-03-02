#ifndef TSP_PROBLEM_H
#define TSP_PROBLEM_H

#include <vector>
#include <string>
#include <cmath>

/*
* Class representing a Traveling Salesman Problem instance.
* Stores city coordinates and provices distance calculations.
 */

class TSPProblem{
public:
    // Constructor with number of cities
    TSPProblem(int numCities = 0);

    // Constructor from coordinates
    TSPProblem(const std::vector<std::pair<double, double>>& coorinates);

    // Load problem from a TSPLIB file
    bool loadFromFile(const std::string& filename);

    // Generate a random problem instance
    void generateRandom(int numCities, double minCoord = 0.0, double maxCoord = 100.0);

    // Getters
    int getNumCities() const { return numCities; }
    double getDistance(int city1, int city2) const;
    const std::vector<std::pair<double, double>>& getCoordinates() const {return coordinates; }

    // Calculate total tour length for a given permutation of cities
    double calculateTourLength(const std::vector<int>& tour) const;

private:
    int numCities;
    std::vector<std::pair<double, double>> coordinates;
    std::vector<std::vector<double>> distanceMatrix;

    // Initialize the distance matrix
    void initializeDistanceMatrix();

    // Calculate Euclidean distance between two cities
    double calculateEuclideanDistance(int city1, int city2) const;
};

#endif // TSP_PROBLEM_H