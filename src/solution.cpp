#include "solution.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include <unordered_set>

Solution::Solution() : problem(nullptr), cost(0.0) {
}

Solution::Solution(const TSPProblem& problem) : problem(&problem), cost(0.0) {
    int numCities = problem.getNumCities();
    tour.resize(numCities);
    
    // Initialize with identity permutation
    for (int i = 0; i < numCities; ++i) {
        tour[i] = i;
    }
    
    updateCost();
}

void Solution::randomize(std::mt19937& rng) {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a TSP problem");
    }
    
    // Initialize with sequential order
    int numCities = problem->getNumCities();
    tour.resize(numCities);
    for (int i = 0; i < numCities; ++i) {
        tour[i] = i;
    }
    
    // Shuffle the tour
    std::shuffle(tour.begin(), tour.end(), rng);
    
    // Update the cost
    updateCost();
}

void Solution::applyRandomMove(std::mt19937& rng) {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a TSP problem");
    }
    
    int numCities = problem->getNumCities();
    
    // Choose a random move type: 0 = 2-opt, 1 = insertion, 2 = swap
    std::uniform_int_distribution<int> moveTypeDist(0, 2);
    int moveType = moveTypeDist(rng);
    
    std::uniform_int_distribution<int> cityDist(0, numCities - 1);
    
    switch (moveType) {
        case 0: {  // 2-opt
            int i = cityDist(rng);
            int j = cityDist(rng);
            while (i == j) {
                j = cityDist(rng);
            }
            if (i > j) std::swap(i, j);
            twoOpt(i, j);
            break;
        }
        case 1: {  // Insertion
            int removePos = cityDist(rng);
            int insertPos = cityDist(rng);
            while (removePos == insertPos) {
                insertPos = cityDist(rng);
            }
            insertMove(removePos, insertPos);
            break;
        }
        case 2: {  // Swap
            int i = cityDist(rng);
            int j = cityDist(rng);
            while (i == j) {
                j = cityDist(rng);
            }
            swapMove(i, j);
            break;
        }
    }
    
    // Update the tour cost
    updateCost();
}

void Solution::setTour(const std::vector<int>& newTour) {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a TSP problem");
    }
    
    if (newTour.size() != problem->getNumCities()) {
        throw std::invalid_argument("Tour size doesn't match problem size");
    }
    
    // Check if all cities are present exactly once
    std::unordered_set<int> cities;
    for (int city : newTour) {
        if (city < 0 || city >= problem->getNumCities()) {
            throw std::invalid_argument("Invalid city index in tour");
        }
        cities.insert(city);
    }
    
    if (cities.size() != problem->getNumCities()) {
        throw std::invalid_argument("Tour must visit each city exactly once");
    }
    
    tour = newTour;
    updateCost();
}

bool Solution::isValid() const {
    if (!problem) {
        return false;
    }
    
    if (tour.size() != problem->getNumCities()) {
        return false;
    }
    
    std::unordered_set<int> cities;
    for (int city : tour) {
        if (city < 0 || city >= problem->getNumCities()) {
            return false;
        }
        cities.insert(city);
    }
    
    return cities.size() == problem->getNumCities();
}

void Solution::twoOpt(int i, int j) {
    // Reverse the segment between positions i and j
    std::reverse(tour.begin() + i, tour.begin() + j + 1);
}

void Solution::insertMove(int removePos, int insertPos) {
    int city = tour[removePos];
    
    // Remove city from its current position
    tour.erase(tour.begin() + removePos);
    
    // Adjust insert position if it was after remove position
    if (insertPos > removePos) {
        insertPos--;
    }
    
    // Insert city at new position
    tour.insert(tour.begin() + insertPos, city);
}

void Solution::swapMove(int i, int j) {
    std::swap(tour[i], tour[j]);
}

void Solution::updateCost() {
    if (!problem) {
        throw std::runtime_error("Solution not associated with a TSP problem");
    }
    
    cost = problem->calculateTourLength(tour);
}