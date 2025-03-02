#ifndef SOLUTION_H
#define SOLUTION_H

#include "tsp_problem.h"
#include <vector>
#include <random>
#include <memory>

/**
 * Class representing a solution (tour) for a TSP problem.
 * Provides methods for manipulating and evaluating tours.
 */
class Solution {
public:
    // Constructors
    Solution();
    Solution(const TSPProblem& problem);
    
    // Initialize with a random tour
    void randomize(std::mt19937& rng);
    
    // Apply a random move operator to modify the solution
    void applyRandomMove(std::mt19937& rng);
    
    // Get the tour cost (total distance)
    double getCost() const { return cost; }
    
    // Get the tour
    const std::vector<int>& getTour() const { return tour; }
    
    // Set a specific tour
    void setTour(const std::vector<int>& newTour);
    
    // Check if the solution is valid
    bool isValid() const;

private:
    const TSPProblem* problem;
    std::vector<int> tour;
    double cost;
    
    // Different move operators
    void twoOpt(int i, int j);
    void insertMove(int removePos, int insertPos);
    void swapMove(int i, int j);
    
    // Recalculate the cost after changes
    void updateCost();
};

#endif // SOLUTION_H