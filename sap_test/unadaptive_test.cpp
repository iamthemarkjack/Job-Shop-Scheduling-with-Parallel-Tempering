#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>

class ParallelTempering {
public:
    ParallelTempering(int n_chains, int n_iter, double min_temp);
    void run();
    void saveHistory() const;

private:
    int n_chains, n_iter;
    double min_temp;
    std::vector<double> chains;
    std::vector<double> betas;
    std::vector<double> solution_proximity_history;
    std::mt19937 rng;

    constexpr static double TRUE_OPTIMUM = -1.92798;
    
    double E(double x); // One dimensional Energy function
    double gamma(int t); // Decay function
    void metropolis_step(int chain_idx);
    void attempt_swap();
    void update_betas(int iteration);
};

ParallelTempering::ParallelTempering(int n_chains, int n_iter, double min_temp)
    : n_chains(n_chains), n_iter(n_iter), min_temp(min_temp), rng(std::random_device{}()) {
    
    chains.resize(n_chains);
    betas.assign(n_chains, 1.0 / min_temp);
    
    std::normal_distribution<double> norm(100.0, 10.0);
    for (int i = 0; i < n_chains; ++i) {
        chains[i] = norm(rng);
        if (i > 0) betas[i] = betas[i - 1] / 10.0;
    }
}

double ParallelTempering::E(double x) {
    return std::sin(3*x) + std::cos(5*x) + (x*x) / 10;
}

void ParallelTempering::metropolis_step(int chain_idx) {
    std::normal_distribution<double> norm(0.0, 1.0);
    double curr_x = chains[chain_idx];
    double prop_x = curr_x + norm(rng);
    double dE = E(prop_x) - E(curr_x);
    if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < std::exp(-dE * betas[chain_idx])) {
        chains[chain_idx] = prop_x;
    }
}

void ParallelTempering::attempt_swap() {
    std::uniform_int_distribution<int> dist(0, n_chains - 2);
    int l = dist(rng);
    double E1 = E(chains[l]), E2 = E(chains[l + 1]);
    double b1 = betas[l], b2 = betas[l + 1];
    double sap = std::min(1.0, std::exp((E1 - E2) * (b1 - b2)));
    
    if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < sap) {
        std::swap(chains[l], chains[l + 1]);
    }
}

void ParallelTempering::run() {
    for (int iteration = 0; iteration < n_iter; ++iteration) {
        for (int i = 0; i < n_chains; ++i) {
            metropolis_step(i);
        }
        attempt_swap();

        double best_solution = *std::min_element(chains.begin(), chains.end(), 
        [this](double a, double b) { return E(a) < E(b); });

        double proximity = 1.0 / (1.0 + std::abs(E(best_solution) - TRUE_OPTIMUM));
        solution_proximity_history.push_back(proximity);

        if (iteration % (n_iter / 10) == 0) {
            std::cout << "Iteration " << iteration << "/" << n_iter << "\n";
        }
    }
    saveHistory();
}

void ParallelTempering::saveHistory() const {
    std::ofstream proximity_file("unadap_proximity.csv");
    proximity_file << "Iteration,Proximity" << "\n";

    for (size_t iter = 0; iter < solution_proximity_history.size(); ++iter){
        proximity_file << iter << "," << solution_proximity_history[iter] << "\n";
    }
    proximity_file.close();
}

int main() {
    ParallelTempering pt(4, 10000, 0.01);
    pt.run();
    return 0;
}