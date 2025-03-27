#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

class ParallelTempering {
public:
    ParallelTempering(int n_chains, int n_iter, double min_temp, double target_sap, double c, double eta);
    void run();
    void saveHistory() const;

private:
    int n_chains, n_iter;
    double min_temp, target_sap, c, eta;
    std::vector<double> chains, betas, rhos;
    std::vector<int> swap_attempted, swap_accepted;
    std::vector<double> ratio;
    std::vector<std::vector<double>> sap_history;
    std::vector<std::vector<double>> temp_history;
    std::mt19937 rng;
    
    double E(double x); // One dimensional Energy function
    double gamma(int t); // Decay function
    void metropolis_step(int chain_idx);
    void attempt_swap();
    void update_betas(int iteration);
};

ParallelTempering::ParallelTempering(int n_chains, int n_iter, double min_temp, double target_sap, double c, double eta)
    : n_chains(n_chains), n_iter(n_iter), min_temp(min_temp), target_sap(target_sap), c(c), eta(eta), rng(std::random_device{}()) {
    
    chains.resize(n_chains);
    betas.assign(n_chains, 1.0 / min_temp);
    rhos.assign(n_chains - 1, 1.0);
    swap_attempted.assign(n_chains - 1, 0);
    swap_accepted.assign(n_chains - 1, 0);
    ratio.assign(n_chains - 1, 0.0);
    sap_history.resize(n_chains - 1);
    temp_history.resize(n_chains);
    
    std::normal_distribution<double> norm(0.0, 1.0);
    for (int i = 0; i < n_chains; ++i) {
        chains[i] = norm(rng);
        if (i > 0) betas[i] = betas[i - 1] / 10.0;
    }
}

double ParallelTempering::E(double x) {
    return x * x * x * x - 3 * x * x + 5;
}

double ParallelTempering::gamma(int t) {
    return c * std::pow(t + 1, -eta);
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
    
    swap_attempted[l]++;
    if (std::uniform_real_distribution<>(0.0, 1.0)(rng) < sap) {
        std::swap(chains[l], chains[l + 1]);
        swap_accepted[l]++;
    }

    for (int i = 0; i < n_chains - 1; ++i){
        double val = swap_attempted[i] > 0 ? static_cast<double>(swap_accepted[i]) / swap_attempted[i] : 0.0;
        ratio[i] = val;
    }
}

void ParallelTempering::update_betas(int iteration) {
    for (int k = 0; k < n_chains - 1; ++k) {
        double E1 = E(chains[k]), E2 = E(chains[k + 1]);
        double sap = std::min(1.0, std::exp((E1 - E2) * (betas[k] - betas[k + 1])));
        rhos[k] += std::min(10.0, std::max(-10.0, gamma(iteration) * (target_sap - sap)));
        betas[k + 1] = betas[k] * std::exp( -rhos[k] );
    }
}

void ParallelTempering::run() {
    for (int iteration = 0; iteration < n_iter; ++iteration) {
        for (int i = 0; i < n_chains; ++i) {
            metropolis_step(i);
        }
        attempt_swap();
        update_betas(iteration);

        for (int i = 0; i < n_chains - 1; ++i) {
            sap_history[i].push_back(ratio[i]);
        }

        for (int i = 0; i < n_chains; ++i) {
            temp_history[i].push_back(1.0 / betas[i]);
        }

        if (iteration % (n_iter / 10) == 0) {
            std::cout << "Iteration " << iteration << "/" << n_iter << "\n";

        }
    }
    saveHistory();
}

void ParallelTempering::saveHistory() const {
    std::ofstream file("sap_history.csv");
    file << "Iteration";
    for (int i = 0; i < n_chains - 1; ++i) {
        file << ",Pair" << i << "_" << i + 1;
    }
    file << "\n";
    
    for (size_t iter = 0; iter < sap_history[0].size(); ++iter) {
        file << iter;
        for (int i = 0; i < n_chains - 1; ++i) {
            file << "," << sap_history[i][iter];
        }
        file << "\n";
    }
    file.close();

    std::ofstream temp_file("temp_history.csv");
    temp_file << "Iteration";
    for (int i = 0; i < n_chains; ++i) {
        temp_file << ",Chain" << i;
    }
    temp_file << "\n";
    
    for (size_t iter = 0; iter < temp_history[0].size(); ++iter) {
        temp_file << iter;
        for (int i = 0; i < n_chains; ++i) {
            temp_file << "," << temp_history[i][iter];
        }
        temp_file << "\n";
    }
    temp_file.close();
}

int main() {
    ParallelTempering pt(5, 100000, 0.01, 0.234, 1, 2.0 / 3.0);
    pt.run();
    return 0;
}