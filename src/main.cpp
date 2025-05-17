#include "parallel_tempering.h"
#include "jssp_problem.h"
#include "solution.h"

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <fstream>

void printHelp(){
    std::cout << "JSSP Solver using Automatic Parallel Tempering\n";
    std::cout << "Usage:\n";
    std::cout << "  --processing-times <filename>  Load processing times of JSSP instance from the file\n";
    std::cout << "  --machines <filename>          Load machines of JSSP instance from the file\n";
    std::cout << "  --replicas <num>               Number of temperature replicas (default: 8)\n";
    std::cout << "  --min-temp <value>             Minimum temperature (default: 1.0)\n";
    std::cout << "  --rho <value>                  Initial Values of rho (default: 1.0)\n";
    std::cout << "  --nsweep <num>                 MCMC sampling iterations before swap (default: 100)\n";
    std::cout << "  --time <value>                 Time limit for the solver in seconds (default: 60)\n";
    std::cout << "  --target-sap <value>           Target Swap Acceptance Probability (default: 0.234)\n";
    std::cout << "  --C <value>                    The scale of the polynomial decay (default: 1)\n";
    std::cout << "  --eta <value>                  The decay rate of the polynomial decay (default : 2/3)\n";
    std::cout << "  --seed <num>                   Random seed\n";
    std::cout << "  --his-prefix <string>          Prefix of history files (default: 'pt_run')\n";
    std::cout << "  --help                         Display this help message\n";
}

int main(int argc, char* argv[]){
    // Default parameters
    std::string processingTimesFile;
    std::string machinesFile;
    int numReplicas = 8;
    double minTemp = 1.0;
    double rho = 1.0;
    double maxTimeSeconds = 60.0;
    int N_sweep= 100;
    double target_SAP = 0.234;
    double C = 1;
    double eta = 2/3;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::string hisPrefix = "pt_run";

    // Parse command line arguments
    for (int i = 1; i < argc; ++i){
        std::string arg = argv[i];

        if (arg == "--help"){
            printHelp();
            return 0;
        } else if (arg == "--processing-times" && i + 1 < argc){
            processingTimesFile = argv[++i];
        } else if (arg == "--machines" && i + 1 < argc){
            machinesFile = argv[++i];
        } else if (arg == "--replicas" && i + 1 < argc){
            numReplicas = std::stoi(argv[++i]);
        } else if (arg == "--min-temp" && i + 1 < argc){
            minTemp = std::stod(argv[++i]);
        } else if (arg == "--rho" && i + 1 < argc){
            rho = std::stod(argv[++i]);
        } else if (arg == "--nsweep" && i + 1 < argc){
            N_sweep = std::stoi(argv[++i]);
        } else if (arg == "--time" && i + 1 < argc){
            maxTimeSeconds = std::stod(argv[++i]);
        } else if (arg == "--target-sap" && i + 1 < argc){
            target_SAP = std::stod(argv[++i]);
        } else if (arg == "--C" && i + 1 < argc){
            C = std::stod(argv[++i]);
        } else if (arg == "--eta" && i + 1 < argc){
            eta = std::stod(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc){
            seed = std::stoul(argv[++i]);
        } else if (arg == "--his-prefix" && i + 1 < argc){
            hisPrefix = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            printHelp();
            return 1;
        }
    }

    if (processingTimesFile.empty() || machinesFile.empty()){
        std::cerr << "Error: Must specify both processing times file and machines file\n";
        printHelp();
        return 1;
    }

    JSSPProblem problem;

    try {
        std::cout << "Loading JSSP instance from the processing times file: " << processingTimesFile << " and machines file: " << machinesFile << "\n";
        if (!problem.loadFromFiles(processingTimesFile, machinesFile)) {
            std::cerr << "Failed to load JSSP instance." << std::endl;
            return 1;
        }
        std::cout << "Loaded the JSSP instance.\n";

        std::cout << "Initializing parallel tempering solver with " << numReplicas << " replicas (threads).\n";
        std::cout << "Minimum Temperature: " << minTemp << "\n";

        AdaptiveParallelTempering solver(problem, numReplicas, minTemp, rho, maxTimeSeconds, N_sweep, target_SAP, C, eta, hisPrefix);

        std::cout << "Starting parallel tempering algorithm...\n";
        JSSPSolution bestSolution = solver.solve();

        std::cout << "Best tour length: " << std::fixed << std::setprecision(2) << bestSolution.getMakespan() << "\n";

    } catch (const std::exception& e){
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}