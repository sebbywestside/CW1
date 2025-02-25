#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <cmath>

// Include custom headers for debugging and graphics displaying 
#include "debug.h"
#include "terminal_graphics.h"

// Struct for eddy current parameters
// stores amplitude and desired gradient 
struct EddyCurrentParams {
    float amplitude;
    float rateConstant;
};

// Global variables
std::vector<EddyCurrentParams> eddy_components;
std::vector<float> desired_gradient;

// Load eddy current parameters from file
void loadParameters(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Failed to open parameters file: " + filename);
    }
    float amplitude, rate;
    while (infile >> amplitude >> rate) {
        eddy_components.push_back({amplitude, rate});
    }
    if (eddy_components.empty()) {
        throw std::runtime_error("No eddy current parameters found in: " + filename);
    }
}

// Load gradient time course from file
void loadGradient(const std::string& filename, std::vector<float>& gradient) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Failed to open gradient file: " + filename);
    }
    float value;
    gradient.clear();
    while (infile >> value) {
        gradient.push_back(value);
    }
    if (gradient.empty()) {
        throw std::runtime_error("No gradient data found in: " + filename);
    }
}

// Compute the predicted gradient time course based on eddy currents
std::vector<float> computePredictedGradient(const std::vector<float>& input_gradient) {
    std::vector<float> predicted_gradient = input_gradient;
    std::vector<std::vector<float>> currents(eddy_components.size(), std::vector<float>(input_gradient.size(), 0.0f));

    for (size_t t = 1; t < input_gradient.size(); ++t) {
        float dG = input_gradient[t] - input_gradient[t - 1];
        for (size_t n = 0; n < eddy_components.size(); ++n) {
            currents[n][t] = currents[n][t - 1] + dG - eddy_components[n].rateConstant * currents[n][t - 1];
            predicted_gradient[t] -= eddy_components[n].amplitude * currents[n][t];
        }
    }
    return predicted_gradient;
}

// Compute maximum absolute deviation between two vectors
float computeMaxDeviation(const std::vector<float>& desired, const std::vector<float>& predicted) {
    if (desired.size() != predicted.size()) {
        throw std::runtime_error("Gradient vectors size mismatch");
    }
    float max_dev = 0.0f;
    for (size_t i = 0; i < desired.size(); ++i) { 
        float dev = std::abs(desired[i] - predicted[i]);
        max_dev = std::max(max_dev, dev);
    }
    return max_dev;
}

// Display gradients using terminal graphics with dynamic ylim based on all data
void displayGradients(const std::vector<float>& input, const std::vector<float>& predicted) {
    if (input.empty() || predicted.empty()) {
        std::cerr << "Gradient data is empty. Nothing to display.\n";
        return;
    }

    // Find the min and max across desired, input, and predicted gradients
    std::vector<float> all_values;
    all_values.insert(all_values.end(), desired_gradient.begin(), desired_gradient.end());
    all_values.insert(all_values.end(), input.begin(), input.end());
    all_values.insert(all_values.end(), predicted.begin(), predicted.end());

    auto [min_it, max_it] = std::minmax_element(all_values.begin(), all_values.end());
    float min_val = *min_it;
    float max_val = *max_it;

    // Add a small buffer (e.g., 10% of the range) to ensure all data is visible
    float range = max_val - min_val;
    float buffer = range * 0.2;
    float new_min = min_val - buffer;
    float new_max = max_val + buffer;

    TG::plot(1600, 600)
        .set_ylim(new_min, new_max)  // Dynamic y-axis limits with buffer
        .add_line(input, 2)             // yellow
        .add_line(predicted, 3);         // magenta
}

// Save the compensated gradient to a file
void saveGradient(const std::string& filename, const std::vector<float>& gradient) {
    std::ofstream outfile(filename);
    if (!outfile) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
    for (float value : gradient) {
        outfile << std::scientific << value << "\n";
    }
}

// Main program logic with 0-based iteration counter, showing only initial and last iterations graphically
void run(std::vector<std::string>& args) {
    // Handle verbose flag
    debug::verbose = std::erase(args, "-v");

    // Parse number of iterations
    int num_iterations = 10; // Default
    auto n_it = std::find(args.begin(), args.end(), "-n");
    if (n_it != args.end() && std::next(n_it) != args.end()) {
        num_iterations = std::stoi(*std::next(n_it));
        args.erase(n_it, n_it + 2);
    }

    // Check minimum arguments
    if (args.size() < 3) {
        throw std::runtime_error("Missing arguments - expected at least 2 arguments: <program> <params_file> <gradient_file> [output_file]");
    }

    // Load input data
    std::cerr << "Reading parameters file \"" << args[1] << "\"\n";
    loadParameters(args[1]);
    std::cerr << "Reading gradient file \"" << args[2] << "\"\n";
    loadGradient(args[2], desired_gradient);

    // Initialize input gradient as a copy of desired gradient
    std::vector<float> input_gradient = desired_gradient;

    // Compute and display initial state 
    std::vector<float> predicted_gradient = computePredictedGradient(input_gradient);
    float max_dev = computeMaxDeviation(desired_gradient, predicted_gradient);
    std::cerr << "iteration " << 0 << ", maximum absolute deviation = " << max_dev << "\n";
    displayGradients(input_gradient, predicted_gradient); // Show initial plot

    // Iterative pre-emphasis with 0-based iteration counter
    for (int iter = 0; iter < num_iterations - 1; ++iter) { 
        // Compute difference and update input gradient
        for (size_t i = 0; i < desired_gradient.size(); ++i) {
            float diff = desired_gradient[i] - predicted_gradient[i];
            input_gradient[i] += diff;
        }
        // Compute predicted gradient 
        predicted_gradient = computePredictedGradient(input_gradient);
    }

    // Display and report the last iteration (e.g., iteration 10 for -n 10)
    // Compute difference and update input gradient for the last iteration
    for (size_t i = 0; i < desired_gradient.size(); ++i) {
        float diff = desired_gradient[i] - predicted_gradient[i];
        input_gradient[i] += diff;
    }

    predicted_gradient = computePredictedGradient(input_gradient);
    max_dev = computeMaxDeviation(desired_gradient, predicted_gradient);
    std::cerr << "iteration " << num_iterations << ", maximum absolute deviation = " << max_dev << "\n";
    displayGradients(input_gradient, predicted_gradient); // Show final plot

    // Save output if requested
    if (args.size() > 3) {
        std::cerr << "Writing compensated gradient to \"" << args[3] << "\"\n";
        saveGradient(args[3], input_gradient);
    }
}

int main(int argc, char* argv[]) {
    try {
        std::vector<std::string> args(argv, argv + argc);
        run(args);
    }
    catch (std::exception& excp) {
        std::cerr << "ERROR: " << excp.what() << " - aborting\n";
        return 1;
    }
    catch (...) {
        std::cerr << "ERROR: unknown exception thrown - aborting\n";
        return 1;
    }
    return 0;
}