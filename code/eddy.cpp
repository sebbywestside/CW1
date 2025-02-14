#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "debug.h"
#include "terminal_graphics.h"


// This function contains our program's core functionality:

// Creating struct for eddy current data
struct eddyCurrentParams {
  float amplitude;
  float gradient;
};

// Load the data
std::vector<eddyCurrentParams> eddy_components;
std::vector<float> desired_gradient;


void loadGradient (const std::string& filename, std::vector<float>& gradient){
  std::ifstream infile(filename);
  if (!infile){
    throw std::runtime_error("Failed to open file: " + filename);
  }

  float value;
  while (infile >> value){
    gradient.push_back(value);
  }
}

void displayGradient(const std::vector<float>& gradient) {
  if (gradient.empty()) {
      std::cerr << "Gradient data is empty. Nothing to display.\n";
      return;
  }
  // Compute the minimum and maximum values in the gradient
  auto [min_it, max_it] = std::minmax_element(gradient.begin(), gradient.end());
  float min_val = *min_it;
  float max_val = *max_it;

  // Display the gradient using TG::plot()
  TG::plot(1200, 600)
    .set_ylim(min_val * 2, max_val * 2)
    .add_line(gradient)
    .show();
}


void run (std::vector<std::string>& args)
{
  debug::verbose = std::erase(args, "-v");

  if (args.size() < 3)
    throw std::runtime_error ("missing arguments - expected at least 2 argument");

  std::cerr << "reading file \"" << args[1] << "\"\n";

  std::cerr << "reading file \"" << args[2] << "\"\n";
  loadGradient(args[2], desired_gradient);

  displayGradient(desired_gradient);
}


// skeleton main() function, whose purpose is now to pass the arguments to
// run() in the expected format, and catch and handle any exceptions that may
// be thrown.

int main (int argc, char* argv[])
{
  try {
    std::vector<std::string> args (argv, argv+argc);
    run (args);
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
