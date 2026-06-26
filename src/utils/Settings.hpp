#pragma once
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include "globals.hpp"
#include "helpString.hpp"

namespace element_validity {
class Settings {
	public:
	std::string filePath = "";
	std::string resultsPath = "";
	bool globalQuery = false;
	bool staticCheck = false;
	fp_t staticCheckTime = 1.;
	int maxIterations = 0;
	bool preCheck = false;
	int preCheckMaxIter = 0;
	int numElem = 0;
	int firstElem = 0;
	int numThreads = 1;
	fp_t precision = .01;
	bool abort = false;

	// Argument parser
	Settings(int argc, char** argv);
};
}