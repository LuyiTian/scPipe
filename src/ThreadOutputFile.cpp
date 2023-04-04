#include <mutex>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "Fragments.h"
#include "ThreadOutputFile.h"

ThreadOutputFile::ThreadOutputFile() {
	// Rcpp::Rcout << "Created ThreadOutputFile\n";
}

ThreadOutputFile::ThreadOutputFile (std::string path, int i) {
	this->_path = path;
	// Rcpp::Rcout << "Created ThreadOutputFile: " << i << "\n";
}

void ThreadOutputFile::setFile (std::string file) {
	this->_path = file;
}

void ThreadOutputFile::open() {
	if (!this->outputFile.is_open()) {
		// open the file
		this->outputFile.open(this->_path, std::ofstream::out | std::ofstream::app);
	}
}

void ThreadOutputFile::write(const std::vector<FragmentStruct> &fragments) {
	// Rcpp::Rcout << "writing to out file.\n";
	if (!this->outputFile.is_open()) {
		// open the file
		this->outputFile.open(this->_path, std::ofstream::out | std::ofstream::app);
	}
	// print all the fragments to the file, 
	// con	catenating each attribute (except for complete flag)
	// and joining with tab characters
	for (auto frag : fragments) {
		this->outputFile << frag.chromosome << "\t"
			<< frag.start << "\t"
			<< frag.end << "\t"
			<< frag.cell_barcode << "\t"
			<< frag.sum << "\n";
	}
}

void ThreadOutputFile::debugWrite(const std::string &line) {
	this->open();
	this->outputFile << line << "\n";
}

void ThreadOutputFile::debugWrite(std::vector<const char *> lines) {
	this->open();
	
	for (auto s : lines) {
		this->outputFile << s << "\t";
	}
	this->outputFile << "\n";
}

std::string ThreadOutputFile::getPath() const {
	return this->_path;
}

ThreadOutputFile::~ThreadOutputFile() {
	// Rcpp::Rcout << "Destroyed ThreadOutputFile\n";
	if (this->outputFile.is_open()) this->outputFile.close();
}