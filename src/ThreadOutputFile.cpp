#include <mutex>
#include <string>
#include <vector>

#include "Fragments.h"
#include "ThreadOutputFile.h"


ThreadOutputFile::ThreadOutputFile (const std::string &path) : _path(path) {
	// open the file
	this->outputFile.open(path);
}

void ThreadOutputFile::write(const std::vector<FragmentStruct> &fragments) {
	std::lock_guard<std::mutex> lock(_writerMutex);

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

ThreadOutputFile::~ThreadOutputFile() {
	this->outputFile.close();
}