#ifndef THREAD_OUTPUT_FILE_H
#define THREAD_OUTPUT_FILE_H

#include <vector>
#include <string>
#include <mutex>
#include <fstream>

#include "Fragments.h"


class ThreadOutputFile {
	public:
		ThreadOutputFile (const std::string &);

		void write(const std::vector<FragmentStruct> &);

		~ThreadOutputFile();
	
	private:
		std::string _path;
		std::mutex _writerMutex;
		std::ofstream outputFile;
};

#endif
