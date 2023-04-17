#ifndef THREAD_OUTPUT_FILE_H
#define THREAD_OUTPUT_FILE_H

#include <vector>
#include <string>
#include <mutex>
#include <fstream>

#include "Fragments.h"


class ThreadOutputFile {
	public:
		ThreadOutputFile();
		ThreadOutputFile (std::string, int);

		void open();
		void setFile(std::string);
        void write(const std::vector<FragmentStruct> &);
		void debugWrite(const std::string &);
		void debugWrite(std::vector<const char *>);
		~ThreadOutputFile();
	
		std::string getPath() const;
	private:
		std::string _path;
		std::ofstream outputFile;
};

#endif
