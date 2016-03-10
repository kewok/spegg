#ifndef INPUT_FILE_CHECKER_H
#define INPUT_FILE_CHECKER_H

#include <algorithm>
#include <istream>
#include <fstream>

int count_lines_in_file(const char *filename)
	{
	int answer;
	std::ifstream inFile(filename);
	answer = std::count(std::istreambuf_iterator<char>(inFile),
		   	    std::istreambuf_iterator<char>(),'\n');
	return answer;
	}

#endif
