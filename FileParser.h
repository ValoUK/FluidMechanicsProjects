#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <map> //std::map
#include <iostream> 
#include <string> //std::string
#include <vector> //std::vector
#include <fstream> //std::ifstream
#include <sstream> // std::istringstream

class FileParser
{
public:
	std::map<std::string, std::string> readFile(std::string &filePath);
private:
	std::string _filePath;
};

#endif