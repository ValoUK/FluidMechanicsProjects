#include "FileParser.h"
#include <iostream> 
#include <fstream> 
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

std::map<std::string, std::string> FileParser::readFile(std::string &filePath)
{
	std::map<std::string, std::string> map;
	std::ifstream inputFile(filePath.c_str());
	std::string line,key,value;
	if (inputFile.is_open())
	{
		while (std::getline(inputFile, line))
		{
			if (boost::algorithm::contains(line, "#"))
			{
				key = line;
				boost::erase_all(key, "#");
				boost::trim(key);
			}
			else
			{
				value = line;
				boost::trim(value);
				map[key] = value;
			}
		}
	} 
	return map;
}
