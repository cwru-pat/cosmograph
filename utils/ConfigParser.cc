
#include "ConfigParser.h"

#include <iostream>
#include <iomanip>
#include <fstream>

namespace cosmo
{
ConfigParser::ConfigParser()
{
}

ConfigParser::ConfigParser(std::string fname)
{
  parse(fname);
}

void ConfigParser::parse(std::string fname)
{
  fileName = fname;
  
  /* param = val */
  std::string param;
  std::string eq;
  std::string val;

  std::ifstream fin(fname.c_str());
  if(!fin) {
    std::cerr << "error opening " << fname
      << "! Make sure you have specified a valid config file.";
    std::cerr << std::endl;
    throw -1;
    return;
  }

  fin >> param >> eq >> val;
  while(fin) {
    std::cout << "config[" << param << "] = " << val << std::endl;
    config[param] = val;
    fin >> param >> eq >> val;
  }
}

std::string ConfigParser::getFileName()
{
  return fileName;
}

std::string ConfigParser::operator[](std::string param)
{
  if( config.find(param) == config.end() )
  {
    std::cout << "Error: param `" << param
      << "` is required, but was not set!"
      << " Please set this in the configuration file, " << fileName << ".\n";
    throw -1;
  }

  return config[param];
}

} /* namespace */
