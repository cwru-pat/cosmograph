
#include "ConfigParser.h"

namespace cosmo
{

ConfigParser::ConfigParser(std::string fname)
{
  /* param = val */
  std::string param;
  std::string eq;
  std::string val;

  std::ifstream fin(fname.c_str());
  if(!fin) {
    std::cerr << "error opening " << fname << ", using default values ";
    std::cerr << std::endl;
    return;
  }

  fin >> param >> eq >> val;
  while(fin) {
    std::cout << "params[" << param << "] = " << val << std::endl;
    config[param] = val;
    fin >> param >> eq >> val;
  }
}

std::string ConfigParser::operator[](std::string param)
{
  return config[param];
}

} /* namespace */
