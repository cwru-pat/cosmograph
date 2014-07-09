#ifndef COSMO_CONFIGPARSER_H
#define COSMO_CONFIGPARSER_H

#include <string>
#include <map>

namespace cosmo
{

class ConfigParser
{
public:
  ConfigParser();
  ConfigParser(std::string fname);
  void parse(std::string fname);

  std::string operator[](std::string param);

private:
  std::map<std::string, std::string> config;
};

} /* namespace */
#endif
