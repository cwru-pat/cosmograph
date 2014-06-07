#ifndef COSMO_CONFIGPARSER_H
#define COSMO_CONFIGPARSER_H

#include "../cosmo.h"
namespace cosmo
{

class ConfigParser
{
public:
  ConfigParser(std::string fname);

  std::string operator[](std::string param);

private:
  std::map<std::string, std::string> config;
};

extern ConfigParser _config;

} /* namespace */
#endif
