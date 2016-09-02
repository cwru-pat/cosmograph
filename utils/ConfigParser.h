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
  std::string getFileName();

  std::string operator[](std::string param);
  std::string operator()(std::string param, std::string default_val);

private:
  std::map<std::string, std::string> config;
  std::string fileName;
};

} /* namespace */
#endif
