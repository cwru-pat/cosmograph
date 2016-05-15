#ifndef COSMO_UTILS_IO_H
#define COSMO_UTILS_IO_H

#include <string>
#include <iostream>
#include <fstream>

#define COSMO_IO_VERBOSITY_OFF 0
#define COSMO_IO_VERBOSITY_ON 1
#define COSMO_IO_VERBOSITY_DEBUG 2

#define COSMO_IO_DEBUG_COUT if(verbosity >= COSMO_IO_VERBOSITY_DEBUG) std::cout
#define COSMO_IO_COUT if(verbosity >= COSMO_IO_VERBOSITY_ON) std::cout

namespace cosmo
{

/** I/O class **/
class IO
{
  private:
    int verbosity;
    std::string output_dir;
    std::ofstream log;

  public:
    IO(std::string output_dir_in)
    {
    
    }

    ~IO()
    {

    }

};

#endif
