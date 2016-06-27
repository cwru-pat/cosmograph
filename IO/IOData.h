#ifndef COSMO_UTILS_IODATA_H
#define COSMO_UTILS_IODATA_H

#include <string>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <fstream>

#define COSMO_IODATA_VERBOSITY_OFF 0
#define COSMO_IODATA_VERBOSITY_ON 1
#define COSMO_IODATA_VERBOSITY_DEBUG 2

namespace cosmo
{

/** I/O class **/
class IOData
{
  private:
    int verbosity;
    std::string output_dir;
    std::ofstream logfile;

    idx_t slice_output_interval;
    idx_t grid_output_interval;
    idx_t meta_output_interval;
    idx_t spec_output_interval;

    void _init(std::string output_dir_in, int verbosity_in)
    {
      output_dir = output_dir_in;
      verbosity = verbosity_in;
      size_t len_dir_name = output_dir.length();

      std::string log_filename = "log.txt";

      // default output directory name
      if(len_dir_name == 0)
        output_dir = "output";

      // temporarily remove trailing slash
      if(output_dir[len_dir_name - 1] == '/')
      {
        output_dir = output_dir.substr(0, len_dir_name - 1);
      }

      // check for conflicting data in directory
      if(std::ifstream(output_dir + "/" + log_filename))
      {
        std::cout << "Data files in output directory seem to already exist!";
        int s=1;
        while(std::ifstream(
            output_dir + "." + std::to_string(s) + "/" + log_filename
          ))
          s += 1;

        output_dir = output_dir + "." + std::to_string(s);
        std::cout << " Using '" << output_dir << "' instead.\n";
      }

      mkdir(output_dir.c_str(), 0755);

      // add in trailing slash
      output_dir += '/';

      // open log file
      logfile.open(output_dir + log_filename);
      log("Log file open.");
    }

  public:
    IOData(std::string output_dir_in)
    {
      _init(output_dir_in, COSMO_IODATA_VERBOSITY_ON);
    }

    IOData(std::string output_dir_in, int verbosity_in)
    {
      _init(output_dir_in, verbosity_in);
    }

    ~IOData()
    {
      logfile.close();
    }

    /**
     * @brief Write message to log file; console out if desired 
     */
    void log(std::string message)
    {
      if(verbosity >= COSMO_IODATA_VERBOSITY_ON)
      {
        std::cout << message << "\n" << std::flush;
      }
      logfile << message << "\n";
      logfile.flush();
    }

    /**
     * @brief Only write if in debug mode 
     */
    void debug(std::string message)
    {
      if(verbosity >= COSMO_IODATA_VERBOSITY_DEBUG)
      {
        std::cout << message << "\n" << std::flush;
        logfile << message << "\n";
        logfile.flush();
      }
    }

    /**
     * @brief Copy a file to the output dir.
     */
    void backupFile(std::string file)
    {
      std::ifstream source(file, std::ios::binary);
      std::ofstream dest(output_dir + file, std::ios::binary);
      dest << source.rdbuf();
      source.close();
      dest.close();
    }

    /**
     * @brief Return output directory
     */
    std::string dir()
    {
      return output_dir;
    }

};

template <typename T>
std::string stringify(const T value)
{
    std::ostringstream out;
    out << std::setprecision(18) << value;
    return out.str();
}

} // namespace cosmo

#endif
