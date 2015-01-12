#ifndef NWEB23
#define NWEB23

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <thread>
#include <string>
#define VERSION   23
#define BUFSIZE   8096
#define ERROR     42
#define LOG       44
#define FORBIDDEN 403
#define NOTFOUND  404
#define PORT      8080

namespace cosmo
{

/** very simple web server class, based on NWeb **/
class webServer
{
  real_t running;
  std::string logfile;

  void logger(int type, std::string s1, std::string s2, int socket_fd);
  void serve(int fd, int hit);

public:

  webServer();
  ~webServer();

  void server_start();

};

}

#endif
