
#include "nweb23.h"

namespace cosmo
{

void webServer::logger(int type, std::string s1, std::string s2, int socket_fd)
{
  int fd ;
  char logbuffer[BUFSIZE*2];

  switch (type) {
  case ERROR: (void)sprintf(logbuffer,"ERROR: %s:%s Errno=%d exiting pid=%d", s1.c_str(), s2.c_str(), errno, getpid());
    break;
  case FORBIDDEN: 
    (void)write(socket_fd, "HTTP/1.1 403 Forbidden\nContent-Length: 185\nConnection: close\nContent-Type: text/html\n\n<html><head>\n<title>403 Forbidden</title>\n</head><body>\n<h1>Forbidden</h1>\nThe requested URL, file type or operation is not allowed on this simple static file webserver.\n</body></html>\n",271);
    (void)sprintf(logbuffer,"FORBIDDEN: %s:%s", s1.c_str(), s2.c_str()); 
    break;
  case NOTFOUND: 
    (void)write(socket_fd, "HTTP/1.1 404 Not Found\nContent-Length: 136\nConnection: close\nContent-Type: text/html\n\n<html><head>\n<title>404 Not Found</title>\n</head><body>\n<h1>Not Found</h1>\nThe requested URL was not found on this server.\n</body></html>\n",224);
    (void)sprintf(logbuffer,"NOT FOUND: %s:%s",s1.c_str(), s2.c_str()); 
    break;
  case LOG: (void)sprintf(logbuffer," INFO: %s:%s:%d", s1.c_str(), s2.c_str(), socket_fd); break;
  } 
  /* No checks here, nothing can be done with a failure anyway */
  if((fd = open(logfile, O_CREAT| O_WRONLY | O_APPEND,0644)) >= 0) {
    (void)write(fd,logbuffer,strlen(logbuffer)); 
    (void)write(fd,"\n",1);      
    (void)close(fd);
  }
  if(type == ERROR || type == NOTFOUND || type == FORBIDDEN) exit(3);
}

/* this is a child web server process, so we can exit on errors */
void webServer::serve(int fd, int hit)
{
  int j, buflen;
  long i, ret, len;
  char * fstr;
  static char buffer[BUFSIZE+1]; /* static so zero filled */

  ret =read(fd,buffer,BUFSIZE);   /* read Web request in one go */
  if(ret == 0 || ret == -1) { /* read failure stop now */
    logger(FORBIDDEN,"failed to read browser request","",fd);
  }
  if(ret > 0 && ret < BUFSIZE)  /* return code is valid chars */
    buffer[ret]=0;    /* terminate the buffer */
  else buffer[0]=0;
  for(i=0;i<ret;i++)  /* remove CF and LF characters */
    if(buffer[i] == '\r' || buffer[i] == '\n')
      buffer[i]='*';
  logger(LOG,"request",buffer,hit);
  if( strncmp(buffer,"GET ",4) && strncmp(buffer,"get ",4) ) {
    logger(FORBIDDEN,"Only simple GET operation supported",buffer,fd);
  }

  /* write out something... */
  std::string out = "WELL HELLO THERE";
  (void)write(fd, out.c_str(), strlen(out.c_str()));
  
  sleep(1); /* allow socket to drain before signalling the socket is closed */
  close(fd);
  exit(1);
}

void webServer::server_start()
{
  int i, port, pid, listenfd, socketfd, hit;
  socklen_t length;
  static struct sockaddr_in cli_addr; /* static = initialised to zeros */
  static struct sockaddr_in serv_addr; /* static = initialised to zeros */

  /* Become deamon + unstopable and no zombies children (= no wait()) */
  if(fork() != 0)
    return 0; /* parent returns OK to shell */
  (void)signal(SIGCLD, SIG_IGN); /* ignore child death */
  (void)signal(SIGHUP, SIG_IGN); /* ignore terminal hangups */
  for(i=0;i<32;i++)
    (void)close(i);   /* close open files */
  (void)setpgrp();    /* break away from process group */
  logger(LOG, "nweb starting", std::to_string(PORT), getpid());

  /* setup the network socket */
  if((listenfd = socket(AF_INET, SOCK_STREAM,0)) <0)
    logger(ERROR, "system call", "socket",0);
  port = PORT;
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  serv_addr.sin_port = htons(port);
  if(bind(listenfd, (struct sockaddr *)&serv_addr,sizeof(serv_addr)) <0)
    logger(ERROR,"system call","bind",0);
  if( listen(listenfd,64) <0)
    logger(ERROR,"system call","listen",0);

  // run webserver...
  for(hit=1; ;hit++) {
    length = sizeof(cli_addr);
    if((socketfd = accept(listenfd, (struct sockaddr *)&cli_addr, &length)) < 0)
      logger(ERROR,"system call","accept",0);
    if((pid = fork()) < 0) {
      logger(ERROR,"system call","fork",0);
    }
    else {
      if(pid == 0) {  /* child */
        (void)close(listenfd);
        serve(socketfd, hit); /* never returns */
      } else {  /* parent */
        (void)close(socketfd);
      }
    }
  }
}

}
