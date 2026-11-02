// Linux/WSL measurement helper: exact wait4 user/sys and peak RSS for the child.
// Its own monitoring time is excluded. Every output is created exclusively.
#include <chrono>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fcntl.h>
#include <signal.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {
volatile sig_atomic_t stop_signal = 0;
void request_stop(int signal) { stop_signal = signal; }
}

int main(int argc, char ** argv) {
  unsigned long timeout_seconds = 900;
  if (argc > 1 && std::strcmp(argv[1], "--timeout-seconds") == 0) {
    if (argc < 6 || !argv[2][0]) return 2;
    for (const char *p=argv[2]; *p; ++p) if (*p<'0' || *p>'9') return 2;
    errno=0;
    timeout_seconds=std::strtoul(argv[2],nullptr,10);
    if (errno) return 2;
    argc-=2; argv+=2;
  }
  if (argc < 4) { std::fprintf(stderr,"usage: measure [--timeout-seconds N] timing log command [args...] (N=0: unlimited; default: 900)\n"); return 2; }
  const int timing = open(argv[1], O_WRONLY|O_CREAT|O_EXCL, 0644);
  if (timing < 0) { std::perror("timing output"); return 2; }
  const int log = open(argv[2], O_WRONLY|O_CREAT|O_EXCL, 0644);
  if (log < 0) { std::perror("log output"); close(timing); return 2; }
  const auto start = std::chrono::steady_clock::now();
  struct sigaction action{};
  action.sa_handler=request_stop;
  sigemptyset(&action.sa_mask);
  sigaction(SIGINT,&action,nullptr);
  sigaction(SIGTERM,&action,nullptr);
  const pid_t pid = fork();
  if (pid < 0) { std::perror("fork"); close(log); close(timing); return 2; }
  if (!pid) {
    signal(SIGINT,SIG_DFL);
    signal(SIGTERM,SIG_DFL);
    setpgid(0,0);
    close(timing);
    if (dup2(log,STDOUT_FILENO)<0 || dup2(log,STDERR_FILENO)<0) _exit(126);
    close(log);
    execvp(argv[3],argv+3);
    std::perror("execvp"); _exit(127);
  }
  close(log);
  setpgid(pid,pid);
  int status=0;
  rusage usage{};
  bool timeout=false;
  bool stopping=false;
  double stop_started=0;
  for (;;) {
    const auto result=wait4(pid,&status,WNOHANG,&usage);
    if (result==pid) break;
    if (result<0 && errno!=EINTR) {
      std::perror("wait4"); kill(-pid,SIGKILL); close(timing); return 2;
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if (!stopping && (stop_signal || (timeout_seconds && seconds>timeout_seconds))) {
      timeout=!stop_signal;
      stopping=true;
      stop_started=seconds;
      kill(-pid,stop_signal ? stop_signal : SIGTERM);
    }
    if (stopping && seconds>stop_started+5) kill(-pid,SIGKILL);
    usleep(10000);
  }
  const double real=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  const double user=usage.ru_utime.tv_sec+usage.ru_utime.tv_usec/1e6;
  const double sys=usage.ru_stime.tv_sec+usage.ru_stime.tv_usec/1e6;
  const int code=timeout ? 124 : stop_signal ? 128+stop_signal : WIFEXITED(status) ? WEXITSTATUS(status) : 128+WTERMSIG(status);
  dprintf(timing,"real_seconds=%.6f user_seconds=%.6f sys_seconds=%.6f max_rss_kib=%ld exit_code=%d\n",
          real,user,sys,usage.ru_maxrss,code);
  close(timing);
  return code;
}
