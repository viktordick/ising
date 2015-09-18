#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/random.hpp>
#include <vector>
#include <sstream>

#include <math.h>
#include <signal.h>
#include <sys/stat.h>


typedef float floatT;

off_t getFilesize(const char *path){
    struct stat fStat;
    if (!stat(path, &fStat)) 
        return fStat.st_size;
    else
        return 0;
};

int dec(int x, int N=L) { //decrease index with periodic boundaries
    return (x==0)?(N-1):(x-1);
}
int inc(int x, int N=L){
    return (x==N-1)?0:(x+1);
}

template <class T>
void write(std::ostream &out, T x) {
    out.write((char*)&x, sizeof(T));
}

#ifdef GDEV
class Random {
    private:
        boost::mt19937 base;
    public:
        typedef uint64_t result_type;
        void seed(result_type _seed) {
            base.seed(_seed);
        }
        Random(result_type seed):
            base(seed) {};
        result_type operator()() {
            return result_type(base())<<32 | base();
        }
};
#else
typedef boost::mt19937_64 Random;
#endif

static volatile bool keepRunning = true;
void intHandler(int dummy) {
    keepRunning = false;
}

#include "line.h"
#include "ising.h"
#include "random.h"


int main(int argc, char** argv)
{
    signal(SIGINT, intHandler);
    signal(SIGTERM, intHandler);

    std::stringstream s;
    for (int i=1; i<argc; i++)
        s << argv[i] << ' ';
    std::string sig;
    int nmeas;
    s >> sig >> nmeas;
    if (s.fail()) {
        std::cerr << "usage: " << argv[0] << " sig nmeas" << std::endl;
        return 1;
    }

    time_t rawtime;
    tm *timeinfo;
    const unsigned seed = time(&rawtime);
    timeinfo = localtime(&rawtime);

    std::cout << std::setfill('0')
        << std::setw(2) << timeinfo->tm_hour << ':'
        << std::setw(2) << timeinfo->tm_min << ':'
        << std::setw(2) << timeinfo->tm_sec << ' '
        << std::setfill(' ');

    run(sig,seed,nmeas);
};

