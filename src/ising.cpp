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

class mt19937_64 {
    private:
        boost::mt19937 base;
    public:
        typedef uint64_t result_type;
        mt19937_64(result_type seed):
            base(seed) {};
        result_type operator()() {
            return result_type(base())<<32 | base();
        }
};
class Random {
    private:
        typedef mt19937_64 Gen;
        typedef Gen::result_type result_type;
        Gen r;
    public:
        Random(result_type seed)
        :r(seed)
        { }
        result_type uniform() { //each bit set with 50% chance
            return r();
        }
        result_type e4beta(); //each bit set with exp(-4*beta) chance
};

#include "random.h"
#include "line.h"
#include "ising.h"


template<typename T>
void write(std::ostream &f, const T &x) {
    f.write((char*)&x, sizeof(T));
}


static volatile bool keepRunning = true;
void intHandler(int dummy) {
    keepRunning = false;
}


int main(int argc, char** argv)
{
    signal(SIGINT, intHandler);
    signal(SIGTERM, intHandler);

    srand48(time(NULL));
    std::stringstream s;
    for (int i=1; i<argc; i++)
        s << argv[i] << ' ';
    int nmeas;
    s >> nmeas;
    if (s.fail()) {
        std::cerr << "usage: " << argv[0] << " [measurements]" << std::endl;
        return 1;
    }
    const int extent = L;


    const unsigned seed = time(NULL);
    Ising ising(seed);
    std::ofstream::openmode m;

    std::stringstream out_file_name; 
    mkdir("data", 0775);
    out_file_name << "data/" << std::setfill('0') << std::setw(3) << extent;
    mkdir(out_file_name.str().c_str(), 0775);
    out_file_name << "/" << std::fixed << std::setprecision(10) << beta;
    int measured = getFilesize(out_file_name.str().c_str())/sizeof(floatT);

    if (ising.load()) {
        std::cout << "# Resume "; 
        m = std::ofstream::app;
    } else {
        std::cout << "#   Init ";
        m = std::ofstream::trunc;
        if (beta < 0.44)
            ising.randomize();
    }
    std::cout << "L=" << extent << ", beta=" << std::fixed << beta << ", Nmeas=";

    if (m == std::ofstream::app)
        std::cout << std::max(nmeas-measured,0) << '/';
    std::cout << nmeas << std::endl;

    std::ofstream outfile(out_file_name.str().c_str(), m);
    if (outfile.fail()) {
        std::cerr << "Output file could not be opened. Aborting." << std::endl;
        return -1;
    }
    for (; keepRunning && measured<nmeas; measured++) {
        for (int j=0; j<L; j++)
            ising.sweep();
        floatT M = ising.magnetization();
        write(outfile,M);
    }
    ising.save();
};

