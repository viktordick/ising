#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <sstream>
#include <boost/filesystem.hpp>
#include <math.h>
#include <signal.h>

#include "extent.h"

namespace fs = boost::filesystem;

int dec(int x, int N=L) { //decrease index with periodic boundaries
    return (x==0)?(N-1):(x-1);
}
int inc(int x, int N=L){
    return (x==N-1)?0:(x+1);
}

class Random {
    private:
        typedef std::mt19937_64 Gen;
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

    auto now = std::chrono::system_clock::now();
    auto tt = std::chrono::system_clock::to_time_t(now);
    auto tm = std::tm{0};
    localtime_r(&tt, &tm);
    std::cout.fill('0');
    std::cout << std::setw(2) << tm.tm_hour 
        << ':' << std::setw(2) << tm.tm_min 
        << ':' << std::setw(2) << tm.tm_sec << ' ';
    const unsigned seed = now.time_since_epoch().count();
    Ising ising(seed);
    std::ofstream::openmode m;
    if (ising.load()) {
        std::cout << "# Resuming "; 
        m = std::ofstream::app;
    } else {
        std::cout << "# Initialized ";
        m = std::ofstream::trunc;
        ising.randomize();
    }
    std::cout << "ising with L=" << extent << ", beta=" << std::fixed << beta
        << ", Nmeas=" << nmeas << std::endl;

    std::stringstream out_file_name; 
    out_file_name << "data/" << std::setfill('0') << std::setw(3) << extent;
    fs::create_directories(out_file_name.str());
    out_file_name << "/" << std::fixed << std::setprecision(10) << beta;
    std::ofstream outfile(out_file_name.str().c_str(), m);
    if (outfile.fail()) {
        std::cerr << "Output file could not be opened. Aborting." << std::endl;
        return -1;
    }
    if (m == std::ofstream::trunc)
        write(outfile, beta);
    for (int i=0; keepRunning && i<nmeas; i++) {
        for (int j=0; j<10; j++)
            ising.sweep();
        double M = ising.magnetization();
        write(outfile,M);
    }
    ising.save();
};

