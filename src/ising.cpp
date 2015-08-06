#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <sstream>
#include <boost/filesystem.hpp>

#include <omp.h>
#include <math.h>
#include "extent.h"

namespace fs = boost::filesystem;


int dec(int x, int N=L) { //decrease index with periodic boundaries
    return (x==0)?(N-1):(x-1);
}
int inc(int x, int N=L){
    return (x==N-1)?0:(x+1);
}

static std::default_random_engine *pRND;
#ifdef _OPENMP
#pragma omp threadprivate(pRND)
#endif

class Random {
    private:
        static bool initialized;
        static std::uniform_real_distribution<double> dst;
        Random();
    public:
    static void init(std::default_random_engine::result_type seed) {
        if (initialized)
            return;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            pRND = new std::default_random_engine(seed+omp_get_thread_num());
            pRND->discard(10);
        }
        initialized = true;
    }
    static double get() {
        if (!initialized)
            throw "forgot to call Random::init()";
        return dst(*pRND);
    }
};

bool Random::initialized = false;
std::uniform_real_distribution<double> Random::dst{0,1};

struct Ising {
    char data[L][L];
    double beta; 
    double weights[2]; //exp(-4*n*beta)

    Ising(unsigned seed, double _beta) {
        Random::init(seed);
        beta = _beta;
        weights[0] = exp(-4*beta); 
        weights[1] = pow(weights[0],2);
    }
    void init() {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int x=0; x<L; x++)
            for (int y=0; y<L; y++)
                data[x][y] = (Random::get()<0.5)?1:-1;
    }
    void sweep() {
        for (int eo=0; eo<2; eo++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int x=0; x<L; x++) {
                for (int y=(x+eo)%2; y<L; y+=2) {
                    char &cur = data[x][y];
                    char N = 
                        data[dec(x)][y] + data[inc(x)][y] +
                        data[x][dec(y)] + data[x][inc(y)];
                    N *= cur;
                    if (N <= 0)
                        cur = -cur;
                    else {
                        double r = Random::get();
                        if (r<weights[N/2-1])
                            cur = -cur;
                    }
                }
            }
            }
        }
        double measure() {
            int result = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
            for (int x=0; x<L; x++)
                for (int y=0; y<L; y++)
                    result += data[x][y];
            return fabs(result)/(L*L);
        }
        void save() const {
            std::stringstream fname;
            fname << ".state/" << L;
            fs::create_directories(fname.str());
            fname << "/" << std::fixed << std::setprecision(5) << beta;
            std::ofstream f(fname.str(), std::ofstream::ate);
            for (int x=0; x<L; x++) {
                for (int y=0; y<L; y++)
                    f << ((data[x][y]==1)?'1':'0');
                f << std::endl;
            }
        }
        bool load() { 
            std::stringstream fname;
            fname << ".state/" << L << "/" << std::fixed << std::setprecision(5) << beta;
            std::ifstream f(fname.str());
            if (f.fail()) {
                return false;
            }
            char c;
            for (int x=0; x<L; x++) {
                for (int y=0; y<L; y++) {
                    f.get(c);
                    if (c!='1' && c!= '0')
                        return false;
                    data[x][y] = (c=='1')?1:-1;
                }
                f.get(c); //endline
            }
            return !f.fail();
        }

    };


    template<typename T>
        void write(std::ostream &f, const T &x) {
            f.write((char*)&x, sizeof(T));
        }
    int main(int argc, char** argv)
    {
        srand48(time(NULL));
        std::stringstream s;
        for (int i=1; i<argc; i++)
            s << argv[i] << ' ';
        double beta;
        int nmeas;
        s >> beta >> nmeas;
        if (s.fail()) {
            std::cerr << "usage: " << argv[0] << " [beta] [measurements]" << std::endl;
            return 1;
        }
        const int extent = L;

        const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        Ising ising(seed,beta);
        std::ofstream::openmode m;
        if (ising.load()) {
            std::cout << "# Resuming "; 
            m = std::ofstream::app;
        } else {
            std::cout << "# Initialized ";
            m = std::ofstream::trunc;
            ising.init();
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
        for (int i=0; i<nmeas; i++) {
            for (int j=0; j<10; j++)
                ising.sweep();
            double M = ising.measure();
            write(outfile,M);
        }
        ising.save();

    };

