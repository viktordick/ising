#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <sstream>
#include <boost/filesystem.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include "extent.h"
#include "random.h"

namespace fs = boost::filesystem;


int dec(int x, int N=L) { //decrease index with periodic boundaries
    return (x==0)?(N-1):(x-1);
}
int inc(int x, int N=L){
    return (x==N-1)?0:(x+1);
}

static std::mt19937_64 *pRND;
#ifdef _OPENMP
#pragma omp threadprivate(pRND)
#endif

class Random {
    private:
        static bool initialized;
        //         static std::uniform_real_distribution<double> dst;
        Random();
    public:
        static void init(std::mt19937_64::result_type seed) {
            if (initialized)
                return;
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
#ifdef _OPENMP
                seed += omp_get_thread_num();
#endif
                pRND = new std::mt19937_64(seed);
                pRND->discard(10);
            }
            initialized = true;
        }
        static double get() {
            return (*pRND)()/(double(std::mt19937_64::max())+1);
        }
        static auto bits() {
            return random_bits(*pRND);
        }
};

bool Random::initialized = false;
// std::uniform_real_distribution<double> Random::dst{0,1};

#if 1

struct Line { //one line of spin variables on half lattice
    private:
        static const int LH = L/2;
        //use shortest unsigned integer type that has enough bits to hold LH
        typedef unsigned char INT0;
        typedef unsigned int INT1;
        typedef unsigned long int INT2;
        typedef unsigned long long int INT3;

        typedef std::conditional< 8*sizeof(INT0)>=LH, INT0,
                std::conditional< 8*sizeof(INT1)>=LH, INT1, 
                std::conditional< 8*sizeof(INT2)>=LH, INT2, 
                INT3>::type>::type>::type T;
        static const int B = 8*sizeof(T); //how many bits in one element
        static const int N = 1+(LH-1)/B; //how many elements we need for one line
        //constant line where every relevant bit is set
        static const T FullLine = (~T(0)<<(N*B-LH))>>(N*B-LH);
        //bit for last column is least significant, i.e. value '1'
        T dat[N];
    public:
        void init(std::mt19937_64 &rnd) {
            for (int i=0; i<N; i++)
                dat[i] = rnd();
            dat[0] &= FullLine;
        }
        //only keep 'up' bits with probability exp(-beta dE) 
        void randomize() {
            for (int i=0; i<N; i++)
                dat[i] &= ~Random::bits();
        }
        Line() {
            for (int i=0; i<N; i++)
                dat[i] = 0;
        }
        void operator^=(const Line &a) {
            for (int i=0; i<N; i++)
                dat[i] ^= a.dat[i];
        }
        bool get(int i) {
            const int j = i/B;
            i -= j*B;
            return (dat[j]>>(i-1))&1;
        }
        void set(int i, bool up) {
            const int j = i/B;
            i -= j*B;
            if (up)
                dat[j] |= (T(1)<<i);
            else
                dat[j] &= ~(T(1)<<i);
        }
        Line operator&(const Line &l) const {
            Line result;
            for (int i=0; i<N; i++)
                result.dat[i] = dat[i]&l.dat[i];
            return result;
        }
        Line operator|(const Line &l) const {
            Line result;
            for (int i=0; i<N; i++)
                result.dat[i] = dat[i]|l.dat[i];
            return result;
        }
        Line operator^(const Line &l) const {
            Line result;
            for (int i=0; i<N; i++)
                result.dat[i] = dat[i]^l.dat[i];
            return result;
        }
        Line operator~() const { //only relevant bits are negated
            Line result;
            result.dat[0] = FullLine & ~dat[0];
            for (int i=1; i<N; i++)
                result.dat[i] = ~dat[i];
            return result;
        }
        Line operator|=(const Line &l) {
            for (int i=0; i<N; i++)
                dat[i] |= l.dat[i];
            return *this;
        }
        Line operator&=(const Line &l) {
            for (int i=0; i<N; i++)
                dat[i] &= l.dat[i];
            return *this;
        }
        // bit for column 0 becomes bit for column 2
        Line shift_right() const { 
            Line result;
            result.dat[0] = (dat[0]>>1) | ((dat[N-1]&1)<<((LH-1)%B));
            for (int i=1; i<N; i++)
                result.dat[i] = (dat[i]>>1) | ((dat[i-1]&1)<<(B-1));
            return result;
        }
        Line shift_left() const {
            Line result;
            for (int i=0; i<N-1; i++)
                result.dat[i] = (dat[i]<<1) | (dat[i+1]>>(B-1));
            result.dat[N-1] = (dat[N-1]<<1) | (dat[0]>>((LH-1)%B));
            result.dat[0] &= FullLine;
            return result;
        }
        Line shift(bool right) const {
            return right?shift_right():shift_left();
        }

        int count() {
            int result = 0;
            for (int i=0; i<N; i++)
                for (T mask = T(1); mask; mask <<= 1)
                    if (dat[i] & mask)
                        result++;
            return result;
        }
        friend std::ostream &operator<<(std::ostream &os, const Line &l) {
            for (T mask = T(1)<<((LH-1)%B); mask; mask>>=1)
                os << ((l.dat[0]&mask)?'1':'0');
            for (int i=1; i<N; i++)
                for (T mask = T(1)<<(B-1); mask; mask>>=1)
                    os << ((l.dat[i]&mask)?'1':'0');
            return os;
        };
        friend std::istream &operator>>(std::istream &is, Line &l) {
            char c;
            l.dat[0] = 0;
            for (int j=0; j<=(LH-1)%B; j++) {
                is.get(c);
                l.dat[0] <<= 1;
                if (c=='1')
                    l.dat[0] |= 1;
            }
            for (int i=1; i<N; i++) {
                l.dat[i] = 0;
                for (int j=0; j<B; j++) {
                    is.get(c);
                    l.dat[i] <<= 1;
                    if (c=='1')
                        l.dat[i] |= 1;
                }
            }
            return is;
        }

};

class Ising {
    private:
        Line dat[2][L];
    public:
        Ising(unsigned seed)
        {
            Random::init(seed);
        }
        void init() {
            std::mt19937_64 rnd;
            for (int eo=0; eo<2; eo++)
                for (int x=0; x<L; x++)
                    dat[eo][x].init(rnd);
        }
        void sweep() {
            for (int eo=0; eo<2; eo++) { //sublattice, 0=even, 1=odd

                //depends on if that line represents columns 024... or 135...
                //maybe it should be !eo, but that just amounts to redefining which lattice is
                //even and which odd
                bool right = eo; 

                for (int x=0; x<L; x++) {
                    Line &c = dat[eo][x];
                    Line *o = dat[1-eo];
                    //each line represents one neighbor (nb.), with bit set if that neighbor
                    //is antiparallel (a.p.) to c
                    Line n[4] = {
                        c^o[dec(x)], //line above
                        c^o[x],  //same line, without shift
                        c^o[x].shift(right),  //same line, with shift
                        c^o[inc(x)] //line below
                    };
                    //count how many nb. are a.p. - first counting 'at least 1' and 'at least 2'
                    Line c1, c234;
                    for (int i=0; i<4; i++) {
                        c234 |= c1&n[i]; //has at least two a.p. neighbors
                        c1 |= n[i]; // has at least one a.p. neighbor
                    }
                    //now compute 'exactly 0' and 'exactly 1'
                    Line c0 = ~c1;
                    c1 &= ~c234;

                    c1.randomize(); //we want to flip these bits with some probability
                    c0.randomize();
                    c0.randomize(); //here the probability is exp(-8beta), so we apply twice
                    c ^= c0|c1|c234;
                    right = !right;
                }
            }
        }
        double measure() {
            int result = 0;
            for (int eo=0; eo<2; eo++)
                for (int x=0; x<L; x++)
                    result += dat[eo][x].count();
            return fabs(double(result*2)/(L*L)-1);
        }
        void save() const {
            std::stringstream fname;
            fname << ".state/" << L;
            fs::create_directories(fname.str());
            fname << "/" << std::fixed << std::setprecision(5) << beta;
            std::ofstream f(fname.str(), std::ofstream::ate);
            for (int x=0; x<L; x++) 
                f << dat[0][x] << ' ' << dat[1][x] << std::endl;
        }
        bool load() {
            std::stringstream fname;
            fname << ".state/" << L << "/" << std::fixed << std::setprecision(5) << beta;
            std::ifstream f(fname.str());
            if (f.fail()) {
                return false;
            }
            for (int x=0; x<L; x++) {
                f >> dat[0][x];
                f.get();
                f >> dat[1][x];
                f.get();
            }
            return !f.fail();
        }

        friend std::ostream &operator<<(std::ostream &os, const Ising &ising) {
            for (int x=0; x<L; x++)
                os << ising.dat[0][x] << ' ' << ising.dat[1][x] << std::endl;
            return os;
        }
};

#else
struct Ising {
    char data[L][L];
    double weights[2]; //exp(-4*n*beta)

    Ising(unsigned seed) {
        Random::init(seed);
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

#endif

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
    int nmeas;
    s >> nmeas;
    if (s.fail()) {
        std::cerr << "usage: " << argv[0] << " [measurements]" << std::endl;
        return 1;
    }
    const int extent = L;

    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    Ising ising(seed);
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

