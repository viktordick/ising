#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <sstream>
#include <boost/filesystem.hpp>

#include <math.h>
#include "extent.h"

namespace fs = boost::filesystem;

const uint LH = L/2;

int dec(int x, int N=L) { //decrease index with periodic boundaries
    return (x==0)?(N-1):(x-1);
}
int inc(int x, int N=L){
    return (x==N-1)?0:(x+1);
}

#if 0
//use shortest unsigned integer type that has enough bits to hold LH
typedef unsigned char INT0;
typedef unsigned int INT1;
typedef unsigned long int INT2;
typedef unsigned long long int INT3;

typedef std::conditional< 8*sizeof(INT0)>=LH, INT0,
        std::conditional< 8*sizeof(INT1)>=LH, INT1, 
        std::conditional< 8*sizeof(INT2)>=LH, INT2, 
        INT3>::type>::type>::type LineT;
//constant line where every relevant bit is set
const LineT FullLine = (~LineT(0)<<(8*sizeof(LineT)-LH))>>(8*sizeof(LineT)-LH);

class Random {
    private:
        std::default_random_engine &r;
        std::bernoulli_distribution d;
    public:
        Random(std::default_random_engine &r, double p)
            :r(r)
             ,d(p)
    { }
        bool operator()() { return d(r); }
};

struct Line { //one line of spin variables on half lattice
    typedef LineT T;
    //bit for last column is least significant, i.e. value '1'
    T dat;

    //only keep 'up' bits with probability exp(-beta dE) 
    void randomize(Random &r) {
        T mask = 1;
        T cp = dat;
        while (cp) {
            if ((cp & 1) && !r())
                dat ^= mask;
            cp >>= 1;
            mask <<= 1;
        }
    }
    Line() :dat(0) {}
    operator T() { return dat; } //cast to T
    Line(T d) : dat(d) {}
    void operator^=(const Line &a) {
        dat ^= a.dat;
    }
    bool get(int i) {
        return (dat>>(i-1))&1;
    }
    void set(int i, bool up) {
        if (up)
            dat |= (T(1)<<i);
        else
            dat &= ~(T(1)<<i);
    }
    Line operator&(const Line &l) const {
        return dat&l.dat;
    }
    Line operator|=(const Line &l) {
        dat |= l.dat;
        return *this;
    }
    Line operator&=(const Line &l) {
        dat &= l.dat;
        return *this;
    }
    // bit for column 0 becomes bit for column 2
    Line shift_right() const { 
        return (dat>>1) | ((dat&1)<<(LH-1));
    }
    Line shift_left() const {
        const int cut = 8*sizeof(T)-LH;
        return ((dat<<(cut+1))>>cut) | (dat>>(LH-1));
    }
    Line shift(bool right) const {
        return right?shift_right():shift_left();
    }

    int count() {
        int result = 0;
        for (T mask = T(1)<<(LH-1); mask; mask >>= 1)
            if (dat & mask)
                result++;
        return result;
    }
    friend std::ostream &operator<<(std::ostream &os, const Line &l) {
        for (T mask = T(1)<<(LH-1); mask; mask>>=1)
            os << ((l.dat&mask)?'1':'0');
        return os;
    };
    friend std::istream &operator>>(std::istream &is, Line &l) {
        char c;
        for (int i=LH-1; i>=0; i--) {
            is.get(c);
            l.set(i,c=='1');
        }
        is.get(c);
        return is;
    }

};

class Ising {
    private:
        Line dat[2][L];
        std::default_random_engine rnd;
        double beta;
        Random e4beta, e8beta;
    public:
        Ising(unsigned seed, double beta) :
            rnd(seed),
            beta(beta),
            e4beta(rnd,exp(-4*beta)),
            e8beta(rnd,exp(-8*beta))
    {
    }
        void init() {
            std::bernoulli_distribution d;
            for (int eo=0; eo<2; eo++)
                for (int x=0; x<L; x++)
                    for (int y=0; y<LH; y++)
                        dat[eo][x].set(y,d(rnd));
            //                         dat[eo][x].set(y,(x==3 && y==2 && eo==0));
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
                    Line c0 = FullLine & ~c1;
                    c1 &= ~c234;

                    c1.randomize(e4beta); //we want to flip these bits with some probability
                    c0.randomize(e8beta);
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
            for (int eo=0; eo<2; eo++)
                for (int x=0; x<L; x++)
                    f << dat[eo][x] << std::endl;
        }
        bool load() {
            std::stringstream fname;
            fname << ".state/" << L << "/" << std::fixed << std::setprecision(5) << beta;
            std::ifstream f(fname.str());
            if (f.fail()) {
                return false;
            }
            for (int eo=0; eo<2; eo++)
                for (int x=0; x<L; x++)
                    f >> dat[eo][x];
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
    std::default_random_engine rnd;
    std::uniform_real_distribution<double> dst;
    double beta, e4beta, e8beta;

    Ising(unsigned seed, double _beta):
        rnd(seed) {
            beta = _beta;
            e4beta = exp(-4*beta);
            e8beta = exp(-8*beta);
        }
    void init() {
        std::bernoulli_distribution d;
        for (int x=0; x<L; x++)
            for (int y=0; y<L; y++)
                data[x][y] = d(rnd)?1:-1;
    }
    void sweep() {
        for (int eo=0; eo<2; eo++) {
            for (int x=0; x<L; x++) {
                for (int y=(x+eo)%2; y<L; y+=2) {
                    char &cur = data[x][y];
                    char N = 
                        data[dec(x)][y] +
                        data[inc(x)][y] +
                        data[x][dec(y)] +
                        data[x][inc(y)];
                    N *= cur;
                    if (N <= 0)
                        cur = -cur;
                    else {
                        double r = dst(rnd);
                        if ((N==2 && r<e4beta) || (N==4 && r<e8beta))
                            cur = -cur;
                    }
                }
            }
            }
        }
        double measure() {
            int result = 0;
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
        std::cout << "ising with L=" << extent << ", beta=" << beta
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
#if 0
        std::cout << ising << std::endl;
        ising.sweep();
        std::cout << ising << std::endl;
#else
        for (int i=0; i<nmeas; i++) {
            for (int j=0; j<10; j++)
                ising.sweep();
            double M = ising.measure();
            write(outfile,M);
        }
        ising.save();
#endif

    };

