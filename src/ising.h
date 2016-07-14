#ifndef __ISING_H
#define __ISING_H

const int SWEEP_PER_MEAS = 0x100;
struct Lattice {
    Random r;
    std::string sig;
    std::vector<Line> dat[2];
    std::ofstream out;
    int measured;

    Lattice(std::string _sig, Random::result_type seed, int nmeas) 
        :r(seed),sig(_sig)
    {
        dat[0].resize(L);
        dat[1].resize(L);

        // set p from sig
        floatT p = 0;
        floatT val = 1;
        for (int i=0; i<sig.size(); i++) {
            val /= 2;
            if (sig[i] == '1')
                p += val;
        }

        std::ofstream::openmode m;
        std::stringstream out_file_name;
        mkdir("data",0775);
        out_file_name << "data/" << std::setfill('0') << std::setw(3) << L;
        mkdir(out_file_name.str().c_str(), 0775);
        out_file_name << "/" << sig;
        measured = getFilesize(out_file_name.str().c_str())/sizeof(floatT);
        if (load()) {
            std::cout << "# Resume "; 
            m = std::ofstream::app;
        } else {
            std::cout << "#   Init ";
            m = std::ofstream::trunc;
            if (p < 0.828)
                for (int eo=0; eo<2; eo++)
                    for (int x=0; x<L; x++)
                        dat[eo][x].init_random(r);
        }
        // output information
        std::cout 
            << "L=" << L 
            << ", p=" << std::fixed << p
            << ", beta=" << std::fixed << (-0.25*log(1-p))
            << ", Nmeas=";
        if (m == std::ofstream::app)
            std::cout << std::max(nmeas-measured,0) << '/';
        std::cout << nmeas << std::endl;

        out.open(out_file_name.str().c_str(), m);
        if (out.fail()) {
            std::cerr << "Output file could not be opened!" << std::endl;
            return;
        }
    }
    void measure() {
        int count = 0;
        for (int eo=0; eo<2; eo++)
            for (int x=0; x<L; x++)
                count += dat[eo][x].count();
        float M = fabs(float(count *2)/(L*L)-1);
        write(out,M);
    }
    ~Lattice() {
        out.close();
        std::stringstream fname;
        fname << ".state";
        mkdir(fname.str().c_str(), 0755);
        fname << '/' << L;
        mkdir(fname.str().c_str(), 0755);
        fname << "/" << sig;
        std::ofstream f(fname.str().c_str(), std::ofstream::ate);
        for (int eo=0; eo<2; eo++)
            for (int x=0; x<L; x++)
                dat[eo][x].write(f);
    }
    bool load() {
        std::stringstream fname;
        fname << ".state/" << L << "/" << sig;
        std::ifstream f(fname.str().c_str());
        if (f.fail()) {
            return false;
        }
        for (int eo=0; eo<2; eo++)
            for (int x=0; x<L; x++)
                dat[eo][x].read(f);
        return !f.fail();
    }

};

template <class R>
class Ising {
    private:
        Lattice lat;
    public:
        Ising(std::string sig, Random::result_type seed, int nmeas) : lat(sig,seed,nmeas) {
            for (; keepRunning && lat.measured<nmeas; lat.measured++) {
                for (int j=0; j<SWEEP_PER_MEAS; j++)
                    sweep();
                lat.measure();
            }
        }
        void sweep() {
            for (int eo=0; eo<2; eo++) { //sublattice, 0=even, 1=odd

                //depends on if that line represents columns 024... or 135...
                //maybe it should be !eo, but that just amounts to redefining which lattice is
                //even and which odd
                bool right = eo; 

                for (int x=0; x<L; x++) {
                    Line &c = lat.dat[eo][x];
                    const std::vector<Line> &o = lat.dat[1-eo];
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

                    c1.template randomize<R>(lat.r); //we want to flip these bits with some probability
                    c0.template randomize<R>(lat.r);
                    c0.template randomize<R>(lat.r); //here the probability is exp(-8beta), so we apply twice
                    c ^= c0|c1|c234;
                    right = !right;
                }
            }
        }
};
#endif
