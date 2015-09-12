#ifndef __ISING_H
#define __ISING_H

template <class R>
class Ising {
    private:
        Random r;
        std::vector<Line> dat[2];
        std::ofstream out;
        int measured;
    public:
        Ising(Random::result_type seed, int nmeas) :r(seed) {
            dat[0].resize(L);
            dat[1].resize(L);
            std::ofstream::openmode m;
            std::stringstream out_file_name;
            mkdir("data",0775);
            out_file_name << "data/" << std::setfill('0') << std::setw(3) << L;
            mkdir(out_file_name.str().c_str(), 0775);
            out_file_name << "/" << std::fixed << std::setprecision(20) << R::p;
            measured = getFilesize(out_file_name.str().c_str())/sizeof(floatT);
            if (load()) {
                std::cout << "# Resume "; 
                m = std::ofstream::app;
            } else {
                std::cout << "#   Init ";
                m = std::ofstream::trunc;
                if (R::p < 0.83)
                    randomize();
            }
            std::cout 
                << "L=" << L 
                << ", p=" << std::fixed << R::p
                << ", beta=" << std::fixed << -0.25*log(1-R::p) 
                << ", Nmeas=";
            if (m == std::ofstream::app)
                std::cout << std::max(nmeas-measured,0) << '/';
            std::cout << nmeas << std::endl;

            out.open(out_file_name.str().c_str(), m);
            if (out.fail()) {
                std::cerr << "Output file could not be opened!" << std::endl;
                return;
            }
            for (; keepRunning && measured<nmeas; measured++) {
                std::cout << measured << std::endl;
                for (int j=0; j<10; j++)
                    sweep();
                measure();
            }

        }
        void randomize() {
            for (int eo=0; eo<2; eo++)
                for (int x=0; x<L; x++)
                    dat[eo][x].init_random(r);
        }
        void sweep() {
            for (int eo=0; eo<2; eo++) { //sublattice, 0=even, 1=odd

                //depends on if that line represents columns 024... or 135...
                //maybe it should be !eo, but that just amounts to redefining which lattice is
                //even and which odd
                bool right = eo; 

                for (int x=0; x<L; x++) {
                    Line &c = dat[eo][x];
                    const std::vector<Line> &o = dat[1-eo];
//                     Line *o = dat[1-eo];
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

                    c1.template randomize<R>(r); //we want to flip these bits with some probability
                    c0.template randomize<R>(r);
                    c0.template randomize<R>(r); //here the probability is exp(-8beta), so we apply twice
                    c ^= c0|c1|c234;
                    right = !right;
                }
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
        ~Ising() {
            out.close();
            std::stringstream fname;
            fname << ".state";
            mkdir(fname.str().c_str(), 0755);
            fname << '/' << L;
            mkdir(fname.str().c_str(), 0755);
            fname << "/" << std::fixed << std::setprecision(20) << R::p;
            std::ofstream f(fname.str().c_str(), std::ofstream::ate);
            for (int x=0; x<L; x++) 
                f << dat[0][x] << ' ' << dat[1][x] << std::endl;
        }
        bool load() {
            std::stringstream fname;
            fname << ".state/" << L << "/" << std::fixed << std::setprecision(20) << R::p;
            std::ifstream f(fname.str().c_str());
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
#endif
