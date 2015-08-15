#ifndef __ISING_H
#define __ISING_H

class Ising {
    private:
        Random r;
        Line dat[2][L];
    public:
        Ising(unsigned seed)
            :r(seed)
        { }
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

                    c1.randomize(r); //we want to flip these bits with some probability
                    c0.randomize(r);
                    c0.randomize(r); //here the probability is exp(-8beta), so we apply twice
                    c ^= c0|c1|c234;
                    right = !right;
                }
            }
        }
        double magnetization() {
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
#endif
