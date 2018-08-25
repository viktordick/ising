#ifndef __ISING_H
#define __ISING_H

const int SWEEP_PER_MEAS = 0x10;

/** The 2D Ising lattice
 * It contains a random number generator and an even and an odd sublattice,
 * which each consist of Lines of data.
 *
 * All code that is not dependent on the flip probability (i.e., beta), is
 * collected in the struct Lattice. The part that can benefit by the compiler
 * already knowing which probability is going to be used is put into the
 * template Ising.
 * This separation is done to reduce compile time without impacting run time.
 * All code in Ising is built separately for each given probability, with all
 * optimizations the compiler can perform. Code in Lattice is only built once.
 */
struct Lattice {
    Random r;
    std::string sig;
    std::vector<Line> dat[2];
    std::ofstream out;
    int measured;

    /** Initialize Lattice. 
     * Load last state from file if found, otherwise set to either zero or
     * random, depending on the signature and therefore on beta (hot or cold
     * start).
     */
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
    /** Measure magnetization as float value and write (binary!) to output file.
     */
    void measure() {
        int count = 0;
        for (int eo=0; eo<2; eo++)
            for (int x=0; x<L; x++)
                count += dat[eo][x].count();
        float M = fabs(float(count *2)/(L*L)-1);
        write(out,M);
    }

    /** Destructor. Store the current state into .state/<L>/<sig>
     */
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
    /** Load the state that was last used when running this
     * signature/prob./beta value on this lattice size.
     */
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

/** The Ising template.
 * As a template argument, it takes a Random struct that provides the method
 * get(). This method should return a bit pattern where each bit is set with a
 * given probability.
 *
 * The code collected here can benefit from the compiler already knowing the
 * probability and therefore which bitwise and/or operations are performed on
 * the raw random numbers by R::get().
 */
template <class R>
class Ising {
    private:
        Lattice lat;
    public:
        /** Initialize and run */
        Ising(std::string sig, Random::result_type seed, int nmeas) : lat(sig,seed,nmeas) {
            for (; keepRunning && lat.measured<nmeas; lat.measured++) {
                for (int j=0; j<SWEEP_PER_MEAS; j++)
                    sweep();
                lat.measure();
            }
        }
        /** Do one sweep of the lattice.
         * For the even and odd sublattices, iterate through Lines and
         * 1) collect Lines that represent neighbors, each XORed with the
         *    current Line so each 1 means an antiparallel neighbor
         * 2) Create Lines c234, c1 and c0, which have a 1 in each bit if
         *    - at least two neighbors are antiparallel
         *    - exactly one neighbor is antiparallel
         *    - no neighbor is antiparallel
         * 3) Flip the ones of c0 and c1 with some probabilities depending on
         *    beta (exp(-4beta) for c1, exp(-8beta) for c0). 
         * 4) Join c0, c1 and c234 by OR, so we obtain a Line where the bits
         *    are set if the corresponding bit in the original Line should be
         *    flipped (either there are more than two antiparallel neighbors or
         *    there is one antiparallel neighbor and an event with probability
         *    exp(-4beta) occured or there is no antiparallel neighbor and an
         *    event with probability exp(-8beta) occured).
         * 5) Flip all bits of the current Line where the result of 4) has a
         *    set bit.
         */
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
