#ifndef __LINE_H
#define __LINE_H

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
        static const T FullLine = T((~T(0))<<(N*B-LH))>>(N*B-LH);
        //bit for last column is least significant, i.e. value '1'
        T dat[N];
    public:
        void init_random(Random &r) {
            for (int i=0; i<N; i++)
                dat[i] = r.uniform();
            dat[0] &= FullLine;
        }
        //only keep 'up' bits with probability exp(-beta dE) 
        void randomize(Random &r) {
            for (int i=0; i<N; i++)
                dat[i] &= ~r.e4beta();
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

#endif
