mod rnd;
mod generated;

use std::io::{stdout,Write};
use std::env::args;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::generated::{rand_p,N};

const L: usize = N*256;

#[derive(Clone, Copy)]
struct Line {
    c: [u128; N],
}

impl Line {
    fn new() -> Self {
        Line {c: [0; N]}
    }
    fn shift(&self, right: bool) -> Self {
        let mut c = [0u128; N];
        if right {
            c[0] = (self.c[0] >> 1) | (self.c[N-1]<<127);
            for i in 1..N {
                c[i] = (self.c[i]>>1) | (self.c[i-1]<<127);
            }
        } else {
            c[N-1] = (self.c[0] << 1) | (self.c[N-1]>>127);
            for i in 0..N-1 {
                c[i] = (self.c[i] << 1) | (self.c[i+1]>>127);
            }
        }
        Line {c: c}
    }

    /// Update line according to randomness and neighboring lines
    fn update(&mut self, r: &mut SmallRng, nb: [Line; 4]) {
        for i in 0..N {
            // Count how many neighbors are antiparallel
            let mut n = [0u128; 4];
            for j in 0..4 {
                n[j] = self.c[i] ^ nb[j].c[i];
            }
            let mut c0 = !(n[0]|n[1]|n[2]|n[3]);
            let mut c1 = n[0] & !n[1] & !n[2] & !n[3] |
                    !n[0] &  n[1] & !n[2] & !n[3] |
                    !n[0] & !n[1] &  n[2] & !n[3] |
                    !n[0] & !n[1] & !n[2] &  n[3];
            let c234 = !c0 & !c1;
            // Keep bits in c0 with probability exp(-8beta) and c1 with exp(-4beta)
            c0 &= rand_p(r) & rand_p(r);
            c1 &= rand_p(r);
            // Flip each bit if it is in one of the masks
            self.c[i] ^= c0|c1|c234;
        }
    }
}

struct HalfLattice {
    odd: bool,
    row: [Line; L],
}

impl HalfLattice {
    fn new(odd: bool) -> Self {
        let line = Line::new();
        Self {odd: odd, row: [line; L]}
    }

    fn update(&mut self, r: &mut SmallRng, other: &HalfLattice) {
        for row in 0..L {
            let nb = [
                other.row[(row+L-1)%L],
                other.row[row].shift(self.odd),
                other.row[(row+1)%L],
                other.row[row],
            ];
            self.row[row].update(r, nb);
        }
    }

    fn mag(&self) -> u32 {
        let mut result = 0;
        for i in 0..L {
            for j in 0..N {
                result += self.row[i].c[j].count_ones()
            }
        }
        result
    }
}

fn main() {
    let mut args = args();
    let _ = args.next();
    let nmeas: usize = args.next().unwrap().parse().unwrap();
    let mut rng = SmallRng::from_entropy();
    let mut even = HalfLattice::new(false);
    let mut odd = HalfLattice::new(true);
    let mut stdout = stdout();
    for _ in 0..nmeas {
        for _ in 0..16 {
            even.update(&mut rng, &odd);
            odd.update(&mut rng, &even);
        }
        let mag = even.mag() + odd.mag();
        //println!("{mag}");
        stdout.write_all(&mag.to_ne_bytes()).unwrap();
    }
}
