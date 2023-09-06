mod rnd;
mod generated;

use std::io::{Write, Error, Read};
use std::env::args;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::fs::{File, create_dir_all, OpenOptions};
use std::path::Path;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use crate::generated::{rand_p,N,SIG};

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
        let mut shift_right = self.odd;
        for row in 0..L {
            let nb = [
                other.row[(row+L-1)%L],
                other.row[row].shift(shift_right),
                other.row[(row+1)%L],
                other.row[row],
            ];
            self.row[row].update(r, nb);
            shift_right = !shift_right;
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

    fn store(&self, file: &mut File) {
        for i in 0..L {
            for j in 0..N {
                file.write_all(&self.row[i].c[j].to_ne_bytes()).unwrap();
            }
        }
    }

    fn load(&mut self, file: &mut File) {
        for i in 0..L {
            for j in 0..N {
                let mut target = [0u8; 16];
                assert!(file.read(&mut target).unwrap() == 16);
                self.row[i].c[j] = u128::from_ne_bytes(target);
            }
        }
    }
}

struct Ising {
    even: HalfLattice,
    odd: HalfLattice,
    rng: SmallRng,
    terminating: Arc<AtomicBool>,
    measured: u64,
    file: File,
}

impl Ising {
    fn new() -> Result<Ising, Error> {
        let terminating = Arc::new(AtomicBool::new(false));
        signal_hook::flag::register_conditional_shutdown(
            signal_hook::consts::SIGINT,
            1,
            Arc::clone(&terminating)
        )?;
        let _ = signal_hook::flag::register(
            signal_hook::consts::SIGINT,
            Arc::clone(&terminating)
        );
        let path = Path::new("data").join(format!("{}", L)).join(SIG);
        let mut opt = OpenOptions::new();
        let (measured, opt) = match path.metadata() {
            Err(_) => (0, opt.write(true).create(true)),
            Ok(md) => (md.len()/4, opt.append(true)),
        };
        Ok(Ising {
            even: HalfLattice::new(false),
            odd: HalfLattice::new(true),
            rng: SmallRng::from_entropy(),
            terminating: terminating,
            measured: measured,
            file: opt.open(&path).unwrap(),
        })
    }

    fn sweep(&mut self) {
        for _ in 0..16 {
            self.even.update(&mut self.rng, &self.odd);
            self.odd.update(&mut self.rng, &self.even);
        }
        let mag = self.even.mag() + self.odd.mag();
        let mag = (((2*mag) as f32)/((L*L) as f32)-1.0).abs();
        self.file.write_all(&mag.to_ne_bytes()).unwrap();
    }

    fn run(&mut self, nmeas: u64) {
        let path = Path::new(".state").join(format!("{}", L));
        create_dir_all(&path).unwrap();
        let path = path.join(SIG);
        if path.exists() {
            let mut file = File::open(&path).unwrap();
            self.even.load(&mut file);
            self.odd.load(&mut file);
        }

        for _ in self.measured..nmeas {
            self.sweep();
            if self.terminating.load(Ordering::Relaxed) {
                break;
            }
        }

        let mut file = OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(&path).unwrap();
        self.even.store(&mut file);
        self.odd.store(&mut file);
    }
}

fn main() -> Result<(), Error> {
    let mut args = args();
    let _ = args.next();
    let nmeas: u64= args.next().unwrap().parse().unwrap();
    let mut ising = Ising::new()?;
    ising.run(nmeas);
    Ok(())
}
