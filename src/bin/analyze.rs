use std::env::args;
use std::path::Path;
use std::fs::File;

use memmap::MmapOptions;

// Compute sum and sum of squares over slice
fn sum_and_sqrsum(values: &[f32]) -> [f32; 2] {
    let mut result = [0.0, 0.0];
    for value in values {
        result[0] += value;
        result[1] += value*value;
    };
    result
}
fn sum(values: &[[f32; 2]]) -> [f32; 2] {
    let mut result = [0.0, 0.0];
    for value in values {
        result[0] += value[0];
        result[1] += value[1];
    };
    result
}

fn main() {
    // Parse args
    let mut args = args();
    let _ = args.next();
    let blockcount: usize = args.next().unwrap().parse().unwrap();
    let therm: f32 = args.next().unwrap().parse().unwrap();
    let path = &args.next().unwrap();

    // Parse path, data/<extent>/<sig>
    let mut parts = path.split("/");
    let _ = parts.next();
    let extent: u32 = parts.next().unwrap().parse().unwrap();
    let sig = parts.next().unwrap();

    // Convert signature to probability
    let mut p: f32 = 0.0;
    let mut weight = 1.0;
    for b in sig.bytes() {
        weight /= 2.0;
        if b == b'1' {
            p += weight;
        }
    }
    let beta = -0.25 * (1.0-p).ln();

    // Read values
    let size: usize = Path::new(path).metadata().unwrap().len() as usize /4;
    let file = File::open(&path).unwrap();
    let mmap = unsafe { MmapOptions::new().map(&file).unwrap() };

    let skip = (therm*size as f32) as usize;
    let blocksize = (size-skip) / blockcount;
    let count = blocksize * blockcount;

    // Jackknife algorithm. Estimate average and variance on the full data (after thermalization
    // skip) as well as on the data that is obtained by leaving out one block. The pseudo-jackknife
    // values are given by subtracting these after multiplying them with the respective sample
    // size.
    let mut partial = Vec::with_capacity(blockcount);
    for i in 0..blockcount {
        let mut sum = 0.0;
        let mut sumsqr = 0.0;
        for j in 0..blocksize {
            let offset = 4*skip + 4*i*blocksize + 4*j;
            let value = f32::from_ne_bytes(mmap[offset..offset+4].try_into().unwrap());
            sum += value;
            sumsqr += value*value;
        }
        partial.push([sum, sumsqr]);
    }
    let total = sum(&partial);
    let std_avr = total[0]/count as f32;
    let std_var = total[1] / count as f32 - std_avr*std_avr;

    let mut jackknife = Vec::with_capacity(blockcount);
    for i in 0..blockcount {
        let mag = (total[0] - partial[i][0]) / (count - blocksize) as f32;
        let var = (total[1] - partial[i][1]) / (count - blocksize) as f32 - mag*mag;
        jackknife.push(blockcount as f32 * std_var - (blockcount as f32 - 1.0) * var);
    }

    // Compute susceptibility and its variance from average and variance of these
    let jksum = sum_and_sqrsum(&jackknife);
    let scale = (extent*extent) as f32 * beta;
    let sus_avr = jksum[0] / blockcount as f32;
    let sus_var = (jksum[1] / blockcount as f32 - sus_avr*sus_avr) / (blockcount as f32 - 1.0);
    eprintln!("L = {extent}  p = {p:<10}  beta = {beta:<12}  N = {count:>10}");

    println!("{extent} {p} {beta} {} {}", scale*sus_avr, scale*sus_var.sqrt());
}
