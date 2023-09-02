use rand::RngCore;
pub use rand::rngs::SmallRng;

/// Generate random 128 bit number from SmallRng
pub fn rng(r: &mut SmallRng) -> u128 {
    ((r.next_u64() as u128) << 64) | (r.next_u64() as u128)
}
