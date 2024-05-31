use rand::{rngs::OsRng, RngCore};

pub type Seed = [u8; 32]; // <ChaCha20Rng as SeedableRng>::Seed;
pub fn generate_secure_random_seed() -> Seed {
    let mut seed = [0u8; 32];
    OsRng.fill_bytes(&mut seed);
    seed
}
