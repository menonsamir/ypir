#![cfg_attr(feature = "explicit_avx512", feature(stdarch_x86_avx512))]

pub mod bits;
pub mod client;
pub mod constants;
pub mod convolution;
pub mod lwe;
pub mod measurement;
pub mod modulus_switch;
pub mod noise_analysis;
pub mod params;
pub mod seed;
pub mod serialize;
pub mod transpose;
pub mod util;

#[cfg(feature = "server")]
pub mod kernel;
#[cfg(feature = "server")]
pub mod matmul;
#[cfg(feature = "server")]
pub mod packing;
#[cfg(feature = "server")]
pub mod scheme;
#[cfg(feature = "server")]
pub mod server;

#[cfg(feature = "test_data")]
pub mod data;
