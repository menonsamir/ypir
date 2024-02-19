# YPIR

This is an implementation of the YPIR scheme for single-server private information retrieval,
introduced in ["YPIR: High-Throughput Single-Server PIR with Silent Preprocessing"](https://eprint.iacr.org/2024/270).
This is joint work with [David Wu](https://www.cs.utexas.edu/~dwu4/).

## Running

To build and run this code:
1. Ensure you are running on Ubuntu, and that AVX-512 is available on the CPU (you can run `lscpu` and look for the `avx512f` flag).
Our benchmarks were collected using the AWS `r6i.16xlarge` instance type, which has all necessary CPU features.
2. Run `sudo apt-get update && sudo apt-get install -y build-essential`.
2. [Install Rust using rustup](https://www.rust-lang.org/tools/install) using `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`.
  - Select `1) Proceed with installation (default)` when prompted
  - After installation, configure the current shell as instructed by running `source "$HOME/.cargo/env"`
3. Run `git clone https://github.com/menonsamir/ypir.git` and `cd ypir`.
4. Run `cargo run --release -- 1073741824` to run YPIR on a random database consisting of 1073741824 bits (~134 MB).
The first time you run this command, Cargo will download and install the necessary libraries to build the code (~2 minutes);
later calls will not take as long. Stability warnings can be safely ignored. 
See below for details on how to interpret the measurements.

We have tested the above steps on a fresh AWS `r6i.16xlarge` Ubuntu 22.04 instance and confirmed they work.

### Options
To pass arguments, make sure to run `cargo run --release -- <ARGS>` (the ` -- ` is important).
Passing `--verbose` or setting the environment variable `RUST_LOG=debug`
will enable detailed logging. All PIR results are checked for correctness.
The full command-line parameters are as follows:

```
Usage: cargo run --release -- [OPTIONS] <NUM_ITEMS> [ITEM_SIZE_BITS] [NUM_CLIENTS] [TRIALS] [OUT_REPORT_JSON]

Arguments:
  <NUM_ITEMS>        Number of items in the database
  [ITEM_SIZE_BITS]   Size of each item in bits (optional, default 1), values over 8 are unsupported
  [NUM_CLIENTS]      Number of clients (optional, default 1) to perform cross-client batching over
  [TRIALS]           Number of trials (optional, default 5) to run the YPIR scheme 
                     and average performance measurements over (with one additional warmup trial excluded)
  [OUT_REPORT_JSON]  Output report file (optional) where results will be written in JSON

Options:
  -v, --verbose  Verbose mode (optional) if set, the program will print debug logs to stderr
  -h, --help     Print help
  -V, --version  Print version
```

### Interpreting measurements
This is an annotated version of the output, detailing what each measurement means:
```js
{
  "offline": {
    // Bytes uploaded by the client in the offline phase
    "uploadBytes": 0,

    // Bytes downloaded by the client in the offline phase
    "downloadBytes": 0,
    
    // Server computation time, in milliseconds, in the offline phase. 
    // Includes any precomputation that must be performed on the plaintext database.
    "serverTimeMs": 3965,
    
    // Not used.
    "clientTimeMs": 0,
    
    // Time spent precomputing just the SimplePIR hint.
    "simplepirPrepTimeMs": 2539,

    // Bytes that the client *would* have to download, in the offline phase,
    // if they were performing SimplePIR (rather than YPIR) 
    // using this implementation (SimplePIR* in the paper).
    "simplepirHintBytes": 29360128,

    // Similarly, bytes that the client *would* have to download, 
    // in the offline phase DoublePIR (DoublePIR* in the paper).
    "doublepirHintBytes": 14680064
  },
  "online": {
    // Bytes uploaded by a single client in the online phase.
    "uploadBytes": 604160,
    
    // Bytes downloaded by a single client in the online phase.
    "downloadBytes": 12288,

    // Bytes that the client *would* have to download, in the online phase,
    // if they were performing SimplePIR (SimplePIR* in the paper).
    "simplepirRespBytes": 28672,

    // Bytes that the client *would* have to download, in the online phase,
    // if they were performing DoublePIR (DoublePIR* in the paper).
    "doublepirRespBytes": 12288,

    // Server computation time, in milliseconds, in the online phase.
    // This is the average time over 5 trials, after a warmup trial.
    "serverTimeMs": 402,

    // Time that the client took to generate the query.
    "clientQueryGenTimeMs": 530,

    // Time that the client took to decode the response (may round down to 0ms).
    "clientDecodeTimeMs": 0,

    // Time spent in the first pass of YPIR (the 'SimplePIR' phase)
    "firstPassTimeMs": 9,

    // Time spent in the second pass of YPIR (the 'DoublePIR' phase)
    "secondPassTimeMs": 3,

    // Time spent performing LWE-to-RLWE conversion
    "ringPackingTimeMs": 387,

    // Not used.
    "sqrtNBytes": 8192,

    // The full set of measured server computation times.
    "allServerTimesMs": [
      401,
      403,
      402,
      401,
      401
    ],
    // The standard deviation of the measured server computation times.
    "stdDevServerTimeMs": 0.8
  }
}
```

### Acknowledgements

YPIR is based on [DoublePIR](https://eprint.iacr.org/2022/949), and this implementation
uses matrix-vector multiplication routines based on the ones in [ahenzinger/simplepir](https://github.com/ahenzinger/simplepir).
We also use the [menonsamir/spiral-rs](https://github.com/menonsamir/spiral-rs) library for Spiral to handle RLWE ciphertexts.

### Citing

```
@misc{cryptoeprint:2024/270,
  author = {Samir Jordan Menon and David J. Wu},
  title = {YPIR: High-Throughput Single-Server PIR with Silent Preprocessing},
  howpublished = {Cryptology ePrint Archive, Paper 2024/270},
  year = {2024},
  note = {\url{https://eprint.iacr.org/2024/270}},
  url = {https://eprint.iacr.org/2024/270}
}
```

(to be updated)
