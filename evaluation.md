# Evaluation

This is a step-by-step guide on how to run YPIR to evaluate its key results.

## Getting Started
We recommend running these instructions in an Ubuntu (or similar Linux) environment.

1. Download and install [Docker](https://docs.docker.com/engine/install/ubuntu/)
2. Run `sudo docker run --security-opt seccomp:unconfined --cpus=1 ghcr.io/menonsamir/ypir 2147483648 1`
  - The `seccomp:unconfined` option disables container sandboxing that Docker runs by default which can degrade performance
  - The `--cpus=1` option runs a single-threaded container (this is not crucial - the implementation will always use only 1 thread)
  - The `2147483648 1` arguments indicate running on a database of 2^31 items, each 1 bit (256 MB)
    - The number of items must be a power of two, and for items larger than 8 bits, you must switch to YPIR+SP by passing in `--is-simplepir` as an additional argument


## Building
To build from source, just run `cargo build --profile release-with-debug`. To run a basic test, run `cargo run --profile release-with-debug --bin run 2147483648 --verbose`
