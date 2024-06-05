FROM rust:1-bookworm as builder
RUN apt-get update && apt-get install -y clang-14 clang++-14
WORKDIR /usr/src/ypir
COPY . .
RUN CC=clang-14 CXX=clang++-14 cargo install --features explicit_avx512 --bin run --path .

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y clang-14 clang++-14 libssl-dev pkg-config && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/cargo/bin/run /usr/local/bin/run
ENV RUST_LOG=debug
ENTRYPOINT ["run"]