FROM rust:1.67 as builder
WORKDIR /usr/src/ypir
COPY . .
RUN cargo install --bin run --path .

FROM debian:bullseye-slim
RUN apt-get update && apt-get install -y libssl-dev pkg-config && rm -rf /var/lib/apt/lists/*
COPY --from=builder /usr/local/cargo/bin/run /usr/local/bin/run
ENV RUST_LOG=debug
ENTRYPOINT ["run"]