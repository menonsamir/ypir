use clap::Parser;
use reqwest::blocking::Client;
use spiral_rs::params::*;
use std::error::Error;
use ypir::bits::u64s_to_contiguous_bytes;
use ypir::client::YPIRClient;
use ypir::params::*;
use ypir::scheme::*;
use ypir::serialize::ToBytes;
use ypir::server::*;

/// Run the YPIR server with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Row to fetch
    target_row: usize,
    /// Number of items in the database
    num_items: usize,
    /// Size of each item in bits (optional, default 1), values over 8 are unsupported
    item_size_bits: usize,
    /// If set, run using SimplePIR instead of Double
    #[clap(long, short, action)]
    is_simplepir: bool,
    /// Port
    #[clap(long, short, default_value = "8080")]
    port: u16,
    /// Verbose mode (optional)
    /// if set, the program will print debug logs to stderr
    #[clap(long, short, action)]
    verbose: bool,
}

fn main() {
    let args = Args::parse();
    let Args {
        target_row,
        num_items,
        item_size_bits,
        verbose,
        is_simplepir,
        port,
    } = args;

    if !is_simplepir {
        panic!("Must use YPIR-SP for now.");
    }

    if verbose {
        println!("Running in verbose mode.");
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Debug)
            .write_style(env_logger::WriteStyle::Always)
            .init();
    } else {
        env_logger::init();
    }

    let client = YPIRClient::from_db_sz(num_items as u64, item_size_bits as u64, is_simplepir);
    assert!(target_row < client.params().db_rows());

    let (query, client_seed) = client.generate_query_simplepir(target_row);
    let query_bytes = query.to_bytes();
    let response_data: Vec<u8> =
        make_http_request(&format!("http://localhost:{}/query", port), query_bytes).unwrap();
    let result = client.decode_response_simplepir(client_seed, &response_data);

    let result_bytes = u64s_to_contiguous_bytes(&result, client.params().pt_modulus_bits());

    println!(
        "Result: {:?}..{:?}",
        &result_bytes[..32],
        &result_bytes[result_bytes.len() - 32..]
    );
}

fn make_http_request(url: &str, query_bytes: Vec<u8>) -> Result<Vec<u8>, Box<dyn Error>> {
    let client = Client::new();
    let response = client.post(url).body(query_bytes).send()?.bytes()?;

    Ok(response.to_vec())
}
