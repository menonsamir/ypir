use clap::Parser;
use reqwest::blocking::Client;
use sha1::{Digest, Sha1};
use std::error::Error;
use std::time::Instant;
use ypir::client::YPIRClient;
use ypir::params::*;
use ypir::serialize::ToBytes;

/// Run the YPIR server with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Number of items in the database
    num_items: usize,
    /// Size of each item in bits (optional, default 1), values over 8 are unsupported
    item_size_bits: usize,
    /// Row to fetch
    #[clap(long)]
    target_row: Option<usize>,
    /// Item to check for inclusion
    #[clap(long)]
    target_item: Option<String>,
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

const SHA1_HASH_BYTES: usize = 20;

fn main() {
    let args = Args::parse();
    let Args {
        target_row,
        target_item,
        num_items,
        item_size_bits,
        verbose,
        is_simplepir,
        port,
    } = args;

    let log2_num_items = (num_items as f64).log2().ceil() as usize;

    let (target_row, item_hash) = if let Some(target_row) = target_row {
        (target_row, None)
    } else {
        let target_item = target_item.expect("Must provide either target_row or target_item");
        let mut hasher: sha1::digest::core_api::CoreWrapper<sha1::Sha1Core> = Sha1::new();
        hasher.update(target_item.as_bytes());
        let item_hash = hasher.finalize();

        let top_idx = u32::from_be_bytes(item_hash[0..4].try_into().unwrap());
        let bucket = top_idx >> (32 - log2_num_items);
        println!("Bucket: {}", bucket);
        (bucket as usize, Some(item_hash))
    };

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
    let now = Instant::now();
    let response_data: Vec<u8> =
        make_http_request(&format!("http://localhost:{}/query", port), query_bytes).unwrap();
    println!("Query time: {:?}", now.elapsed().as_secs_f64());
    let result = client.decode_response_simplepir(client_seed, &response_data);

    println!(
        "Result: {:?}..{:?}",
        &result[..32],
        &result[result.len() - 32..]
    );

    let mut end = result.len() - 1;
    while end > 0 && result[end] == 0 {
        end -= 1;
    }
    println!("Result[{}-36..]: {:?}", end, &result[end - 36..end + 36]);

    println!("Result[32940..]: {:?}", &result[32940..32940 + 36],);

    if let Some(item_hash) = item_hash {
        // search for item_hash in result
        // top floor(log2(num_items) / 8) bytes of each hash in result are omitted
        let omitted_bytes = log2_num_items / 8;
        let hash_bytes = SHA1_HASH_BYTES - omitted_bytes;
        println!("Hash bytes: {}", hash_bytes);
        let looking_for = &item_hash[omitted_bytes..];
        println!("looking_for: {:?}", looking_for);

        let mut found = false;
        for chunk in result.chunks_exact(hash_bytes) {
            if chunk == looking_for {
                found = true;
                break;
            }
        }

        if found {
            println!("Item found!");
        } else {
            println!("Item not found.");
        }
    }
}

fn make_http_request(url: &str, query_bytes: Vec<u8>) -> Result<Vec<u8>, Box<dyn Error>> {
    let client = Client::new();
    let response = client.post(url).body(query_bytes).send()?.bytes()?;

    Ok(response.to_vec())
}
