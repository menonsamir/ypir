use std::time::Instant;

use actix_cors::Cors;
use actix_web::HttpServer;
use actix_web::{get, post, web, App};
use clap::Parser;
use spiral_rs::params::*;
use ypir::bits::u64s_to_contiguous_bytes;
use ypir::params::*;
use ypir::scheme::*;
use ypir::serialize::FilePtIter;
use ypir::server::*;

/// Run the YPIR server with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Number of items in the database
    num_items: usize,
    /// Size of each item in bits (optional, default 1), values over 8 are unsupported
    item_size_bits: Option<usize>,
    /// If set, run using SimplePIR instead of Double
    #[clap(long, short, action)]
    is_simplepir: bool,
    /// Port
    #[clap(long, short, default_value = "8080")]
    port: u16,
    /// Read the database from an input file
    #[clap(long, short)]
    inp_file: Option<String>,
    /// Use a random database
    #[clap(long, short, action)]
    random: bool,
    /// Verbose mode (optional)
    /// if set, the program will print debug logs to stderr
    #[clap(long, short, action)]
    verbose: bool,
}

type T = u16; // TODO: make this generic

#[derive(Clone)]
struct ServerState {
    params: &'static Params,
    server: YServer<'static, T>,
    offline_values: OfflinePrecomputedValues<'static>,
}

#[post("/query")]
async fn query(
    body: web::Bytes,
    data: web::Data<ServerState>,
) -> Result<Vec<u8>, actix_web::error::Error> {
    let req_body = body.to_vec();
    let response = data
        .server
        .perform_full_online_computation_simplepir(&data.offline_values, &req_body);
    Ok(response)
}

#[get("/")]
async fn index(data: web::Data<ServerState>) -> String {
    format!("Hello {}!", data.params.poly_len)
}

#[get("/info")]
async fn info() -> String {
    format!("Info!")
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    use actix_web::web::Data;

    let args = Args::parse();
    let Args {
        num_items,
        item_size_bits,
        inp_file,
        random,
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

    let item_size_bits = item_size_bits.unwrap_or(16384 * 8);

    if item_size_bits > 8 && !is_simplepir {
        panic!("Items can be at must be at most 8 bits.");
    } else if is_simplepir && item_size_bits < 2048 {
        panic!("YPIR-SP requires items to be at least 2048 bits.");
    }

    println!(
        "Starting a YPIR ({}) server on a database of {} bits.",
        if is_simplepir {
            "w/ SimplePIR"
        } else {
            "w/ DoublePIR"
        },
        num_items * item_size_bits,
    );

    let params = if is_simplepir {
        params_for_scenario_simplepir(num_items as u64, item_size_bits as u64)
    } else {
        params_for_scenario(num_items as u64, item_size_bits as u64)
    };
    let pt_modulus = params.pt_modulus;
    let leaked_params = Box::leak(Box::new(params));

    let server = if random {
        let pt_iter = std::iter::repeat_with(|| (u16::sample() as u64 % pt_modulus) as u16);
        YServer::<u16>::new(leaked_params, pt_iter, true, false, true)
    } else {
        assert!(inp_file.is_some());
        let inp_file = inp_file.unwrap();
        let pt_bits = (pt_modulus as f64).log2().ceil() as usize;
        let pt_iter = FilePtIter::from_file(
            &inp_file,
            item_size_bits / 8,
            leaked_params.db_cols_simplepir(),
            pt_bits,
        );
        YServer::<u16>::new(leaked_params, pt_iter, true, false, true)
    };
    println!("Performing precomputation...");
    let now = Instant::now();
    let offline_values = server.perform_offline_precomputation_simplepir(None);
    println!("Done. ({} s)", now.elapsed().as_secs());

    let corr_result_item_1 = server
        .get_row(1)
        .iter()
        .map(|x| x.to_u64())
        .collect::<Vec<_>>();
    let ci_bytes = u64s_to_contiguous_bytes(&corr_result_item_1, leaked_params.pt_modulus_bits());
    println!("item_1: {:?}", &ci_bytes[..32]);

    let state = ServerState {
        params: leaked_params,
        server,
        offline_values,
    };

    println!("Listening on http://127.0.0.1:{}", port);
    HttpServer::new(move || {
        App::new()
            .wrap(Cors::permissive())
            .app_data(Data::new(state.clone()))
            .app_data(web::PayloadConfig::new(1usize << 32))
            .service(index)
            .service(query)
            .service(info)
    })
    .workers(1)
    .bind(("127.0.0.1", port))
    .unwrap()
    .run()
    .await
}

#[cfg(not(feature = "http_server"))]
fn main() {
    panic!("This binary is only available with the 'server' feature enabled.");
}
