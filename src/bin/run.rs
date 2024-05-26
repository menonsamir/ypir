#[cfg(feature = "server")]
use ypir::scheme::run_ypir_batched;

use clap::Parser;

/// Run the YPIR scheme with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Number of items in the database
    num_items: usize,
    /// Size of each item in bits (optional, default 1), values over 8 are unsupported
    item_size_bits: Option<usize>,
    /// Number of clients (optional, default 1)
    /// to perform cross-client batching over
    num_clients: Option<usize>,
    /// Number of trials (optional, default 5)
    /// to run the YPIR scheme and average performance measurements over (a warmup trial is excluded)
    trials: Option<usize>,
    /// Verbose mode (optional)
    /// if set, run as SimplePIR
    #[clap(long, short, action)]
    is_simplepir: bool,
    /// Output report file (optional)
    /// where results will be written in JSON.
    out_report_json: Option<String>,
    /// Verbose mode (optional)
    /// if set, the program will print debug logs to stderr.
    #[clap(long, short, action)]
    verbose: bool,
}

#[cfg(feature = "server")]
fn main() {
    let args = Args::parse();
    let Args {
        num_items,
        item_size_bits,
        num_clients,
        trials,
        out_report_json,
        verbose,
        is_simplepir,
    } = args;

    if verbose {
        println!("Running in verbose mode.");
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Debug)
            .write_style(env_logger::WriteStyle::Always)
            .init();
    } else {
        env_logger::init();
    }

    let item_size_bits = item_size_bits.unwrap_or(1);
    let num_clients = num_clients.unwrap_or(1);
    let trials = trials.unwrap_or(5);

    if item_size_bits > 8 && !is_simplepir {
        panic!("Items can be at must be at most 8 bits.");
    }

    if is_simplepir {
        assert_eq!(num_clients, 1, "SimplePIR variant only supports 1 client.");
    }

    println!(
        "Running YPIR ({}) on a database of {} bits, and performing cross-client batching over {} clients. \n\
        The server performance measurement will be averaged over {} trials.",
        if is_simplepir { "w/ SimplePIR" } else { "w/ DoublePIR" },
        num_items * item_size_bits,
        num_clients,
        trials
    );

    let measurement =
        run_ypir_batched(num_items, item_size_bits, num_clients, is_simplepir, trials);
    println!(
        "Measurement completed. See the README for details on what the following fields mean."
    );
    println!("Result:");
    println!("{}", serde_json::to_string_pretty(&measurement).unwrap());

    if let Some(out_report_json) = out_report_json {
        println!("Writing report to {}", out_report_json);
        let mut file = std::fs::File::create(out_report_json).unwrap();
        serde_json::to_writer_pretty(&mut file, &measurement).unwrap();
        println!("Report written.");
    }
}

#[cfg(not(feature = "server"))]
fn main() {
    panic!("This binary is only available with the 'server' feature enabled.");
}
