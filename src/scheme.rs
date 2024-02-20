use std::time::Instant;

use log::debug;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use serde_json::Value;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{arith::*, client::*, params::*, poly::*};

use super::{client::*, measurement::*, server::*};

pub const STATIC_PUBLIC_SEED: [u8; 32] = [0u8; 32];
pub const SEED_0: u8 = 0;
pub const SEED_1: u8 = 1;

static DEFAULT_MODULI: [u64; 2] = [268369921u64, 249561089u64];
const DEF_MOD_STR: &str = "[\"268369921\", \"249561089\"]";

fn ext_params_from_json(json_str: &str) -> Params {
    let v: Value = serde_json::from_str(json_str).unwrap();

    let n = v["n"].as_u64().unwrap() as usize;
    let db_dim_1 = v["nu_1"].as_u64().unwrap() as usize;
    let db_dim_2 = v["nu_2"].as_u64().unwrap() as usize;
    let instances = v["instances"].as_u64().unwrap_or(1) as usize;
    let p = v["p"].as_u64().unwrap();
    let q2_bits = u64::max(v["q2_bits"].as_u64().unwrap(), MIN_Q2_BITS);
    let t_gsw = v["t_gsw"].as_u64().unwrap() as usize;
    let t_conv = v["t_conv"].as_u64().unwrap() as usize;
    let t_exp_left = v["t_exp_left"].as_u64().unwrap() as usize;
    let t_exp_right = v["t_exp_right"].as_u64().unwrap() as usize;
    let do_expansion = v.get("direct_upload").is_none();

    let mut db_item_size = v["db_item_size"].as_u64().unwrap_or(0) as usize;
    if db_item_size == 0 {
        db_item_size = instances * n * n;
        db_item_size = db_item_size * 2048 * log2_ceil(p) as usize / 8;
    }

    let version = v["version"].as_u64().unwrap_or(0) as usize;

    let poly_len = v["poly_len"].as_u64().unwrap_or(2048) as usize;
    let moduli = v["moduli"]
        .as_array()
        .map(|x| {
            x.as_slice()
                .iter()
                .map(|y| {
                    y.as_u64()
                        .unwrap_or_else(|| y.as_str().unwrap().parse().unwrap())
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or(DEFAULT_MODULI.to_vec());
    let noise_width = v["noise_width"].as_f64().unwrap_or(6.4);

    Params::init(
        poly_len,
        &moduli,
        noise_width,
        n,
        p,
        q2_bits,
        t_conv,
        t_exp_left,
        t_exp_right,
        t_gsw,
        do_expansion,
        db_dim_1,
        db_dim_2,
        instances,
        db_item_size,
        version,
    )
}

fn internal_params_for(
    nu_1: usize,
    nu_2: usize,
    p: u64,
    q2_bits: usize,
    t_exp_left: usize,
    moduli: &str,
) -> Params {
    ext_params_from_json(&format!(
        r#"
        {{
            "n": 1,
            "nu_1": {},
            "nu_2": {},
            "p": {},
            "q2_bits": {},
            "t_gsw": 3,
            "t_conv": 4,
            "t_exp_left": {},
            "t_exp_right": 2,
            "instances": 1,
            "db_item_size": 0,
            "moduli": {},
            "noise_width": 16.042421
        }}
        "#,
        nu_1, nu_2, p, q2_bits, t_exp_left, moduli
    ))
}

pub fn params_for_scenario(num_items: usize, item_size_bits: usize) -> Params {
    let total_db_bytes = num_items * item_size_bits / 8;
    let lwe_pt_word_bytes = 1;
    let num_items = total_db_bytes / lwe_pt_word_bytes;
    let num_tiles = num_items as f64 / (2048. * 2048.);
    let num_tiles_usize = num_tiles.ceil() as usize;
    let num_tiles_log2 = (num_tiles_usize as f64).log2().ceil() as usize;

    let (nu_1, nu_2) = if num_tiles_log2 % 2 == 0 {
        (num_tiles_log2 / 2, num_tiles_log2 / 2)
    } else {
        ((num_tiles_log2 + 1) / 2, (num_tiles_log2 - 1) / 2)
    };

    debug!("chose nu_1: {}, nu_2: {}", nu_1, nu_2);

    let p = 32768;
    let q2_bits = 28;
    let t_exp_left = 3;

    internal_params_for(nu_1, nu_2, p, q2_bits, t_exp_left, DEF_MOD_STR)
}

pub fn run_ypir<const K: usize>(
    num_items: usize,
    item_size_bits: usize,
    trials: usize,
) -> Measurement {
    let params = params_for_scenario(num_items, item_size_bits);
    run_ypir_on_params::<K>(params, trials)
}

pub fn run_ypir_batched(
    num_items: usize,
    item_size_bits: usize,
    num_clients: usize,
    trials: usize,
) -> Measurement {
    let params = params_for_scenario(num_items, item_size_bits);
    match num_clients {
        1 => run_ypir_on_params::<1>(params, trials),
        2 => run_ypir_on_params::<2>(params, trials),
        3 => run_ypir_on_params::<3>(params, trials),
        4 => run_ypir_on_params::<4>(params, trials),
        5 => run_ypir_on_params::<5>(params, trials),
        6 => run_ypir_on_params::<6>(params, trials),
        7 => run_ypir_on_params::<7>(params, trials),
        8 => run_ypir_on_params::<8>(params, trials),
        9 => run_ypir_on_params::<9>(params, trials),
        10 => run_ypir_on_params::<10>(params, trials),
        11 => run_ypir_on_params::<11>(params, trials),
        12 => run_ypir_on_params::<12>(params, trials),
        _ => panic!("Unsupported number of clients: {}", num_clients),
    }
}

pub trait Sample {
    fn sample() -> Self;
}

impl Sample for u8 {
    fn sample() -> Self {
        fastrand::u8(..)
    }
}

impl Sample for u16 {
    fn sample() -> Self {
        fastrand::u16(..)
    }
}

pub trait GetQPrime {
    /// The smaller reduced modulus, used on the second row of the encoding
    fn get_q_prime_1(&self) -> u64;

    /// The larger reduced modulus, used on the first row of the encoding
    fn get_q_prime_2(&self) -> u64;
}

impl GetQPrime for Params {
    fn get_q_prime_1(&self) -> u64 {
        1 << 20
    }

    fn get_q_prime_2(&self) -> u64 {
        if self.q2_bits == self.modulus_log2 {
            self.modulus
        } else {
            Q2_VALUES[self.q2_bits as usize]
        }
    }
}

impl GetQPrime for LWEParams {
    fn get_q_prime_1(&self) -> u64 {
        u64::MAX // unsupported
    }

    fn get_q_prime_2(&self) -> u64 {
        if self.q2_bits == (self.modulus as f64).log2().ceil() as usize {
            self.modulus
        } else {
            Q2_VALUES[self.q2_bits as usize]
        }
    }
}

pub fn run_ypir_on_params<const K: usize>(params: Params, trials: usize) -> Measurement {
    let lwe_params = LWEParams::default();

    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    let db_rows_padded = params.db_rows_padded();
    let db_cols = 1 << (params.db_dim_2 + params.poly_len_log2);

    let sqrt_n_bytes = db_cols * (lwe_params.pt_modulus as f64).log2().floor() as usize / 8;

    let mut rng = thread_rng();

    // RLWE reduced moduli
    let rlwe_q_prime_1 = params.get_q_prime_1();
    let rlwe_q_prime_2 = params.get_q_prime_2();

    // LWE reduced moduli
    let lwe_q_prime_bits = lwe_params.q2_bits as usize;

    // The number of bits represented by a plaintext RLWE coefficient
    let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
    // assert_eq!(pt_bits, 16);

    // The factor by which ciphertext values are bigger than plaintext values
    let blowup_factor = lwe_q_prime_bits as f64 / pt_bits as f64;
    debug!("blowup_factor: {}", blowup_factor);

    let mut smaller_params = params.clone();
    smaller_params.db_dim_1 = params.db_dim_2;
    smaller_params.db_dim_2 = ((blowup_factor * (lwe_params.n + 1) as f64) / params.poly_len as f64)
        .log2()
        .ceil() as usize;

    let out_rows = 1 << (smaller_params.db_dim_2 + params.poly_len_log2);
    let rho = 1 << smaller_params.db_dim_2; // rho

    debug!("rho: {}", rho);

    assert_eq!(smaller_params.db_dim_1, params.db_dim_2);
    assert!(out_rows as f64 >= (blowup_factor * (lwe_params.n + 1) as f64));

    // --

    let lwe_q_bits = (lwe_params.modulus as f64).log2().ceil() as usize;

    let rlwe_q_prime_1_bits = (rlwe_q_prime_1 as f64).log2().ceil() as usize;
    let rlwe_q_prime_2_bits = (rlwe_q_prime_2 as f64).log2().ceil() as usize;
    let simplepir_hint_bytes = (lwe_params.n * db_cols * lwe_q_prime_bits) / 8;
    let doublepir_hint_bytes = (params.poly_len * out_rows * rlwe_q_prime_2_bits) / 8;
    let simplepir_query_bytes = db_rows * lwe_q_bits / 8;
    let doublepir_query_bytes = db_cols * params.modulus_log2 as usize / 8;
    let simplepir_resp_bytes = (db_cols * lwe_q_prime_bits) / 8;
    let doublepir_resp_bytes = ((rho * params.poly_len) * rlwe_q_prime_2_bits
        + (rho * params.poly_len) * rlwe_q_prime_1_bits)
        / 8;
    debug!(
        "          \"simplepirHintBytes\": {},",
        simplepir_hint_bytes
    );
    debug!("          \"doublepirHintBytes\": {}", doublepir_hint_bytes);
    debug!(
        "          \"simplepirQueryBytes\": {},",
        simplepir_query_bytes
    );
    debug!(
        "          \"doublepirQueryBytes\": {},",
        doublepir_query_bytes
    );
    debug!(
        "          \"simplepirRespBytes\": {},",
        simplepir_resp_bytes
    );
    debug!(
        "          \"doublepirRespBytes\": {},",
        doublepir_resp_bytes
    );

    // --

    let now = Instant::now();
    let pt_iter = std::iter::repeat_with(|| u8::sample());
    let y_server = YServer::<u8>::new(&params, pt_iter, false, true);
    debug!("Created server in {} us", now.elapsed().as_micros());
    debug!(
        "Database of {} bytes",
        y_server.db().len() * std::mem::size_of::<u8>()
    );
    assert_eq!(
        y_server.db().len() * std::mem::size_of::<u8>(),
        db_rows_padded * db_cols * (lwe_params.pt_modulus as f64).log2().ceil() as usize / 8
    );

    // ================================================================
    // OFFLINE PHASE
    // ================================================================
    let mut measurements = vec![Measurement::default(); trials + 1];

    let start_offline_comp = Instant::now();
    let offline_values = y_server.perform_offline_precomputation(Some(&mut measurements[0]));
    let offline_server_time_ms = start_offline_comp.elapsed().as_millis();

    measurements[0].offline.server_time_ms = offline_server_time_ms as usize;
    measurements[0].offline.simplepir_hint_bytes = simplepir_hint_bytes;
    measurements[0].offline.doublepir_hint_bytes = doublepir_hint_bytes;
    measurements[0].online.simplepir_resp_bytes = simplepir_resp_bytes;
    measurements[0].online.doublepir_resp_bytes = doublepir_resp_bytes;

    let packed_query_row_sz = params.db_rows_padded();
    // let mut all_queries_packed = AlignedMemory64::new(K * packed_query_row_sz);

    for trial in 0..trials + 1 {
        debug!("trial: {}", trial);
        let mut measurement = &mut measurements[trial];
        measurement.offline.simplepir_hint_bytes = simplepir_hint_bytes;
        measurement.offline.doublepir_hint_bytes = doublepir_hint_bytes;
        measurement.online.simplepir_resp_bytes = simplepir_resp_bytes;
        measurement.online.doublepir_resp_bytes = doublepir_resp_bytes;

        // ================================================================
        // QUERY GENERATION PHASE
        // ================================================================
        let mut online_upload_bytes = 0;
        let mut queries = Vec::new();

        let mut clients = (0..K).map(|_| Client::init(&params)).collect::<Vec<_>>();

        for (_batch, client) in (0..K).zip(clients.iter_mut()) {
            let target_idx: usize = rng.gen::<usize>() % (db_rows * db_cols);
            let target_row = target_idx / db_cols;
            let target_col = target_idx % db_cols;
            debug!(
                "Target item: {} ({}, {})",
                target_idx, target_row, target_col
            );

            let start = Instant::now();
            client.generate_secret_keys();
            let sk_reg = &client.get_sk_reg();
            let pack_pub_params = raw_generate_expansion_params(
                &params,
                &sk_reg,
                params.poly_len_log2,
                params.t_exp_left,
                &mut ChaCha20Rng::from_entropy(),
                &mut ChaCha20Rng::from_entropy(),
            );
            let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params) / 2;
            debug!("pub params size: {} bytes", pub_params_size);

            let y_client = YClient::new(client, &params);
            let query_row = y_client.generate_query(SEED_0, params.db_dim_1, false, target_row);
            let query_row_last_row: &[u64] = &query_row[lwe_params.n * db_rows..];
            let mut aligned_query_packed = AlignedMemory64::new(query_row_last_row.len());
            aligned_query_packed
                .as_mut_slice()
                .copy_from_slice(&query_row_last_row);
            let packed_query_row = aligned_query_packed;
            let packed_query_row_u32 = packed_query_row
                .as_slice()
                .iter()
                .map(|x| *x as u32)
                .collect::<Vec<_>>();

            let query_col = y_client.generate_query(SEED_1, params.db_dim_2, true, target_col);
            let query_col_last_row = &query_col[params.poly_len * db_cols..];
            let packed_query_col = pack_query(&params, query_col_last_row);

            let query_size = query_row_last_row.len() * 4 + query_col_last_row.len() * 8;

            measurement.online.client_query_gen_time_ms = start.elapsed().as_millis() as usize;
            debug!("Generated query in {} us", start.elapsed().as_micros());

            online_upload_bytes = query_size + pub_params_size;
            debug!("Query size: {} bytes", online_upload_bytes);

            queries.push((
                y_client,
                target_idx,
                packed_query_row_u32,
                packed_query_col,
                pack_pub_params,
            ));
        }

        let mut all_queries_packed = vec![0u32; K * packed_query_row_sz];
        for (i, chunk_mut) in all_queries_packed
            .as_mut_slice()
            .chunks_mut(packed_query_row_sz)
            .enumerate()
        {
            (&mut chunk_mut[..db_rows]).copy_from_slice(queries[i].2.as_slice());
        }

        let mut offline_values = offline_values.clone();

        // ================================================================
        // ONLINE PHASE
        // ================================================================

        let start_online_comp = Instant::now();
        let response = y_server.perform_online_computation::<K>(
            &mut offline_values,
            &all_queries_packed,
            &queries
                .iter()
                .map(|x| (x.3.as_slice(), x.4.as_slice()))
                .collect::<Vec<_>>(),
            Some(&mut measurement),
        );
        let online_server_time_ms = start_online_comp.elapsed().as_millis();

        let online_download_bytes = get_size_bytes(&response);

        measurement.online.upload_bytes = online_upload_bytes;
        measurement.online.download_bytes = online_download_bytes;
        measurement.online.server_time_ms = online_server_time_ms as usize;
        measurement.online.sqrt_n_bytes = sqrt_n_bytes;
    }

    // discard the first measurement (if there were multiple trials)
    if trials > 1 {
        measurements.remove(0);
    }

    let mut final_measurement = measurements[0].clone();
    final_measurement.online.server_time_ms = mean(
        &measurements
            .iter()
            .map(|m| m.online.server_time_ms)
            .collect::<Vec<_>>(),
    )
    .round() as usize;
    final_measurement.online.all_server_times_ms = measurements
        .iter()
        .map(|m| m.online.server_time_ms)
        .collect::<Vec<_>>();
    final_measurement.online.std_dev_server_time_ms =
        std_dev(&final_measurement.online.all_server_times_ms);

    final_measurement
}

fn mean(xs: &[usize]) -> f64 {
    xs.iter().map(|x| *x as f64).sum::<f64>() / xs.len() as f64
}

fn std_dev(xs: &[usize]) -> f64 {
    let mean = mean(xs);
    let mut variance = 0.;
    for x in xs {
        variance += (*x as f64 - mean).powi(2);
    }
    (variance / xs.len() as f64).sqrt()
}

#[cfg(test)]
mod test {
    use super::*;
    use test_log::test;

    #[test]
    fn test_ypir() {
        run_ypir_batched(1 << 30, 1, 1, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_1gb() {
        run_ypir_batched(1 << 33, 1, 1, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_2gb() {
        run_ypir_batched(1 << 34, 1, 1, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_4gb() {
        run_ypir_batched(1 << 35, 1, 9, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_8gb() {
        run_ypir_batched(1 << 36, 1, 8, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_16gb() {
        run_ypir_batched(1 << 37, 1, 1, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_32gb() {
        run_ypir_batched(1 << 38, 1, 1, 5);
    }

    #[test]
    #[ignore]
    fn test_batched_4_ypir() {
        run_ypir_batched(1 << 30, 1, 8, 5);
    }

    fn pp_poly(params: &Params, poly: &[u64], scale: u64) {
        let mut vals1: Vec<i64> = (&poly).iter().map(|x| *x as i64).collect();
        // let mut vals2: Vec<i64> = (&poly[poly.len() / 2..poly.len() / 2 + params.poly_len / 2])
        //     .iter()
        //     .map(|x| *x as i64)
        //     .collect();

        let m = params.modulus as i64;
        for i in 0..vals1.len() {
            vals1[i] %= m;
            if vals1[i] >= (m / 2) {
                vals1[i] -= m;
            }
            vals1[i] /= scale as i64;
            // if vals2[i] >= m / 2 {
            //     vals2[i] -= m;
            // }
        }

        debug!("m: {}", m);
        debug!("{:?}", vals1);
        // debug!("{:?}... |\n{:?}...", vals1, vals2);
    }

    #[test]
    fn test_automorph() {
        let params = params_for_scenario(1 << 30, 1);
        let mut c0 = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            c0.data[i] = i as u64;
        }

        let mut c1 = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            c1.data[i] = 100000 + i as u64;
        }

        debug!("init:");
        pp_poly(&params, c0.as_slice(), 1);
        pp_poly(&params, c1.as_slice(), 1);
        debug!("");

        let ell = 1;
        let t = (1 << ell) + 1;
        let mut x_to_the = PolyMatrixRaw::zero(&params, 1, 1);
        x_to_the.data[params.poly_len / 2] = 1;
        let term_0 = &c0 + &scalar_multiply_alloc(&x_to_the.ntt(), &c1.ntt()).raw();
        let auot_inp = &c0 + &-(&scalar_multiply_alloc(&x_to_the.ntt(), &c1.ntt()).raw());
        let term_1 = automorph_alloc(&auot_inp, t);
        let res = &term_0 + &term_1;

        // pp_poly(&params, res.as_slice(), 1);
        pp_poly(&params, &res.as_slice()[..10], 2);
        pp_poly(&params, &res.as_slice()[1024..1024 + 10], 2);

        // for ell in 0..params.poly_len_log2 {
        //     let t = (1 << ell) + 1;
        //     let poly_auto = automorph_alloc(&poly, t);
        //     debug!("poly_auto (ell: {}, t: {}):", ell, t,);
        //     pp_poly(&params, poly_auto.as_slice());
        //     debug!("");
        // }
    }
}
