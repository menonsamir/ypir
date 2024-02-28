use std::time::Instant;

use log::debug;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::arith::rescale;
use spiral_rs::poly::{PolyMatrix, PolyMatrixRaw};
use spiral_rs::{client::*, params::*};

use crate::bits::{read_bits, u64s_to_contiguous_bytes};
use crate::modulus_switch::ModulusSwitch;
use crate::noise_analysis::YPIRSchemeParams;
use crate::packing::condense_matrix;

use super::{client::*, lwe::LWEParams, measurement::*, params::*, server::*};

pub const STATIC_PUBLIC_SEED: [u8; 32] = [0u8; 32];
pub const SEED_0: u8 = 0;
pub const SEED_1: u8 = 1;

pub const STATIC_SEED_2: [u8; 32] = [
    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

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
    let measurement = match num_clients {
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
    };
    debug!("{:#?}", measurement);
    measurement
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

    let packed_query_row_sz = params.db_rows_padded();
    // let mut all_queries_packed = AlignedMemory64::new(K * packed_query_row_sz);

    for trial in 0..trials + 1 {
        debug!("trial: {}", trial);
        let mut measurement = &mut measurements[trial];
        measurement.offline.server_time_ms = offline_server_time_ms as usize;
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
                &mut ChaCha20Rng::from_seed(STATIC_SEED_2),
            );
            // let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params) / 2;
            let mut pack_pub_params_row_1s = pack_pub_params.to_vec();
            for i in 0..pack_pub_params.len() {
                pack_pub_params_row_1s[i] =
                    pack_pub_params[i].submatrix(1, 0, 1, pack_pub_params[i].cols);
                pack_pub_params_row_1s[i] = condense_matrix(&params, &pack_pub_params_row_1s[i]);
            }
            let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params_row_1s);
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
                pack_pub_params_row_1s,
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
        let responses = y_server.perform_online_computation::<K>(
            &mut offline_values,
            &all_queries_packed,
            &queries
                .iter()
                .map(|x| (x.3.as_slice(), x.4.as_slice()))
                .collect::<Vec<_>>(),
            Some(&mut measurement),
        );
        let online_server_time_ms = start_online_comp.elapsed().as_millis();
        let online_download_bytes = get_size_bytes(&responses); // TODO: this is not quite right for multiple clients

        // check correctness
        for (response_switched, (y_client, target_idx, _, _, _)) in
            responses.iter().zip(queries.iter())
        {
            let (target_row, target_col) = (target_idx / db_cols, target_idx % db_cols);
            let corr_result = y_server.get_elem(target_row, target_col).to_u64();

            let scheme_params = YPIRSchemeParams::from_params(&params, &lwe_params);
            let log2_corr_err = scheme_params.delta().log2();
            let log2_expected_outer_noise = scheme_params.expected_outer_noise().log2();
            debug!("log2_correctness_err: {}", log2_corr_err);
            debug!("log2_expected_outer_noise: {}", log2_expected_outer_noise);

            let start_decode = Instant::now();

            debug!("rescaling response...");
            let mut response = Vec::new();
            for ct_bytes in response_switched.iter() {
                let ct = PolyMatrixRaw::recover(&params, rlwe_q_prime_1, rlwe_q_prime_2, ct_bytes);
                response.push(ct);
            }

            debug!("decrypting outer cts...");
            let outer_ct = response
                .iter()
                .flat_map(|ct| {
                    decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len)
                        .as_slice()
                        .to_vec()
                })
                .collect::<Vec<_>>();
            assert_eq!(outer_ct.len(), out_rows);
            // debug!("outer_ct: {:?}", &outer_ct[..]);
            let outer_ct_t_u8 = u64s_to_contiguous_bytes(&outer_ct, pt_bits);

            let mut inner_ct = PolyMatrixRaw::zero(&params, 2, 1);
            let mut bit_offs = 0;
            let lwe_q_prime = lwe_params.get_q_prime_2();
            let special_offs =
                ((lwe_params.n * lwe_q_prime_bits) as f64 / pt_bits as f64).ceil() as usize;
            for z in 0..lwe_params.n {
                let val = read_bits(&outer_ct_t_u8, bit_offs, lwe_q_prime_bits);
                bit_offs += lwe_q_prime_bits;
                assert!(
                    val < lwe_q_prime,
                    "val: {}, lwe_q_prime: {}",
                    val,
                    lwe_q_prime
                );
                inner_ct.data[z] = rescale(val, lwe_q_prime, lwe_params.modulus);
            }

            let mut val = 0;
            for i in 0..blowup_factor.ceil() as usize {
                val |= outer_ct[special_offs + i] << (i * pt_bits);
            }
            assert!(
                val < lwe_q_prime,
                "val: {}, lwe_q_prime: {}",
                val,
                lwe_q_prime
            );
            debug!("got b_val of: {}", val);
            inner_ct.data[lwe_params.n] = rescale(val, lwe_q_prime, lwe_params.modulus);

            debug!("decrypting inner ct...");
            // let plaintext = decrypt_ct_reg_measured(y_client.client(), &params, &inner_ct.ntt(), 1);
            // let final_result = plaintext.data[0];
            let inner_ct_as_u32 = inner_ct
                .as_slice()
                .iter()
                .take(lwe_params.n + 1)
                .map(|x| *x as u32)
                .collect::<Vec<_>>();
            let decrypted = y_client.lwe_client().decrypt(&inner_ct_as_u32);
            let final_result = rescale(decrypted as u64, lwe_params.modulus, lwe_params.pt_modulus);

            measurement.online.client_decode_time_ms = start_decode.elapsed().as_millis() as usize;

            debug!("got {}, expected {}", final_result, corr_result);
            // debug!("was correct? {}", final_result == corr_result);
            assert_eq!(final_result, corr_result);
        }

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
    fn test_ypir_basic() {
        run_ypir_batched(1 << 30, 1, 1, 1);
    }

    #[test]
    fn test_ypir_many_clients() {
        run_ypir_batched(1 << 30, 1, 2, 1);
    }

    #[test]
    fn test_ypir_many_clients_and_trials() {
        run_ypir_batched(1 << 30, 1, 2, 5);
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
        run_ypir_batched(1 << 35, 1, 1, 5);
    }

    #[test]
    #[ignore]
    fn test_ypir_8gb() {
        run_ypir_batched(1 << 36, 1, 1, 5);
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
        run_ypir_batched(1 << 30, 1, 4, 5);
    }
}
