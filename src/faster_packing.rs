//! Implements a faster version of packing LWE->RLWE in the preprocessing model.

use std::time::Instant;

use log::debug;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use spiral_rs::{
    gadget::{build_gadget, gadget_invert_alloc},
    params::*,
    poly::*,
};

use crate::client::get_fresh_reg_public_key;

/// Preprocess the input (query-independent) ciphertexts into a hint.
pub fn prep_pack<'a>(
    params: &'a Params,
    v_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
) -> Vec<PolyMatrixNTT<'a>> {
    todo!("Implement faster packing")
}

fn multiply_automorphism_keys<'a>(
    outer_ct: &PolyMatrixNTT<'a>,
    t_outer: usize,
    inner_ct: &PolyMatrixNTT<'a>,
    t_exp: usize,
) -> PolyMatrixNTT<'a> {
    let inner_ct_0 = inner_ct.submatrix(0, 0, 1, inner_ct.cols).raw();
    let inner_ct_0 = automorph_alloc(&inner_ct_0, t_outer);
    let inner_ct_1 = inner_ct.submatrix(1, 0, 1, inner_ct.cols).raw();
    let inner_ct_1 = automorph_alloc(&inner_ct_1, t_outer);

    let w_prime = outer_ct * &gadget_invert_alloc(t_exp, &inner_ct_0).ntt();
    let w_prime_next = &w_prime + &inner_ct_1.pad_top(1).ntt();

    w_prime_next
}

/// Process the given automorphism keys into the full set.
pub fn online_generate_full_automorphism_keys<'a>(
    params: &'a Params,
    pub_params: &[PolyMatrixNTT<'a>],
    auto_amounts: &[usize],
) -> Vec<PolyMatrixNTT<'a>> {
    // For each automorphism key in pub_params, either apply it or don't in a
    // binary tree, where the root is a ct with [g^T \\ 0]. The index of an
    // output (index in range 0..2^(pub_params.len())) has an MSB corresponding
    // to whether the first automorphism key is applied, the next MSB to the
    // second automorphism key, and so on.

    let t_exp = params.t_exp_left;
    let mut init_raw = PolyMatrixRaw::zero(params, 2, t_exp);
    init_raw.copy_into(&build_gadget(params, 1, t_exp), 0, 0);
    let init = init_raw.ntt();
    let mut res = vec![init];
    for i in 0..pub_params.len() {
        let mut new_res = Vec::new();
        for j in 0..res.len() {
            let inp_ct = &res[j];
            let new_ct =
                multiply_automorphism_keys(&pub_params[i], auto_amounts[i], &inp_ct, t_exp);
            new_res.push(inp_ct.clone());
            new_res.push(new_ct);
        }
        res = new_res;
    }

    res
}

/// Online pack the input (query-dependent) ciphertexts, using the hint.
pub fn online_pack<'a>(
    params: &'a Params,
    v_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
    hint: &[PolyMatrixNTT<'a>],
) -> PolyMatrixNTT<'a> {
    if v_cts.len() != 1 {
        todo!("Implement online packing for multiple ciphertexts");
    }
    if hint.len() > 0 {
        todo!("Implement online packing with hint");
    }

    let t_exp = params.t_exp_left;
    let rounds = pub_params.len();
    let auto_amounts = automorphism_amounts(params, rounds);
    println!("auto_amounts: {:?}", auto_amounts);

    let now = Instant::now();
    let full_automorph_keys =
        online_generate_full_automorphism_keys(params, pub_params, &auto_amounts);
    println!(
        "full_automorph_keys took: {:?} us",
        now.elapsed().as_micros()
    );
    let ct = v_cts[0].clone();
    let mut sum = ct.clone();
    for i in 1..(1 << rounds) {
        let mut t = 1;
        for idx_bit in 0..rounds {
            let bit = (i >> idx_bit) & 1;
            if bit == 1 {
                t *= auto_amounts[auto_amounts.len() - 1 - idx_bit];
                t = t % (2 * params.poly_len);
            }
        }
        // println!("i: {}, t: {}", i, t);

        let w = &full_automorph_keys[i];
        let ct_auto = apply_automorphsim_key(&ct, w, t, t_exp);
        sum = &sum + &ct_auto;
    }

    sum
}

pub fn generate_automorphism_keys<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    t_vals: &[usize],
    t_exp: usize,
) -> Vec<PolyMatrixNTT<'a>> {
    let mut rng = ChaCha20Rng::from_entropy();
    let mut rng_pub = ChaCha20Rng::from_entropy();

    let g_exp = build_gadget(params, 1, t_exp);
    debug!("using gadget base {}", g_exp.get_poly(0, 1)[0]);
    let g_exp_ntt = g_exp.ntt();
    let mut res = Vec::new();

    for i in 0..t_vals.len() {
        let t = t_vals[i]; //(params.poly_len / (1 << i)) + 1;
        let tau_sk_reg = automorph_alloc(&sk_reg, t);
        let prod = &tau_sk_reg.ntt() * &g_exp_ntt;

        // let w_exp_i = client.encrypt_matrix_reg(&prod, rng, rng_pub);
        let sample = get_fresh_reg_public_key(params, &sk_reg, t_exp, &mut rng, &mut rng_pub);
        let w_exp_i = &sample + &prod.pad_top(1);
        res.push(w_exp_i);
    }

    res
}

pub fn apply_automorphsim_key<'a>(
    ct: &PolyMatrixNTT<'a>,
    w: &PolyMatrixNTT<'a>,
    t: usize,
    t_exp: usize,
) -> PolyMatrixNTT<'a> {
    let tau_ct = automorph_alloc(&ct.submatrix(0, 0, 1, ct.cols).raw(), t);
    let prod = w * &gadget_invert_alloc(t_exp, &tau_ct).ntt();
    let res = &prod
        + &automorph_alloc(&ct.submatrix(1, 0, 1, ct.cols).raw(), t)
            .ntt()
            .pad_top(1);
    res
}

pub fn automorphism_amounts(params: &Params, rounds: usize) -> Vec<usize> {
    let mut auto_amounts = Vec::new();
    for i in 0..rounds {
        auto_amounts.push(params.poly_len / (1 << i) + 1);
    }
    auto_amounts
}

#[cfg(test)]
mod test {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;
    use spiral_rs::{
        arith::{multiply_uint_mod, rescale},
        client::*,
        discrete_gaussian::DiscreteGaussian,
        number_theory::invert_uint_mod,
        poly::*,
        util::*,
    };
    use test_log::test;

    use super::*;

    fn get_params() -> Params {
        let mut params = get_test_params();
        params.t_exp_left = 8;
        params
    }

    fn generate_lwe_ciphertext<'a>(
        params: &'a Params,
        sk_reg: &PolyMatrixRaw<'a>,
        val: u64,
    ) -> PolyMatrixNTT<'a> {
        let dg = DiscreteGaussian::init(params.noise_width);
        let mut rng = ChaCha20Rng::from_entropy();
        let mut ct = PolyMatrixRaw::zero(params, 2, 1);
        let mut dot_prod = 0;
        for i in 0..params.poly_len {
            let a_val = rng.gen_range(0..params.modulus);
            let sk_val = sk_reg.data[i];
            dot_prod += ((a_val as u128 * sk_val as u128) % (params.modulus as u128)) as u64;
            dot_prod %= params.modulus;
            ct.data[i] = a_val;
        }
        let e = dg.sample(params.modulus, &mut rng);

        let b_val = (val + e + params.modulus - dot_prod) % params.modulus;
        ct.data[params.poly_len] = b_val;

        ct.ntt()
    }

    fn decrypt_lwe_ciphertext<'a>(
        params: &'a Params,
        sk_reg: &PolyMatrixRaw<'a>,
        ct_ntt: &PolyMatrixNTT<'a>,
    ) -> u64 {
        let ct = ct_ntt.raw();
        let mut dot_prod = 0;
        for i in 0..params.poly_len {
            let a_val = ct.data[i];
            let sk_val = sk_reg.data[i];
            dot_prod += ((a_val as u128 * sk_val as u128) % (params.modulus as u128)) as u64;
            dot_prod %= params.modulus;
        }

        let b_val = ct.data[params.poly_len];
        let decrypted_val = (b_val + dot_prod) % params.modulus;
        decrypted_val
    }

    fn decode_result(params: &Params, result: u64) -> u64 {
        // rescale
        let rescaled_val = rescale(result, params.modulus, params.pt_modulus);
        rescaled_val
    }

    #[test]
    fn test_lwe() {
        let params = get_params();
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let sk_reg = client.get_sk_reg();
        let val = 17;
        let scale_k = params.modulus / params.pt_modulus;
        let scaled_val = val * scale_k;
        let ct = generate_lwe_ciphertext(&params, &sk_reg, scaled_val);
        let result = decrypt_lwe_ciphertext(&params, &sk_reg, &ct);
        let decoded_result = decode_result(&params, result);
        assert_eq!(decoded_result, val);
    }

    fn generate_range_ct<'a>(
        params: &'a Params,
        client: &Client<'a>,
        first_val: u64,
    ) -> PolyMatrixNTT<'a> {
        let scale_k = params.modulus / params.pt_modulus;
        let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            let val = if i == 0 {
                first_val
            } else {
                i as u64 % params.pt_modulus
            };
            pt.data[i] = multiply_uint_mod(val, scale_k, params.modulus);
        }
        let ct = client.encrypt_matrix_reg(
            &pt.ntt(),
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_entropy(),
        );

        ct
    }

    fn get_result<'a>(
        params: &'a Params,
        client: &Client<'a>,
        ct: &PolyMatrixNTT<'a>,
    ) -> PolyMatrixRaw<'a> {
        let result = client.decrypt_matrix_reg(&ct);
        let mut result_raw = result.raw();
        for i in 0..params.poly_len {
            result_raw.data[i] = rescale(result_raw.data[i], params.modulus, params.pt_modulus);
        }
        result_raw
    }

    #[test]
    fn test_apply_automorphsim_key() {
        let params = get_params();
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let t_amount = 3073;
        let t_exp = params.t_exp_left;

        let ct = generate_range_ct(&params, &client, 0);
        let w = generate_automorphism_keys(
            &params,
            &client.get_sk_reg(),
            &[t_amount],
            params.t_exp_left,
        );
        let ct_automorph = apply_automorphsim_key(&ct, &w[0], t_amount, t_exp);
        let ct_final = &ct_automorph + &ct;
        let result = get_result(&params, &client, &ct_final);
        println!("result: {:?}", &result.as_slice());
        // assert_eq!(
        //     result.as_slice()[0..10],
        //     [0, 255, 2, 253, 4, 251, 6, 249, 8, 247] // 0 -1 2 -3 4 -5 6 -7 8 -9
        // );
    }

    #[test]
    fn test_multiply_automorphism_keys() {
        let params = get_params();
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let auto_amounts = [2049, 1025, 3073];
        assert_eq!((2049 * 1025) % 4096, 3073);

        let t_exp = params.t_exp_left;

        let ct = generate_range_ct(&params, &client, 0);
        let w = generate_automorphism_keys(
            &params,
            &client.get_sk_reg(),
            &auto_amounts,
            params.t_exp_left,
        );

        let ct_0_automorph = apply_automorphsim_key(&ct, &w[0], auto_amounts[0], t_exp);
        let ct_1_automorph = apply_automorphsim_key(&ct, &w[1], auto_amounts[1], t_exp);
        let ct_2_automorph = apply_automorphsim_key(&ct, &w[2], auto_amounts[2], t_exp);

        let result_0 = get_result(&params, &client, &ct_0_automorph);
        println!("result_0      : {:?}", &result_0.as_slice()[0..10]);
        let result_1 = get_result(&params, &client, &ct_1_automorph);
        println!("result_1      : {:?}", &result_1.as_slice()[0..10]);
        let result_2 = get_result(&params, &client, &ct_2_automorph);
        println!("result_2      : {:?}", &result_2.as_slice()[0..10]);

        let w_2_guess = multiply_automorphism_keys(&w[0], auto_amounts[0], &w[1], t_exp);

        let ct_2_guess = apply_automorphsim_key(&ct, &w_2_guess, auto_amounts[2], t_exp);
        let result_2_guess = get_result(&params, &client, &ct_2_guess);
        println!("result_2_guess: {:?}", &result_2_guess.as_slice()[0..10]);
        assert_eq!(result_2.as_slice(), result_2_guess.as_slice());
    }

    #[test]
    fn test_online_generate_full_automorphism_keys() {
        let params = get_params();
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let auto_amounts = [2049, 1025, 3073];
        assert_eq!((2049 * 1025) % 4096, 3073);

        let t_exp = params.t_exp_left;

        let ct = generate_range_ct(&params, &client, 0);
        let w = generate_automorphism_keys(
            &params,
            &client.get_sk_reg(),
            &auto_amounts,
            params.t_exp_left,
        );

        let ct_0_automorph = apply_automorphsim_key(&ct, &w[0], auto_amounts[0], t_exp);
        let ct_1_automorph = apply_automorphsim_key(&ct, &w[1], auto_amounts[1], t_exp);
        let ct_2_automorph = apply_automorphsim_key(&ct, &w[2], auto_amounts[2], t_exp);

        let result_0 = get_result(&params, &client, &ct_0_automorph);
        println!("result_0      : {:?}", &result_0.as_slice()[0..10]);
        let result_1 = get_result(&params, &client, &ct_1_automorph);
        println!("result_1      : {:?}", &result_1.as_slice()[0..10]);
        let result_2 = get_result(&params, &client, &ct_2_automorph);
        println!("result_2      : {:?}", &result_2.as_slice()[0..10]);

        let full_automorph_keys =
            online_generate_full_automorphism_keys(&params, &w[..2], &auto_amounts[..2]);
        let w_2_guess = &full_automorph_keys[3];

        let ct_2_guess = apply_automorphsim_key(&ct, &w_2_guess, auto_amounts[2], t_exp);
        let result_2_guess = get_result(&params, &client, &ct_2_guess);
        println!("result_2_guess: {:?}", &result_2_guess.as_slice()[0..10]);
        assert_eq!(result_2.as_slice(), result_2_guess.as_slice());
    }

    #[test]
    fn test_online_pack() {
        let params = get_params();
        let mut client = Client::init(&params);
        client.generate_secret_keys();

        let rounds = params.poly_len_log2;
        let auto_amounts = automorphism_amounts(&params, rounds);

        let val = 7;
        let scale_factor = invert_uint_mod(1 << rounds, params.modulus).unwrap();
        let val_to_pack = multiply_uint_mod(val, scale_factor, params.modulus);
        let mut ct = generate_range_ct(&params, &client, val_to_pack);
        let mut ct_raw = ct.raw();
        // destroy all but the first value
        for i in 1..params.poly_len {
            ct_raw.get_poly_mut(1, 0)[i] = 0;
        }
        ct = ct_raw.ntt();
        let w = generate_automorphism_keys(
            &params,
            &client.get_sk_reg(),
            &auto_amounts,
            params.t_exp_left,
        );

        let ct_online = online_pack(&params, &[ct], &w, &[]);
        let result_online = get_result(&params, &client, &ct_online);
        println!("result_online : {:?}", &result_online.as_slice());
    }
}
