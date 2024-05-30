use log::debug;
use rand::rngs::OsRng;
use rand::{CryptoRng, Rng, RngCore, SeedableRng};
use rand_chacha::ChaCha20Rng;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{
    arith::*, client::*, discrete_gaussian::*, gadget::*, number_theory::*, params::*, poly::*,
};

use crate::measurement::get_vec_pm_size_bytes;
use crate::packing::condense_matrix;
use crate::serialize::pack_vec_pm;

use super::convolution::negacyclic_matrix_u32;
use super::{constants::*, lwe::*, noise_analysis::measure_noise_width_squared, util::*};

pub fn rlwe_to_lwe<'a>(params: &'a Params, ct: &PolyMatrixRaw<'a>) -> Vec<u64> {
    let a = ct.get_poly(0, 0);
    let mut negacylic_a = negacyclic_matrix(&a, params.modulus);
    negacylic_a.extend(ct.get_poly(1, 0));

    negacylic_a
}

pub fn pack_query(params: &Params, query: &[u64]) -> AlignedMemory64 {
    let query_packed = query
        .iter()
        .enumerate()
        .map(|(_i, x)| {
            let crt0 = (*x) % params.moduli[0];
            let crt1 = (*x) % params.moduli[1];
            crt0 | (crt1 << 32)
        })
        .collect::<Vec<_>>();
    let mut aligned_query_packed = AlignedMemory64::new(query_packed.len());
    aligned_query_packed
        .as_mut_slice()
        .copy_from_slice(&query_packed);
    aligned_query_packed
}

pub fn get_reg_sample<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    rng: &mut ChaCha20Rng,
    rng_pub: &mut ChaCha20Rng,
) -> PolyMatrixNTT<'a> {
    let a = PolyMatrixRaw::random_rng(params, 1, 1, rng_pub);
    let e = PolyMatrixRaw::noise(
        params,
        1,
        1,
        &DiscreteGaussian::init(params.noise_width),
        rng,
    );
    let b_p = &sk_reg.ntt() * &a.ntt();
    let b = &e.ntt() + &b_p;
    let mut p = PolyMatrixNTT::zero(params, 2, 1);
    p.copy_into(&(-&a).ntt(), 0, 0);
    p.copy_into(&b, 1, 0);
    p
}

pub fn get_fresh_reg_public_key<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    m: usize,
    rng: &mut ChaCha20Rng,
    rng_pub: &mut ChaCha20Rng,
) -> PolyMatrixNTT<'a> {
    let mut p = PolyMatrixNTT::zero(params, 2, m);

    for i in 0..m {
        p.copy_into(&get_reg_sample(params, sk_reg, rng, rng_pub), 0, i);
    }
    p
}

pub fn raw_generate_expansion_params<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    num_exp: usize,
    m_exp: usize,
    rng: &mut ChaCha20Rng,
    rng_pub: &mut ChaCha20Rng,
) -> Vec<PolyMatrixNTT<'a>> {
    let g_exp = build_gadget(params, 1, m_exp);
    debug!("using gadget base {}", g_exp.get_poly(0, 1)[0]);
    let g_exp_ntt = g_exp.ntt();
    let mut res = Vec::new();

    for i in 0..num_exp {
        let t = (params.poly_len / (1 << i)) + 1;
        let tau_sk_reg = automorph_alloc(&sk_reg, t);
        let prod = &tau_sk_reg.ntt() * &g_exp_ntt;

        // let w_exp_i = client.encrypt_matrix_reg(&prod, rng, rng_pub);
        let sample = get_fresh_reg_public_key(params, &sk_reg, m_exp, rng, rng_pub);
        let w_exp_i = &sample + &prod.pad_top(1);
        res.push(w_exp_i);
    }

    res
}

pub fn decrypt_ct_reg_measured<'a>(
    client: &Client<'a>,
    params: &'a Params,
    ct: &PolyMatrixNTT<'a>,
    coeffs_to_measure: usize,
) -> PolyMatrixRaw<'a> {
    let dec_result = client.decrypt_matrix_reg(ct).raw();

    let mut dec_rescaled = PolyMatrixRaw::zero(&params, dec_result.rows, dec_result.cols);
    for z in 0..dec_rescaled.data.len() {
        dec_rescaled.data[z] = rescale(dec_result.data[z], params.modulus, params.pt_modulus);
    }

    // measure noise width
    let s_2 = measure_noise_width_squared(params, client, ct, &dec_rescaled, coeffs_to_measure);
    debug!("log2(measured noise): {}", s_2.log2());

    dec_rescaled
}

pub struct YClient<'a> {
    inner: &'a mut Client<'a>,
    params: &'a Params,
    lwe_client: LWEClient,
}

pub fn get_seed(public_seed_idx: u8) -> [u8; 32] {
    let mut seed = STATIC_PUBLIC_SEED;
    seed[0] = public_seed_idx;
    seed
}

pub fn generate_matrix_ring(
    rng_pub: &mut ChaCha20Rng,
    n: usize,
    rows: usize,
    cols: usize,
) -> Vec<u32> {
    assert_eq!(rows % n, 0);
    assert_eq!(cols % n, 0);
    let rows_outer = rows / n;
    let cols_outer = cols / n;

    let mut out = vec![0u32; rows * cols];
    for i in 0..rows_outer {
        for j in 0..cols_outer {
            let mut a = vec![0u32; n];
            for idx in 0..n {
                a[idx] = rng_pub.sample::<u32, _>(rand::distributions::Standard);
            }

            let mat = negacyclic_matrix_u32(&a);
            for k in 0..n {
                for l in 0..n {
                    let idx = (i * n + k) * cols + (j * n + l);
                    out[idx] = mat[k * n + l];
                }
            }
        }
    }

    out
}

impl<'a> YClient<'a> {
    pub fn new(inner: &'a mut Client<'a>, params: &'a Params) -> Self {
        Self {
            inner,
            params,
            lwe_client: LWEClient::new(LWEParams::default()),
        }
    }

    pub fn lwe_client(&self) -> &LWEClient {
        &self.lwe_client
    }

    fn rlwes_to_lwes(&self, ct: &[PolyMatrixRaw<'a>]) -> Vec<u64> {
        let v = ct
            .iter()
            .map(|ct| rlwe_to_lwe(self.params, ct))
            .collect::<Vec<_>>();
        concat_horizontal(&v, self.params.poly_len + 1, self.params.poly_len)
    }

    pub fn generate_query_impl(
        &self,
        public_seed_idx: u8,
        dim_log2: usize,
        packing: bool,
        index: usize,
    ) -> Vec<PolyMatrixRaw<'a>> {
        // let db_cols = 1 << (self.params.db_dim_2 + self.params.poly_len_log2);
        // let idx_dim1 = index / db_cols;

        let multiply_ct = true;

        let mut rng_pub = ChaCha20Rng::from_seed(get_seed(public_seed_idx));

        // Generate dim1_bits LWE samples under public randomness
        let mut out = Vec::new();

        let scale_k = self.params.modulus / self.params.pt_modulus;

        for i in 0..(1 << dim_log2) {
            let mut scalar = PolyMatrixRaw::zero(self.params, 1, 1);
            let is_nonzero = i == (index / self.params.poly_len);

            if is_nonzero {
                scalar.data[index % self.params.poly_len] = scale_k;
            }

            if packing {
                let factor =
                    invert_uint_mod(self.params.poly_len as u64, self.params.modulus).unwrap();
                scalar = scalar_multiply_alloc(
                    &PolyMatrixRaw::single_value(self.params, factor).ntt(),
                    &scalar.ntt(),
                )
                .raw();
            }

            // if public_seed_idx == SEED_0 {
            //     out.push(scalar.pad_top(1));
            //     continue;
            // }

            let ct = if multiply_ct {
                let factor =
                    invert_uint_mod(self.params.poly_len as u64, self.params.modulus).unwrap();

                self.inner.encrypt_matrix_scaled_reg(
                    &scalar.ntt(),
                    &mut ChaCha20Rng::from_entropy(),
                    &mut rng_pub,
                    factor,
                )
            } else {
                self.inner.encrypt_matrix_reg(
                    &scalar.ntt(),
                    &mut ChaCha20Rng::from_entropy(),
                    &mut rng_pub,
                )
            };

            // let mut ct = self.inner.encrypt_matrix_reg(
            //     &scalar.ntt(),
            //     &mut ChaCha20Rng::from_entropy(),
            //     &mut rng_pub,
            // );

            // if multiply_ct && packing {
            //     let factor =
            //         invert_uint_mod(self.params.poly_len as u64, self.params.modulus).unwrap();
            //     ct = scalar_multiply_alloc(
            //         &PolyMatrixRaw::single_value(self.params, factor).ntt(),
            //         &ct,
            //     );
            // };

            // if multiply_error && is_nonzero && packing {
            //     let factor =
            //         invert_uint_mod(self.params.poly_len as u64, self.params.modulus).unwrap();
            //     ct = scalar_multiply_alloc(
            //         &PolyMatrixRaw::single_value(self.params, factor).ntt(),
            //         &ct,
            //     );
            // }

            let ct_raw = ct.raw();
            // let ct_0_nega = negacyclic_perm(ct_raw.get_poly(0, 0), 0, self.params.modulus);
            // let ct_1_nega = negacyclic_perm(ct_raw.get_poly(1, 0), 0, self.params.modulus);
            // let mut ct_nega = PolyMatrixRaw::zero(self.params, 2, 1);
            // ct_nega.get_poly_mut(0, 0).copy_from_slice(&ct_0_nega);
            // ct_nega.get_poly_mut(1, 0).copy_from_slice(&ct_1_nega);

            // self-test
            // {
            //     let test_ct = self.inner.encrypt_matrix_reg(
            //         &PolyMatrixRaw::single_value(self.params, scale_k * 7).ntt(),
            //         &mut ChaCha20Rng::from_entropy(),
            //         &mut ChaCha20Rng::from_entropy(),
            //     );
            //     let lwe = rlwe_to_lwe(self.params, &test_ct.raw());
            //     let result = self.decode_response(&lwe);
            //     assert_eq!(result[0], 7);
            // }

            out.push(ct_raw);
        }

        out
    }

    pub fn generate_query(
        &self,
        public_seed_idx: u8,
        dim_log2: usize,
        packing: bool,
        index_row: usize,
    ) -> Vec<u64> {
        if public_seed_idx == SEED_0 && !packing {
            let lwe_params = LWEParams::default();
            let dim = 1 << (dim_log2 + self.params.poly_len_log2);

            // lwes must be (n + 1) x (dim) matrix
            let mut lwes = vec![0u64; (lwe_params.n + 1) * dim];

            let scale_k = lwe_params.scale_k() as u32;
            let mut vals_to_encrypt = vec![0u32; dim];
            vals_to_encrypt[index_row] = scale_k;

            let mut rng_pub = ChaCha20Rng::from_seed(get_seed(public_seed_idx));

            for i in (0..dim).step_by(lwe_params.n) {
                let out = self
                    .lwe_client
                    .encrypt_many(&mut rng_pub, &vals_to_encrypt[i..i + lwe_params.n])
                    .iter()
                    .map(|x| *x as u64)
                    .collect::<Vec<_>>();
                assert_eq!(out.len(), (lwe_params.n + 1) * lwe_params.n);
                for r in 0..lwe_params.n + 1 {
                    for c in 0..lwe_params.n {
                        lwes[r * dim + i + c] = out[r * lwe_params.n + c];
                    }
                }
            }

            lwes
        } else {
            let out = self.generate_query_impl(public_seed_idx, dim_log2, packing, index_row);
            let lwes = self.rlwes_to_lwes(&out);
            lwes
        }
    }

    pub fn generate_full_query(
        &self,
        target_idx: usize,
    ) -> (Vec<u32>, AlignedMemory64, AlignedMemory64) {
        // setup
        let db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        let db_cols = 1 << (self.params.db_dim_2 + self.params.poly_len_log2);
        let target_row = target_idx / db_cols;
        let target_col = target_idx % db_cols;
        debug!(
            "Target item: {} ({}, {})",
            target_idx, target_row, target_col
        );

        // generate pub params
        let sk_reg = self.client().get_sk_reg();
        let pack_pub_params = raw_generate_expansion_params(
            self.params,
            &sk_reg,
            self.params.poly_len_log2,
            self.params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(STATIC_SEED_2),
        );
        // let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params) / 2;
        let mut pack_pub_params_row_1s = pack_pub_params.to_vec();
        for i in 0..pack_pub_params.len() {
            pack_pub_params_row_1s[i] =
                pack_pub_params[i].submatrix(1, 0, 1, pack_pub_params[i].cols);
            pack_pub_params_row_1s[i] = condense_matrix(self.params, &pack_pub_params_row_1s[i]);
        }
        let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params_row_1s);
        debug!("pub params size: {} bytes", pub_params_size);
        let pack_pub_params_row_1s_pm = pack_vec_pm(
            self.params,
            1,
            self.params.t_exp_left,
            &pack_pub_params_row_1s,
        );

        // generate query
        let query_row = self.generate_query(SEED_0, self.params.db_dim_1, false, target_row);
        let query_row_last_row: &[u64] = &query_row[self.lwe_params().n * db_rows..];
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

        let query_col = self.generate_query(SEED_1, self.params.db_dim_2, true, target_col);
        let query_col_last_row = &query_col[self.params.poly_len * db_cols..];
        let packed_query_col = pack_query(self.params, query_col_last_row);
        (
            packed_query_row_u32,
            packed_query_col,
            pack_pub_params_row_1s_pm,
        )
    }

    pub fn generate_full_query_simplepir(
        &self,
        target_idx: usize,
    ) -> (AlignedMemory64, AlignedMemory64) {
        // setup
        let db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        let db_cols = self.params.instances * self.params.poly_len;

        let target_row = target_idx / db_cols;
        let target_col = target_idx % db_cols;
        debug!(
            "Target item: {} ({}, {})",
            target_idx, target_row, target_col
        );

        // generate pub params
        let sk_reg = self.client().get_sk_reg();
        let pack_pub_params = raw_generate_expansion_params(
            self.params,
            &sk_reg,
            self.params.poly_len_log2,
            self.params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(STATIC_SEED_2),
        );
        // let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params) / 2;
        let mut pack_pub_params_row_1s = pack_pub_params.to_vec();
        for i in 0..pack_pub_params.len() {
            pack_pub_params_row_1s[i] =
                pack_pub_params[i].submatrix(1, 0, 1, pack_pub_params[i].cols);
            pack_pub_params_row_1s[i] = condense_matrix(self.params, &pack_pub_params_row_1s[i]);
        }
        let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params_row_1s);
        debug!("pub params size: {} bytes", pub_params_size);
        let pack_pub_params_row_1s_pm = pack_vec_pm(
            self.params,
            1,
            self.params.t_exp_left,
            &pack_pub_params_row_1s,
        );

        // generate query
        let query_row = self.generate_query(SEED_0, self.params.db_dim_1, true, target_row);
        assert_eq!(query_row.len(), (self.params.poly_len + 1) * db_rows);
        let query_row_last_row: &[u64] = &query_row[self.params.poly_len * db_rows..];
        assert_eq!(query_row_last_row.len(), db_rows);
        let packed_query_row = pack_query(self.params, query_row_last_row);

        (packed_query_row, pack_pub_params_row_1s_pm)
    }

    fn lwe_params(&self) -> &LWEParams {
        self.lwe_client().lwe_params()
    }

    pub fn decode_response(&self, response: &[u64]) -> Vec<u64> {
        debug!("Decoding response: {:?}", &response[..16]);
        let db_cols = 1 << (self.params.db_dim_2 + self.params.poly_len_log2);

        let sk = self.inner.get_sk_reg().as_slice().to_vec();

        let mut out = Vec::new();
        for col in 0..db_cols {
            let mut sum = 0u128;
            for i in 0..self.params.poly_len {
                let v1 = response[i * db_cols + col];
                let v2 = sk[i];
                sum += v1 as u128 * v2 as u128;
            }

            sum += response[self.params.poly_len * db_cols + col] as u128;

            let result = (sum % self.params.modulus as u128) as u64;
            let result_rescaled = rescale(result, self.params.modulus, self.params.pt_modulus);
            out.push(result_rescaled);
        }

        out
    }

    pub fn client(&self) -> &Client<'a> {
        self.inner
    }
}

// type Seed = [u8; 32]; // <ChaCha20Rng as SeedableRng>::Seed;
// fn generate_secure_random_seed() -> Seed {
//     let mut seed = [0u8; 32];
//     OsRng.fill_bytes(&mut seed);
//     seed
// }

// pub struct YPIRClient {
//     // inner: YClient<'static>,
//     params: &'static Params,
// }

// pub type YPIRQuery = (Vec<u32>, AlignedMemory64, Vec<PolyMatrixNTT<'static>>);
// pub type YPIRSimpleQuery = (AlignedMemory64, Vec<PolyMatrixNTT<'static>>);

// impl YPIRClient {
//     pub fn new(params: &'static Params) -> Self {
//         Self { params }
//     }

//     pub fn generate_full_query(
//         &self,
//         target_idx: usize,
//     ) -> (Vec<u32>, AlignedMemory64, Vec<PolyMatrixNTT<'static>>, Seed) {
//         let client_seed = generate_secure_random_seed();
//         let mut client = Client::init(self.params);
//         client.generate_secret_keys_from_seed(client_seed);
//         let y_client = YClient::new(&mut client, self.params);
//         let query = y_client.generate_full_query(target_idx);
//         (query.0, query.1, query.2, client_seed)
//     }

//     pub fn decode_response(&self, client_seed: Seed, response: &[u64]) -> Vec<u64> {
//         let mut client = Client::init(self.params);
//         client.generate_secret_keys_from_seed(client_seed);
//         let y_client = YClient::new(&mut client, self.params);
//         y_client.decode_response(response)
//     }

//     pub fn params(&self) -> &Params {
//         self.params
//     }
// }

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_lwe() {
        let lwe_params = LWEParams::default();
        let client = LWEClient::new(lwe_params.clone());
        let pt = fastrand::u32(0..lwe_params.pt_modulus as u32);
        let scaled_pt = pt.wrapping_mul(lwe_params.scale_k() as u32);
        let ct = client.encrypt(&mut ChaCha20Rng::from_entropy(), scaled_pt);
        let pt_dec = client.decrypt(&ct);
        let result = rescale(pt_dec as u64, lwe_params.modulus, lwe_params.pt_modulus) as u32;
        assert_eq!(result, pt);
    }
}
