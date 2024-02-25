use std::time::Instant;

use log::debug;

use spiral_rs::{arith::*, gadget::*, ntt::*, params::*, poly::*};

use super::util::*;

fn homomorphic_automorph<'a>(
    params: &'a Params,
    t: usize,
    t_exp: usize,
    ct: &PolyMatrixNTT<'a>,
    pub_param: &PolyMatrixNTT<'a>,
) -> PolyMatrixNTT<'a> {
    assert_eq!(ct.rows, 2);
    assert_eq!(ct.cols, 1);

    let ct_raw = ct.raw();
    let ct_auto = automorph_alloc(&ct_raw, t);

    let mut ginv_ct = PolyMatrixRaw::zero(params, t_exp, 1);
    gadget_invert_rdim(&mut ginv_ct, &ct_auto, 1);
    let mut ginv_ct_ntt = PolyMatrixNTT::zero(params, t_exp, 1);
    for i in 1..t_exp {
        let pol_src = ginv_ct.get_poly(i, 0);
        let pol_dst = ginv_ct_ntt.get_poly_mut(i, 0);
        reduce_copy(params, pol_dst, pol_src);
        ntt_forward(params, pol_dst);
    }
    // let ginv_ct_ntt = ginv_ct.ntt();
    let w_times_ginv_ct = pub_param * &ginv_ct_ntt;

    let mut ct_auto_1 = PolyMatrixRaw::zero(params, 1, 1);
    ct_auto_1
        .data
        .as_mut_slice()
        .copy_from_slice(ct_auto.get_poly(1, 0));
    let ct_auto_1_ntt = ct_auto_1.ntt();

    &ct_auto_1_ntt.pad_top(1) + &w_times_ginv_ct
}

pub fn pack_lwes_inner<'a>(
    params: &'a Params,
    ell: usize,
    start_idx: usize,
    rlwe_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixNTT<'a> {
    assert_eq!(pub_params.len(), params.poly_len_log2);

    if ell == 0 {
        return rlwe_cts[start_idx].clone();
    }

    let step = 1 << (params.poly_len_log2 - ell);
    let even = start_idx;
    let odd = start_idx + step;

    let mut ct_even = pack_lwes_inner(params, ell - 1, even, rlwe_cts, pub_params, y_constants);
    let ct_odd = pack_lwes_inner(params, ell - 1, odd, rlwe_cts, pub_params, y_constants);

    let (y, neg_y) = (&y_constants.0[ell - 1], &y_constants.1[ell - 1]);

    let y_times_ct_odd = scalar_multiply_alloc(&y, &ct_odd);
    let neg_y_times_ct_odd = scalar_multiply_alloc(&neg_y, &ct_odd);

    let mut ct_sum_1 = ct_even.clone();
    add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
    add_into(&mut ct_even, &y_times_ct_odd);

    // let now = Instant::now();
    let ct_sum_1_automorphed = homomorphic_automorph(
        params,
        (1 << ell) + 1,
        params.t_exp_left,
        &ct_sum_1,
        &pub_params[params.poly_len_log2 - 1 - (ell - 1)],
    );
    // debug!("Homomorphic automorph in {} us", now.elapsed().as_micros());

    &ct_even + &ct_sum_1_automorphed
}

fn pack_lwes_inner_non_recursive<'a>(
    params: &'a Params,
    ell: usize,
    _start_idx: usize,
    rlwe_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
    prepared_vals: Option<&[PolyMatrixNTT<'a>]>,
    mut output_prepared_vals: Option<&mut Vec<PolyMatrixNTT<'a>>>,
) -> PolyMatrixNTT<'a> {
    assert!(pub_params.len() == params.poly_len_log2 || pub_params.len() == 0);
    assert_eq!(params.crt_count, 2);

    // let mut working_set = Vec::with_capacity(1 << (ell - 1));
    // let num_out = 1 << (ell - 1);
    // for i in 0..num_out {
    //     let combined = combine(
    //         params,
    //         1,
    //         &rlwe_cts[i],
    //         &rlwe_cts[num_out + i],
    //         pub_params,
    //         y_constants,
    //     );
    //     working_set.push(combined);
    // }

    let mut working_set = rlwe_cts.to_vec();

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 2, 1);

    let mut ct_raw = PolyMatrixRaw::zero(params, 1, 1);
    let mut ct_auto = PolyMatrixRaw::zero(params, 1, 1);
    let mut ginv_ct = PolyMatrixRaw::zero(params, params.t_exp_left, 1);
    let mut ginv_ct_ntt = PolyMatrixNTT::zero(params, params.t_exp_left, 1);
    let mut ct_auto_1_ntt = PolyMatrixNTT::zero(params, 1, 1);
    let mut w_times_ginv_ct = PolyMatrixNTT::zero(params, 2, 1);
    let mut scratch = PolyMatrixNTT::zero(params, 2, 1);
    let scratch_mut_slc = scratch.as_mut_slice();

    let mut total_0 = 0;
    let mut total_1 = 0;
    let mut total_2 = 0;
    let mut total_3 = 0;
    let mut total_4 = 0;

    let mut num_ntts = 0;

    let mut ct_raw_1_auto_ntt;

    for cur_ell in 1..=ell {
        let num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let now = Instant::now();
            let ct_even = &mut first_half[i];
            let ct_odd = &second_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            // if i == 5 {
            //     debug!("neg_y: {:?}", neg_y.as_slice());
            // }

            scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
            scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);

            // if i == 5 && cur_ell == 1 && output_prepared_vals.is_none() {
            //     debug!(
            //         "ct_even[0]: {:?}",
            //         params.crt_compose(ct_even.get_poly(1, 0), 0)
            //     );
            //     debug!(
            //         "ct_odd[0]: {:?}",
            //         params.crt_compose(ct_odd.get_poly(1, 0), 0)
            //     );
            // }

            ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
            add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
            add_into_no_reduce(ct_even, &y_times_ct_odd);
            total_3 += now.elapsed().as_micros();

            {
                let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
                let t = (1 << cur_ell) + 1;
                let t_exp = params.t_exp_left;
                let (cur_ginv_ct_ntt, cur_ct_auto_1_ntt) = if cur_ell == 1
                    && prepared_vals.is_some()
                {
                    let ginv_ct_ntt = &prepared_vals.unwrap()[i];

                    // In this first round, this value is always zero
                    ct_raw_1_auto_ntt = PolyMatrixNTT::zero(params, 1, 1);
                    (ginv_ct_ntt, &ct_raw_1_auto_ntt)
                } else {
                    let now = Instant::now();
                    // let ct_raw = ct.raw();
                    // nb: scratch has 2nd row of ct in uncrtd form,
                    //     ct_raw has only first row
                    from_ntt_scratch(&mut ct_raw, scratch_mut_slc, ct);
                    if cur_ell == 1 {
                        num_ntts += 2;
                    }
                    total_0 += now.elapsed().as_micros();
                    let now = Instant::now();
                    automorph(&mut ct_auto, &ct_raw, t);
                    total_1 += now.elapsed().as_micros();

                    gadget_invert_rdim(&mut ginv_ct, &ct_auto, 1);

                    let skip_first_gadget_dim = true;
                    if skip_first_gadget_dim {
                        for i in 1..t_exp {
                            let pol_src = ginv_ct.get_poly(i, 0);
                            let pol_dst = ginv_ct_ntt.get_poly_mut(i, 0);
                            pol_dst[..params.poly_len].copy_from_slice(pol_src);
                            pol_dst[params.poly_len..].copy_from_slice(pol_src);

                            ntt_forward(params, pol_dst);
                            if cur_ell == 1 {
                                num_ntts += 1;
                            }
                        }
                    } else {
                        to_ntt(&mut ginv_ct_ntt, &ginv_ct);
                        // num_ntts += ginv_ct_ntt.rows * ginv_ct_ntt.cols;
                    }

                    let now = Instant::now();
                    automorph_poly_uncrtd(params, ct_auto_1_ntt.as_mut_slice(), scratch_mut_slc, t);
                    ntt_forward(params, ct_auto_1_ntt.as_mut_slice());
                    // num_ntts += 1;

                    total_4 += now.elapsed().as_micros();

                    (&ginv_ct_ntt, &ct_auto_1_ntt)
                };

                if output_prepared_vals.is_some() {
                    let opv_mut = output_prepared_vals.as_deref_mut();
                    opv_mut.unwrap().push(cur_ginv_ct_ntt.clone());
                    continue;
                }

                let pub_param = &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
                // let ginv_ct_ntt = ginv_ct.ntt();
                // let w_times_ginv_ct = pub_param * &ginv_ct_ntt;
                w_times_ginv_ct.as_mut_slice().fill(0);
                multiply_no_reduce(&mut w_times_ginv_ct, &pub_param, &cur_ginv_ct_ntt, 1);

                // &ct_auto_1_ntt.pad_top(1) + &w_times_ginv_ct
                let now = Instant::now();
                add_into_at_no_reduce(ct_even, &cur_ct_auto_1_ntt, 1, 0);
                add_into(ct_even, &w_times_ginv_ct);
                total_2 += now.elapsed().as_micros();
            };
        }

        if output_prepared_vals.is_some() {
            return PolyMatrixNTT::zero(params, 2, 1);
        }
    }

    if false {
        debug!("num_ntts: {}", num_ntts);
        debug!("total_0: {} us", total_0);
        debug!("total_1: {} us", total_1);
        debug!("total_2: {} us", total_2);
        debug!("total_3: {} us", total_3);
        debug!("total_4: {} us", total_4);
    }

    // let mut res = PolyMatrixNTT::zero(params, 2, 1);
    // // let mut res = working_set[0].clone();
    // add_into(&mut res, &working_set[0]);
    // res

    working_set[0].clone()
}

pub fn precompute_pack<'a>(
    params: &'a Params,
    ell: usize,
    rlwe_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> (PolyMatrixNTT<'a>, Vec<PolyMatrixNTT<'a>>) {
    assert!(pub_params.len() == params.poly_len_log2);
    assert_eq!(params.crt_count, 2);

    let mut working_set = rlwe_cts.to_vec();

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 2, 1);

    let mut ct_raw = PolyMatrixRaw::zero(params, 1, 1);
    let mut ct_auto = PolyMatrixRaw::zero(params, 1, 1);
    let mut ginv_ct = PolyMatrixRaw::zero(params, params.t_exp_left, 1);
    let mut ginv_ct_ntt = PolyMatrixNTT::zero(params, params.t_exp_left, 1);
    let mut ct_auto_1_ntt = PolyMatrixNTT::zero(params, 1, 1);
    let mut w_times_ginv_ct = PolyMatrixNTT::zero(params, 2, 1);
    let mut scratch = PolyMatrixNTT::zero(params, 2, 1);
    let scratch_mut_slc = scratch.as_mut_slice();

    let mut total_0 = 0;
    let mut total_1 = 0;
    let mut total_2 = 0;
    let mut total_3 = 0;
    let mut total_4 = 0;

    let mut num_ntts = 0;

    let mut res = Vec::new();

    for cur_ell in 1..=ell {
        let num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let now = Instant::now();
            let ct_even = &mut first_half[i];
            let ct_odd = &second_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
            scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);

            ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
            add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
            add_into_no_reduce(ct_even, &y_times_ct_odd);
            total_3 += now.elapsed().as_micros();

            {
                let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
                let t = (1 << cur_ell) + 1;
                let t_exp = params.t_exp_left;
                let (cur_ginv_ct_ntt, cur_ct_auto_1_ntt) = {
                    let now = Instant::now();
                    // let ct_raw = ct.raw();

                    // nb: scratch has 2nd row of ct in uncrtd form,
                    //     ct_raw has only first row
                    from_ntt_scratch(&mut ct_raw, scratch_mut_slc, ct);
                    if cur_ell == 1 {
                        num_ntts += 2;
                    }
                    total_0 += now.elapsed().as_micros();
                    let now = Instant::now();
                    automorph(&mut ct_auto, &ct_raw, t);
                    total_1 += now.elapsed().as_micros();

                    gadget_invert_rdim(&mut ginv_ct, &ct_auto, 1);

                    let skip_first_gadget_dim = false;
                    if skip_first_gadget_dim {
                        for i in 1..t_exp {
                            let pol_src = ginv_ct.get_poly(i, 0);
                            let pol_dst = ginv_ct_ntt.get_poly_mut(i, 0);
                            pol_dst[..params.poly_len].copy_from_slice(pol_src);
                            pol_dst[params.poly_len..].copy_from_slice(pol_src);

                            ntt_forward(params, pol_dst);
                            if cur_ell == 1 {
                                num_ntts += 1;
                            }
                        }
                    } else {
                        to_ntt(&mut ginv_ct_ntt, &ginv_ct);
                        // num_ntts += ginv_ct_ntt.rows * ginv_ct_ntt.cols;
                    }

                    let now = Instant::now();
                    automorph_poly_uncrtd(params, ct_auto_1_ntt.as_mut_slice(), scratch_mut_slc, t);
                    ntt_forward(params, ct_auto_1_ntt.as_mut_slice());
                    // num_ntts += 1;

                    total_4 += now.elapsed().as_micros();

                    (&ginv_ct_ntt, &ct_auto_1_ntt)
                };

                // println!(
                //     "ct_auto_1_ntt.raw(): {:?}",
                //     &ct_auto_1_ntt.raw().as_slice()[..30]
                // );

                res.push(cur_ginv_ct_ntt.clone());

                let pub_param = &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
                // let ginv_ct_ntt = ginv_ct.ntt();
                // let w_times_ginv_ct = pub_param * &ginv_ct_ntt;
                w_times_ginv_ct.as_mut_slice().fill(0);
                multiply_no_reduce(&mut w_times_ginv_ct, &pub_param, &cur_ginv_ct_ntt, 0);

                // &ct_auto_1_ntt.pad_top(1) + &w_times_ginv_ct
                let now = Instant::now();
                add_into_at_no_reduce(ct_even, &cur_ct_auto_1_ntt, 1, 0);
                add_into(ct_even, &w_times_ginv_ct);
                total_2 += now.elapsed().as_micros();
            };
        }
    }

    if false {
        debug!("num_ntts: {}", num_ntts);
        debug!("total_0: {} us", total_0);
        debug!("total_1: {} us", total_1);
        debug!("total_2: {} us", total_2);
        debug!("total_3: {} us", total_3);
        debug!("total_4: {} us", total_4);
    }

    (working_set[0].clone(), res)
}

// pub fn client_precomp_w<'a>(
//     params: &'a Params,
//     ell: usize,
//     pub_params: &[PolyMatrixNTT<'a>],
// ) -> Vec<PolyMatrixNTT<'a>> {

// }

pub fn pack_using_precomp_vals<'a>(
    params: &'a Params,
    ell: usize,
    pub_params: &[PolyMatrixNTT<'a>],
    b_values: &[u64],
    precomp_res: &PolyMatrixNTT<'a>,
    precomp_vals: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixNTT<'a> {
    let mut working_set = vec![PolyMatrixNTT::zero(params, 2, 1); 1 << ell];

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 2, 1);

    let mut ct_raw = PolyMatrixRaw::zero(params, 1, 1);
    let mut ct_auto_1_ntt = PolyMatrixNTT::zero(params, 1, 1);
    let mut scratch = PolyMatrixNTT::zero(params, 2, 1);
    let scratch_mut_slc = scratch.as_mut_slice();

    let tables = generate_automorph_tables_brute_force(&params); // TODO: move out

    let mut idx_precomp = 0;
    for cur_ell in 1..=ell {
        let num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let ct_even = &mut first_half[i];
            let ct_odd = &second_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
            scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);

            ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
            add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
            add_into_no_reduce(ct_even, &y_times_ct_odd);

            // --

            let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
            let ct_row_1 = ct.submatrix(1, 0, 1, 1);

            let t = (1 << cur_ell) + 1;
            let t_exp = params.t_exp_left;

            // nb: scratch has 2nd row of ct in uncrtd form,
            //     ct_raw has only first row
            // from_ntt_scratch(&mut ct_raw, scratch_mut_slc, ct);

            let cur_ginv_ct_ntt = &precomp_vals[idx_precomp];
            idx_precomp += 1;

            let pub_param = &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
            let w_times_ginv_ct = pub_param * cur_ginv_ct_ntt;

            // automorph_poly_uncrtd(params, ct_auto_1_ntt.as_mut_slice(), scratch_mut_slc, t);
            let ct_auto_1_ntt = apply_automorph_ntt(params, &tables, &ct_row_1, t);
            // ntt_forward(params, ct_auto_1_ntt.as_mut_slice());

            // zero out first row
            // for col in 0..w_times_ginv_ct.cols {
            //     w_times_ginv_ct.get_poly_mut(0, col).fill(0);
            // }

            add_into_at_no_reduce(ct_even, &ct_auto_1_ntt, 1, 0);
            add_into(ct_even, &w_times_ginv_ct);
        }
    }

    assert_eq!(idx_precomp, precomp_vals.len());

    let resulting_row_1 = &working_set[0].get_poly(1, 0).to_vec();

    let mut res = precomp_res.clone();
    // {
    //     let r = res.raw();
    //     println!("res row 0: {:?}", &r.get_poly(0, 0)[..30]);
    //     println!("res row 1: {:?}", &r.get_poly(1, 0)[..30]);
    // }
    res.get_poly_mut(1, 0).copy_from_slice(resulting_row_1);

    println!(
        "precomp_res     row 1: {:?}",
        &precomp_res.raw().get_poly(1, 0)[..30]
    );
    println!(
        "resulting_row_1 row 1: {:?}",
        &working_set[0].raw().get_poly(1, 0)[..30]
    );

    // let mut res = precomp_res.clone();

    let mut out_raw = res.raw();
    for z in 0..params.poly_len {
        let val = barrett_reduction_u128(params, b_values[z] as u128 * params.poly_len as u128);
        out_raw.get_poly_mut(1, 0)[z] = barrett_u64(params, out_raw.get_poly(1, 0)[z] + val);
    }
    let out = out_raw.ntt();

    // for z in 0..params.poly_len {
    //     let b_value = b_values[z];
    //     let val = barrett_u64(params, res.get_poly(1, 0)[z] + b_value);
    //     res.get_poly_mut(1, 0)[z] = val;
    // }

    out
}

pub fn combine<'a>(
    params: &'a Params,
    cur_ell: usize,
    ct_even: &mut PolyMatrixNTT<'a>,
    ct_odd: &PolyMatrixNTT<'a>,
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) {
    let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

    let y_times_ct_odd = scalar_multiply_alloc(&y, &ct_odd);
    let neg_y_times_ct_odd = scalar_multiply_alloc(&neg_y, &ct_odd);

    let mut ct_sum_1 = ct_even.clone();
    add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
    add_into(ct_even, &y_times_ct_odd);

    let ct_sum_1_automorphed = homomorphic_automorph(
        params,
        (1 << cur_ell) + 1,
        params.t_exp_left,
        &ct_sum_1,
        &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)],
    );

    add_into(ct_even, &ct_sum_1_automorphed);
}

pub fn prep_pack_lwes<'a>(
    params: &'a Params,
    lwe_cts: &[u64],
    cols_to_do: usize,
) -> Vec<PolyMatrixNTT<'a>> {
    let lwe_cts_size = params.poly_len * (params.poly_len + 1);
    assert_eq!(lwe_cts.len(), lwe_cts_size);

    assert!(cols_to_do == params.poly_len);

    let mut rlwe_cts = Vec::new();
    for i in 0..cols_to_do {
        let mut rlwe_ct = PolyMatrixRaw::zero(params, 2, 1);

        // 'a' vector
        // put this in negacyclic order
        let mut poly = Vec::new();
        for j in 0..params.poly_len {
            poly.push(lwe_cts[j * params.poly_len + i])
        }
        let nega = negacyclic_perm(&poly, 0, params.modulus);

        for j in 0..params.poly_len {
            rlwe_ct.get_poly_mut(0, 0)[j] = nega[j];
        }
        // 'b' scalar (skip)

        rlwe_cts.push(rlwe_ct.ntt());
    }

    rlwe_cts
}

pub fn prep_pack_many_lwes<'a>(
    params: &'a Params,
    lwe_cts: &[u64],
    num_rlwe_outputs: usize,
) -> Vec<Vec<PolyMatrixNTT<'a>>> {
    let lwe_cts_size = (params.poly_len + 1) * (num_rlwe_outputs * params.poly_len);
    assert_eq!(lwe_cts.len(), lwe_cts_size);

    let mut vecs = Vec::new();
    for i in 0..num_rlwe_outputs {
        let mut v = Vec::new();
        for j in 0..params.poly_len + 1 {
            v.extend(
                &lwe_cts[j * (num_rlwe_outputs * params.poly_len) + i * params.poly_len..]
                    [..params.poly_len],
            );
        }
        vecs.push(v);
    }

    let mut res = Vec::new();
    for i in 0..num_rlwe_outputs {
        res.push(prep_pack_lwes(params, &vecs[i], params.poly_len));
    }

    res
}

pub fn prepare_packed_vals_pack_lwes<'a>(
    params: &'a Params,
    preped_rlwe_cts: &[PolyMatrixNTT<'a>],
    _cols_to_do: usize,
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<PolyMatrixNTT<'a>> {
    let now = Instant::now();
    let mut output_preped_packed_vals = Vec::new();
    pack_lwes_inner_non_recursive(
        params,
        params.poly_len_log2,
        0,
        &preped_rlwe_cts,
        &[],
        y_constants,
        None,
        Some(&mut output_preped_packed_vals),
    );
    debug!("prepack: {} us", now.elapsed().as_micros());
    output_preped_packed_vals
}

/// Returns the `prep_packed_vals` value.
pub fn prep_pack_many_lwes_packed_vals<'a>(
    params: &'a Params,
    prep_rlwe_cts: &[Vec<PolyMatrixNTT<'a>>],
    num_rlwe_outputs: usize,
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<Vec<PolyMatrixNTT<'a>>> {
    assert_eq!(prep_rlwe_cts.len(), num_rlwe_outputs);
    assert_eq!(prep_rlwe_cts[0].len(), params.poly_len);

    let mut res = Vec::new();
    for i in 0..num_rlwe_outputs {
        res.push(prepare_packed_vals_pack_lwes(
            params,
            &prep_rlwe_cts[i],
            params.poly_len,
            y_constants,
        ));
    }

    res
}

pub fn pack_lwes<'a>(
    params: &'a Params,
    b_values: &[u64],
    preped_rlwe_cts: &[PolyMatrixNTT<'a>],
    preped_packed_vals: &[PolyMatrixNTT<'a>],
    cols_to_do: usize,
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixNTT<'a> {
    assert_eq!(preped_rlwe_cts.len(), cols_to_do);
    assert_eq!(cols_to_do, params.poly_len);
    assert_eq!(b_values.len(), params.poly_len);

    // let mut rlwe_cts = Vec::new();
    // for i in 0..cols_to_do {
    //     let mut rlwe_ct = preped_rlwe_cts[i].clone();

    //     // 'b' scalar
    //     let b_val = 0; //b_values[i];

    //     assert_eq!(params.crt_count, 2);
    //     let parts = rlwe_ct.get_poly_mut(1, 0).split_at_mut(params.poly_len);
    //     // assert!(parts.0.iter().all(|&x| x == 0));
    //     // assert!(parts.1.iter().all(|&x| x == 0));
    //     parts.0.fill(b_val % params.moduli[0]);
    //     parts.1.fill(b_val % params.moduli[1]); // TODO: use barrett reduction

    //     rlwe_cts.push(rlwe_ct);
    // }

    // let now = Instant::now();
    // let mut output_preped_packed_vals = Vec::new();
    // pack_lwes_inner_non_recursive(
    //     params,
    //     params.poly_len_log2,
    //     0,
    //     &preped_rlwe_cts,
    //     &[],
    //     y_constants,
    //     None,
    //     Some(&mut output_preped_packed_vals),
    // );
    // debug!("prepack: {} us", now.elapsed().as_micros());
    let now = Instant::now();
    let preped_packed_val_opt = if preped_packed_vals.len() == 0 {
        None
    } else {
        Some(preped_packed_vals)
    };
    let out = pack_lwes_inner_non_recursive(
        params,
        params.poly_len_log2,
        0,
        &preped_rlwe_cts,
        pub_params,
        y_constants,
        preped_packed_val_opt,
        None,
    );
    let mut out_raw = out.raw();
    for z in 0..params.poly_len {
        let val = barrett_reduction_u128(params, b_values[z] as u128 * params.poly_len as u128);
        out_raw.get_poly_mut(1, 0)[z] = barrett_u64(params, out_raw.get_poly(1, 0)[z] + val);
    }
    let res = out_raw.ntt();
    debug!("True pack took {} us", now.elapsed().as_micros());

    res
}

pub fn pack_many_lwes<'a>(
    params: &'a Params,
    prep_rlwe_cts: &[Vec<PolyMatrixNTT<'a>>],
    prep_packed_vals: &[Vec<PolyMatrixNTT<'a>>],
    b_values: &[u64],
    num_rlwe_outputs: usize,
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<PolyMatrixNTT<'a>> {
    assert_eq!(prep_rlwe_cts.len(), num_rlwe_outputs);
    assert_eq!(prep_rlwe_cts[0].len(), params.poly_len);
    assert_eq!(b_values.len(), num_rlwe_outputs * params.poly_len);

    let mut res = Vec::new();
    for i in 0..num_rlwe_outputs {
        // The last RLWE output cannot be 'prepared' because it changes for each query
        let empty = Vec::new();
        let prep_packed_val = if i == num_rlwe_outputs - 1 {
            &empty
        } else {
            &prep_packed_vals[i]
        };
        res.push(pack_lwes(
            params,
            &b_values[i * params.poly_len..][..params.poly_len],
            &prep_rlwe_cts[i],
            &prep_packed_val,
            params.poly_len,
            pub_params,
            y_constants,
        ));
    }

    res
}

fn swap_midpoint<T>(a: &mut [T]) {
    let len = a.len();
    let (a, b) = a.split_at_mut(len / 2);
    a.swap_with_slice(b);
}

pub fn produce_table(poly_len: usize, chunk_size: usize) -> Vec<usize> {
    let mut cur = (0..poly_len).collect::<Vec<_>>();

    let outer_chunk_size = poly_len / (chunk_size / 2);
    println!("outer_chunk_size {}", outer_chunk_size);

    let mut do_it = true;
    for outer_chunk in cur.chunks_mut(outer_chunk_size) {
        if !do_it {
            do_it = true;
            continue;
        }
        do_it = false;

        for chunk in outer_chunk.chunks_mut(chunk_size) {
            let mut offs = 0;
            let mut to_add_to_offs = (chunk_size / 2).min(chunk.len() / 2); // weird hack
            while to_add_to_offs > 0 {
                swap_midpoint(&mut chunk[offs..]);
                offs += to_add_to_offs;
                to_add_to_offs /= 2;
            }
        }
    }

    cur
}

pub fn automorph_ntt_tables(poly_len: usize, log2_poly_len: usize) -> Vec<Vec<usize>> {
    let mut tables = Vec::new();
    for i in 0..log2_poly_len {
        let chunk_size = 1 << i;
        println!("table {}", i);
        let table = produce_table(poly_len, 2 * chunk_size);
        println!("table {:?}", &table.as_slice()[..32]);
        tables.push(table);
    }

    tables
}

pub fn generate_automorph_tables_brute_force(params: &Params) -> Vec<Vec<usize>> {
    let mut tables = Vec::new();
    for i in (1..=params.poly_len_log2).rev() {
        let t = (1 << i) + 1;

        let poly = PolyMatrixRaw::random(&params, 1, 1);
        let poly_ntt = poly.ntt();

        let poly_auto = automorph_alloc(&poly, t);
        let poly_auto_ntt = poly_auto.ntt();

        let pol_orig = (&poly_ntt.get_poly(0, 0)[..params.poly_len]).to_vec();
        let pol_auto = (&poly_auto_ntt.get_poly(0, 0)[..params.poly_len]).to_vec();

        let mut table_candidate = vec![0usize; params.poly_len];
        for i in 0..params.poly_len {
            table_candidate[i] = pol_orig.iter().position(|&x| x == pol_auto[i]).unwrap();
        }
        tables.push(table_candidate);
    }
    tables
}

pub fn apply_automorph_ntt_raw<'a>(
    params: &Params,
    poly: &[u64],
    out: &mut [u64],
    t: usize,
    tables: &[Vec<usize>],
) {
    let poly_len = params.poly_len;
    // table_idx = log2(poly_len / (t - 1))
    let table_idx = (poly_len / (t - 1)).trailing_zeros() as usize;
    let table = &tables[table_idx];

    for i in 0..poly_len {
        out[i] = poly[table[i]];
    }
}

pub fn apply_automorph_ntt<'a>(
    params: &'a Params,
    tables: &[Vec<usize>],
    mat: &PolyMatrixNTT<'a>,
    t: usize,
) -> PolyMatrixNTT<'a> {
    // run apply_automorph_ntt on each poly in the matrix
    let mut res = PolyMatrixNTT::zero(params, mat.rows, mat.cols);
    for i in 0..mat.rows {
        for j in 0..mat.cols {
            let poly = mat.get_poly(i, j);
            let res_poly = res.get_poly_mut(i, j);
            for (chunk, res_chunk) in poly
                .chunks_exact(params.poly_len)
                .zip(res_poly.chunks_exact_mut(params.poly_len))
            {
                apply_automorph_ntt_raw(params, chunk, res_chunk, t, tables);
            }
        }
    }
    res
}

#[cfg(test)]
mod test {
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use spiral_rs::{client::Client, number_theory::invert_uint_mod, util::get_test_params};

    use crate::{
        client::raw_generate_expansion_params, params::params_for_scenario,
        server::generate_y_constants,
    };

    use super::*;

    #[test]
    fn test_packing() {
        let params = get_test_params();
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let y_constants = generate_y_constants(&params);

        let pack_seed = [1u8; 32];
        let cts_seed = [2u8; 32];
        let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

        let pack_pub_params = raw_generate_expansion_params(
            &params,
            client.get_sk_reg(),
            params.poly_len_log2,
            params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(pack_seed),
        );

        // generate poly_len ciphertexts
        let mut v_ct = Vec::new();
        let mut b_values = Vec::new();
        for i in 0..params.poly_len {
            let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
            let val = i as u64 % params.pt_modulus;
            let scale_k = params.modulus / params.pt_modulus;
            let mod_inv = invert_uint_mod(params.poly_len as u64, params.modulus).unwrap();
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );
            let mut ct_raw = ct.raw();

            // get the b value
            b_values.push(ct_raw.get_poly(1, 0)[0]);

            // zero out all of the second poly
            ct_raw.get_poly_mut(1, 0).fill(0);
            v_ct.push(ct_raw.ntt());
        }

        let now = Instant::now();
        let packed = pack_lwes(
            &params,
            &b_values,
            &v_ct,
            &[],
            params.poly_len,
            &pack_pub_params,
            &y_constants,
        );
        println!("Packing took {} us", now.elapsed().as_micros());

        let packed_raw = packed.raw();
        println!("packed_0: {:?}", &packed_raw.get_poly(0, 0)[..10]);
        assert_eq!(packed_raw.get_poly(0, 0)[0], 48718429063985342u64);

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed);
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            rescaled.data[i] = rescale(dec_raw.data[i], params.modulus, params.pt_modulus);
        }

        println!("rescaled: {:?}", &rescaled.as_slice()[..50]);
        let mut gold = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            gold.data[i] = i as u64 % params.pt_modulus;
        }
        assert_eq!(rescaled.as_slice(), gold.as_slice());
    }

    #[test]
    fn test_precompute_packing() {
        let params = params_for_scenario(1 << 30, 1);
        println!("modulus: {}", params.modulus);
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let y_constants = generate_y_constants(&params);

        let pack_seed = [1u8; 32];
        let cts_seed = [2u8; 32];
        let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

        let pack_pub_params = raw_generate_expansion_params(
            &params,
            client.get_sk_reg(),
            params.poly_len_log2,
            params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(pack_seed),
        );

        // generate poly_len ciphertexts
        let mut v_ct = Vec::new();
        let mut b_values = Vec::new();
        for i in 0..params.poly_len {
            let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
            let val = i as u64 % params.pt_modulus;
            let scale_k = params.modulus / params.pt_modulus;
            let mod_inv = invert_uint_mod(params.poly_len as u64, params.modulus).unwrap();
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );
            let mut ct_raw = ct.raw();

            // get the b value
            b_values.push(ct_raw.get_poly(1, 0)[0]);

            // zero out all of the second poly
            ct_raw.get_poly_mut(1, 0).fill(0);
            v_ct.push(ct_raw.ntt());
        }

        let now = Instant::now();
        let (precomp_res, precomp_vals) = precompute_pack(
            &params,
            params.poly_len_log2,
            &v_ct,
            &pack_pub_params,
            &y_constants,
        );
        println!(
            "Precomputing for packing took {} us",
            now.elapsed().as_micros()
        );

        let now = Instant::now();
        let packed = pack_using_precomp_vals(
            &params,
            params.poly_len_log2,
            &pack_pub_params,
            &b_values,
            &precomp_res,
            &precomp_vals,
            &y_constants,
        );
        println!("Packing took {} us", now.elapsed().as_micros());

        let packed_raw = packed.raw();
        println!("packed_0: {:?}", &packed_raw.get_poly(0, 0)[..10]);
        // assert_eq!(packed_raw.get_poly(0, 0)[0], 17210016925609510u64);

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed);
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            rescaled.data[i] = rescale(dec_raw.data[i], params.modulus, params.pt_modulus);
        }

        println!("rescaled: {:?}", &rescaled.as_slice()[..50]);
        let mut gold = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            gold.data[i] = i as u64 % params.pt_modulus;
        }
        assert_eq!(rescaled.as_slice(), gold.as_slice());
    }

    #[test]
    fn test_automorph_tables() {
        let params = get_test_params();

        let now = Instant::now();
        let tables = generate_automorph_tables_brute_force(&params);
        println!("Generating tables took {} us", now.elapsed().as_micros());

        for t in [3, 9, 17, 33, 65, 129, 257, 513, 1025, 2049] {
            let poly = PolyMatrixRaw::random(&params, 1, 1);
            let poly_ntt = poly.ntt();

            let poly_auto = automorph_alloc(&poly, t);
            let poly_auto_ntt = poly_auto.ntt();
            let poly_auto_ntt_using_tables = apply_automorph_ntt(&params, &tables, &poly_ntt, t);

            assert_eq!(
                &poly_auto_ntt.as_slice(),
                &poly_auto_ntt_using_tables.as_slice()
            );
        }
    }
}
