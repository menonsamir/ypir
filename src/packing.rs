use std::{arch::x86_64::*, time::Instant};

use log::debug;

use spiral_rs::{arith::*, gadget::*, ntt::*, params::*, poly::*};

use crate::server::Precomp;

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
            fast_add_into_no_reduce(ct_even, &y_times_ct_odd);
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
    fake_pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> (PolyMatrixNTT<'a>, Vec<PolyMatrixNTT<'a>>, Vec<Vec<usize>>) {
    assert!(fake_pub_params.len() == params.poly_len_log2);
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
            fast_add_into_no_reduce(ct_even, &y_times_ct_odd);
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

                res.push(condense_matrix(params, cur_ginv_ct_ntt));

                let pub_param = &fake_pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
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

    let tables = generate_automorph_tables_brute_force(&params);

    (working_set[0].clone(), res, tables)
}

pub fn pack_using_precomp_vals<'a>(
    params: &'a Params,
    ell: usize,
    pub_params: &[PolyMatrixNTT<'a>],
    b_values: &[u64],
    precomp_res: &PolyMatrixNTT<'a>,
    precomp_vals: &[PolyMatrixNTT<'a>],
    precomp_tables: &[Vec<usize>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixNTT<'a> {
    // let now = Instant::now();
    // let mut working_set = vec![PolyMatrixNTT::zero(params, 1, 1); 1 << ell];
    let mut working_set = Vec::with_capacity(1 << (ell - 1));
    for _ in 0..(1 << (ell - 1)) {
        working_set.push(PolyMatrixNTT::zero(params, 1, 1));
    }

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 1, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 1, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 1, 1);
    let mut w_times_ginv_ct = PolyMatrixNTT::zero(params, 1, 1);

    // println!("time_-1: {} us", now.elapsed().as_micros());

    let mut time_0 = 0;
    let mut time_1 = 0;
    let mut time_2 = 0;
    let mut time_3 = 0;
    let mut time_4 = 0;

    let mut idx_precomp = 0;
    let mut num_muls = 0;
    for cur_ell in 1..=ell {
        let mut num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        if num_in == params.poly_len {
            num_in = num_out;
        }

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let now = Instant::now();
            let ct_even = &mut first_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            if cur_ell > 1 {
                let ct_odd = &mut second_half[i];
                scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
                scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);
            }

            time_0 += now.elapsed().as_micros();

            let now = Instant::now();
            if cur_ell > 1 {
                ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
                fast_add_into_no_reduce(&mut ct_sum_1, &neg_y_times_ct_odd);
                fast_add_into_no_reduce(ct_even, &y_times_ct_odd);
            }
            time_1 += now.elapsed().as_micros();

            // --

            let now = Instant::now();
            let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
            let t = (1 << cur_ell) + 1;

            let cur_ginv_ct_ntt = &precomp_vals[idx_precomp];
            idx_precomp += 1;

            let w = &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
            // let w = pub_param.submatrix(1, 0, 1, pub_param.cols);
            // let w_times_ginv_ct = &w * cur_ginv_ct_ntt;
            // multiply(&mut w_times_ginv_ct, &w, &cur_ginv_ct_ntt);
            // w_times_ginv_ct.as_mut_slice().fill(0);
            fast_multiply_no_reduce(params, &mut w_times_ginv_ct, &w, &cur_ginv_ct_ntt, 0);
            num_muls += 1;
            time_2 += now.elapsed().as_micros();

            if cur_ell > 1 {
                let now = Instant::now();
                apply_automorph_ntt(params, &precomp_tables, &ct, ct_even, t);

                // fast_add_into_no_reduce(ct_even, &ct_auto_1_ntt);
                time_3 += now.elapsed().as_micros();
                let now = Instant::now();

                // second condition prevents overflow
                if i < num_out / 2 && ((cur_ell - 1) % 5 != 0) {
                    fast_add_into_no_reduce(ct_even, &w_times_ginv_ct);
                } else {
                    // reduction right before or after addition is much faster than at multiplication time
                    fast_add_into(ct_even, &w_times_ginv_ct);
                }
                time_4 += now.elapsed().as_micros();
            } else {
                let now = Instant::now();
                if i < num_out / 2 {
                    fast_add_into_no_reduce(ct_even, &w_times_ginv_ct);
                } else {
                    fast_add_into(ct_even, &w_times_ginv_ct);
                }
                time_4 += now.elapsed().as_micros();
            }
        }
    }
    // let now = Instant::now();

    if false {
        println!("time_0: {} us", time_0);
        println!("time_1: {} us", time_1);
        println!("time_2: {} us", time_2);
        println!("time_3: {} us", time_3);
        println!("time_4: {} us", time_4);
        println!("idx_precomp: {}", idx_precomp);
        println!("num_muls: {}", num_muls);
    }

    assert_eq!(idx_precomp, precomp_vals.len());

    let mut resulting_row_1 = working_set[0].clone();
    fast_reduce(&mut resulting_row_1);

    let resulting_row_1 = resulting_row_1.as_slice();

    let mut res = precomp_res.clone();
    // {
    //     let r = res.raw();
    //     println!("res row 0: {:?}", &r.get_poly(0, 0)[..30]);
    //     println!("res row 1: {:?}", &r.get_poly(1, 0)[..30]);
    // }
    res.get_poly_mut(1, 0).copy_from_slice(resulting_row_1);

    // println!(
    //     "precomp_res     row 1: {:?}",
    //     &precomp_res.raw().get_poly(1, 0)[..30]
    // );
    // println!(
    //     "resulting_row_1 row 1: {:?}",
    //     &working_set[0].raw().get_poly(1, 0)[..30]
    // );

    // let mut res = precomp_res.clone();

    let mut out_raw = res.raw();
    for z in 0..params.poly_len {
        let val = barrett_reduction_u128(params, b_values[z] as u128 * params.poly_len as u128);
        let idx = params.poly_len + z;
        out_raw.data[idx] += val;
        if out_raw.data[idx] >= params.modulus {
            out_raw.data[idx] -= params.modulus;
        }
    }
    let out = out_raw.ntt();
    // println!("time_5: {} us", now.elapsed().as_micros());

    // for z in 0..params.poly_len {
    //     let b_value = b_values[z];
    //     let val = barrett_u64(params, res.get_poly(1, 0)[z] + b_value);
    //     res.get_poly_mut(1, 0)[z] = val;
    // }

    out
}

pub fn pack_single_lwe<'a>(
    params: &'a Params,
    pub_params: &[PolyMatrixNTT<'a>],
    lwe_ct: &PolyMatrixNTT<'a>,
) -> PolyMatrixNTT<'a> {
    // computing:
    // r0 = f
    // r1 = r0 + automorph(r0, ts[0])
    // r2 = r1 + automorph(r1, ts[1])
    // ...
    // r_\log d = ...

    let mut cur_r = lwe_ct.clone();
    for i in 0..params.poly_len_log2 {
        let t = (params.poly_len / (1 << i)) + 1;
        let pub_param = &pub_params[i];
        let tau_of_r = homomorphic_automorph(params, t, params.t_exp_left, &cur_r, pub_param);
        add_into(&mut cur_r, &tau_of_r);
    }
    cur_r
}

// pub fn fast_scalar_multiply_avx(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT, b: &PolyMatrixNTT) {
//     assert_eq!(a.rows, 1);
//     assert_eq!(a.cols, 1);
//     assert_eq!(b.rows, 1);
//     assert_eq!(b.cols, 1);
//     assert_eq!(res.rows, 1);
//     assert_eq!(res.cols, 1);

//     // this is totally custom for cur_ell == 2

//     let params = res.params;
//     // let pol2 = a.get_poly(0, 0);
//     for i in 0..b.rows {
//         for j in 0..b.cols {
//             let res_slc = res.get_poly_mut(i, j);
//             let a_slc = a.get_poly(0, 0);
//             let b_slc = b.get_poly(i, j);
//             unsafe {
//                 let a_ptr = a_slc.as_ptr();
//                 let b_ptr = b_slc.as_ptr();
//                 let res_ptr = res_slc.as_mut_ptr();

//                 for m in 0..8 {
//                     // all the values in the NTT form of the scalar polynomial are the same
//                     let a_val = *a_ptr.add(m * 512);
//                     let x = _mm256_set1_epi64x(a_val as i64);
//                     for z in (0..512).step_by(4) {
//                         let p_y = b_ptr.add(m * 512 + z);
//                         let y = _mm256_load_si256(p_y as *const _);
//                         let product = _mm256_mul_epu32(x, y);

//                         let p_z = res_ptr.add(m * 512 + z);
//                         _mm256_store_si256(p_z as *mut _, product);
//                     }
//                 }
//             }
//         }
//     }
// }

pub fn fast_barrett_raw_u64(input: u64, const_ratio_1: u64, modulus: u64) -> u64 {
    let tmp = (((input as u128) * (const_ratio_1 as u128)) >> 64) as u64;

    // Barrett subtraction
    let res = input - tmp * modulus;

    res
}

pub fn fast_add_into(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT) {
    assert!(res.rows == a.rows);
    assert!(res.cols == a.cols);

    let params = res.params;
    for i in 0..res.rows {
        for j in 0..res.cols {
            let res_poly = res.get_poly_mut(i, j);
            let a_poly = a.get_poly(i, j);
            for c in 0..params.crt_count {
                for i in 0..params.poly_len {
                    let idx = c * params.poly_len + i;
                    unsafe {
                        let p_res = res_poly.as_mut_ptr().add(idx);
                        let p_a = a_poly.as_ptr().add(idx);
                        let val = *p_res + *p_a;
                        let reduced =
                            fast_barrett_raw_u64(val, params.barrett_cr_1[c], params.moduli[c]);
                        *p_res = reduced;
                    }
                }
            }
        }
    }
}

pub fn fast_multiply_no_reduce(
    params: &Params,
    res: &mut PolyMatrixNTT,
    a: &PolyMatrixNTT,
    b: &PolyMatrixNTT,
    _start_inner_dim: usize,
) {
    assert_eq!(res.rows, a.rows);
    assert_eq!(res.cols, b.cols);
    assert_eq!(res.rows, 1);
    assert_eq!(res.cols, 1);

    assert_eq!(a.cols, b.rows);
    assert_eq!(params.crt_count * params.poly_len, 2 * 2048);

    unsafe {
        let a_ptr = a.as_slice().as_ptr();
        let b_ptr = b.as_slice().as_ptr();
        let res_ptr = res.as_mut_slice().as_mut_ptr();
        let pol_sz = params.poly_len;

        for idx in (0..pol_sz).step_by(8) {
            let mut sum_lo = _mm512_setzero_si512();
            let mut sum_hi = _mm512_setzero_si512();
            for k in 0..a.cols {
                let p_x = a_ptr.add(k * 2 * pol_sz + idx);
                let p_y = b_ptr.add(k * 2 * pol_sz + idx);

                let x = _mm512_load_si512(p_x as *const _);
                let x_lo = x;
                let x_hi = _mm512_srli_epi64(x, 32);
                let y = _mm512_load_si512(p_y as *const _);
                let y_lo = y;
                let y_hi = _mm512_srli_epi64(y, 32);

                let product_lo = _mm512_mul_epu32(x_lo, y_lo);
                let product_hi = _mm512_mul_epu32(x_hi, y_hi);

                sum_lo = _mm512_add_epi64(sum_lo, product_lo);
                sum_hi = _mm512_add_epi64(sum_hi, product_hi);
            }

            let p_z = res_ptr.add(idx);
            _mm512_store_si512(p_z as *mut _, sum_lo);
            let p_z = res_ptr.add(pol_sz + idx);
            _mm512_store_si512(p_z as *mut _, sum_hi);
        }
    }
}

pub fn condense_matrix<'a>(params: &'a Params, a: &PolyMatrixNTT<'a>) -> PolyMatrixNTT<'a> {
    let mut res = PolyMatrixNTT::zero(params, a.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            let res_poly = &mut res.get_poly_mut(i, j);
            let a_poly = a.get_poly(i, j);
            for z in 0..params.poly_len {
                res_poly[z] = a_poly[z] | (a_poly[z + params.poly_len] << 32);
            }
        }
    }
    res
}

pub fn uncondense_matrix<'a>(params: &'a Params, a: &PolyMatrixNTT<'a>) -> PolyMatrixNTT<'a> {
    let mut res = PolyMatrixNTT::zero(params, a.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            let res_poly = &mut res.get_poly_mut(i, j);
            let a_poly = a.get_poly(i, j);
            for z in 0..params.poly_len {
                res_poly[z] = a_poly[z] & ((1u64 << 32) - 1);
                res_poly[z + params.poly_len] = a_poly[z] >> 32;
            }
        }
    }
    res
}

pub fn multiply_add_poly_avx(_params: &Params, res: &mut [u64], a: &[u64], b: &[u64]) {
    unsafe {
        let a_ptr = a.as_ptr();
        let b_ptr = b.as_ptr();
        let res_ptr = res.as_mut_ptr();
        for i in (0..res.len()).step_by(8) {
            let p_x = a_ptr.add(i);
            let p_y = b_ptr.add(i);
            let p_z = res_ptr.add(i);

            let x = _mm512_load_si512(p_x as *const _);
            let y = _mm512_load_si512(p_y as *const _);
            let z = _mm512_load_si512(p_z as *const _);

            let product = _mm512_mul_epu32(x, y);
            let out = _mm512_add_epi64(z, product);

            _mm512_store_si512(p_z as *mut _, out);
        }
    }
}

pub fn multiply_poly_avx(_params: &Params, res: &mut [u64], a: &[u64], b: &[u64]) {
    unsafe {
        let a_ptr = a.as_ptr();
        let b_ptr = b.as_ptr();
        let res_ptr = res.as_mut_ptr();

        for i in (0..res.len()).step_by(8) {
            let p_x = a_ptr.add(i);
            let p_y = b_ptr.add(i);
            let p_z = res_ptr.add(i);

            let x = _mm512_load_si512(p_x as *const _);
            let y = _mm512_load_si512(p_y as *const _);

            let product = _mm512_mul_epu32(x, y);

            _mm512_store_si512(p_z as *mut _, product);
        }
    }
}

pub fn fast_add_into_no_reduce(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT) {
    assert!(res.rows == a.rows);
    assert!(res.cols == a.cols);

    let a_slc = a.as_slice();
    let res_slc = res.as_mut_slice();
    for (res_chunk, a_chunk) in res_slc.chunks_exact_mut(8).zip(a_slc.chunks_exact(8)) {
        for i in 0..8 {
            res_chunk[i] += a_chunk[i];
        }
    }
}

pub fn fast_reduce(res: &mut PolyMatrixNTT) {
    let params = res.params;
    let res_slc = res.as_mut_slice();
    for m in 0..params.crt_count {
        for i in 0..params.poly_len {
            let idx = m * params.poly_len + i;
            // res_slc[idx] = barrett_coeff_u64(params, res_slc[idx], m);
            unsafe {
                let p = res_slc.as_mut_ptr().add(idx);
                *p = barrett_coeff_u64(params, *p, m);
            }
        }
    }
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
    precomp: &Precomp<'a>,
    b_values: &[u64],
    num_rlwe_outputs: usize,
    pack_pub_params_row_1s: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<PolyMatrixNTT<'a>> {
    assert_eq!(prep_rlwe_cts.len(), num_rlwe_outputs);
    assert_eq!(prep_rlwe_cts[0].len(), params.poly_len);
    assert_eq!(b_values.len(), num_rlwe_outputs * params.poly_len);

    let mut res = Vec::new();
    for i in 0..num_rlwe_outputs {
        let (precomp_res, precomp_vals, precomp_tables) = &precomp[i];

        let packed = pack_using_precomp_vals(
            &params,
            params.poly_len_log2,
            &pack_pub_params_row_1s,
            &b_values[i * params.poly_len..(i + 1) * params.poly_len],
            &precomp_res,
            &precomp_vals,
            &precomp_tables,
            &y_constants,
        );

        res.push(packed);
    }

    res
}

fn rotation_poly<'a>(params: &'a Params, amount: usize) -> PolyMatrixNTT<'a> {
    let mut res = PolyMatrixRaw::zero(params, 1, 1);
    res.data[amount] = 1;
    res.ntt()
}
pub fn pack_using_single_with_offset<'a>(
    params: &'a Params,
    pub_params: &[PolyMatrixNTT<'a>],
    cts: &[PolyMatrixNTT<'a>],
    offset: usize,
) -> PolyMatrixNTT<'a> {
    let mut res = PolyMatrixNTT::zero(params, 2, 1);
    for i in 0..cts.len() {
        let packed_single = pack_single_lwe(params, pub_params, &cts[i]);
        let rotation = rotation_poly(params, offset + i);
        let rotated = scalar_multiply_alloc(&rotation, &packed_single);
        add_into(&mut res, &rotated);
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
        let mut table_candidate = vec![0usize; params.poly_len];

        // for 2048 balls and 2^28 bins, we will have a collision ~1% of the time
        // so, we redo if necessary
        loop {
            let t = (1 << i) + 1;

            let poly = PolyMatrixRaw::random(&params, 1, 1);
            let poly_ntt = poly.ntt();

            let poly_auto = automorph_alloc(&poly, t);
            let poly_auto_ntt = poly_auto.ntt();

            let pol_orig = (&poly_ntt.get_poly(0, 0)[..params.poly_len]).to_vec();
            let pol_auto = (&poly_auto_ntt.get_poly(0, 0)[..params.poly_len]).to_vec();

            let mut must_redo = false;

            for i in 0..params.poly_len {
                let mut total = 0;
                let mut found = None;
                for j in 0..params.poly_len {
                    if pol_orig[i] == pol_auto[j] {
                        total += 1;
                        found = Some(j);
                    }
                }
                table_candidate[found.unwrap()] = i;
                if total != 1 {
                    must_redo = true;
                    break;
                }
            }

            if !must_redo {
                break;
            }
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
        out[i] += poly[table[i]];
    }
}

pub fn apply_automorph_ntt<'a>(
    params: &'a Params,
    tables: &[Vec<usize>],
    mat: &PolyMatrixNTT<'a>,
    res: &mut PolyMatrixNTT<'a>,
    t: usize,
) {
    // run apply_automorph_ntt on each poly in the matrix
    // let mut res = PolyMatrixNTT::zero(params, mat.rows, mat.cols);
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
    // res
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
        assert_eq!(packed_raw.get_poly(0, 0)[0], 47649720264253743u64);

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

        let mut fake_pack_pub_params = pack_pub_params.clone();
        // zero out all of the second rows
        for i in 0..pack_pub_params.len() {
            for col in 0..pack_pub_params[i].cols {
                fake_pack_pub_params[i].get_poly_mut(1, col).fill(0);
            }
        }

        let mut pack_pub_params_row_1s = pack_pub_params.clone();
        for i in 0..pack_pub_params.len() {
            pack_pub_params_row_1s[i] =
                pack_pub_params[i].submatrix(1, 0, 1, pack_pub_params[i].cols);
            pack_pub_params_row_1s[i] = condense_matrix(&params, &pack_pub_params_row_1s[i]);
        }

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
        let (precomp_res, precomp_vals, precomp_tables) = precompute_pack(
            &params,
            params.poly_len_log2,
            &v_ct,
            &fake_pack_pub_params,
            &y_constants,
        );
        println!(
            "Precomputing for packing took {} us",
            now.elapsed().as_micros()
        );

        println!("t_exp_left: {}", params.t_exp_left);

        let now = Instant::now();
        let packed = pack_using_precomp_vals(
            &params,
            params.poly_len_log2,
            &pack_pub_params_row_1s,
            &b_values,
            &precomp_res,
            &precomp_vals,
            &precomp_tables,
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
    fn test_single_packing() {
        let params = params_for_scenario(1 << 30, 1);
        let mut client = Client::init(&params);
        client.generate_secret_keys();

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

        // generate 1 ciphertext
        let sentinel_val = 99;
        let mut v_ct = Vec::new();
        for _i in 0..1 {
            let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
            let val = sentinel_val % params.pt_modulus;
            let scale_k = params.modulus / params.pt_modulus;
            let mod_inv = invert_uint_mod(params.poly_len as u64, params.modulus).unwrap();
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );
            v_ct.push(ct);
        }
        let now = Instant::now();
        let packed = pack_single_lwe(&params, &pack_pub_params, &v_ct[0]);
        println!("Packing took {} us", now.elapsed().as_micros());

        let packed_raw = packed.raw();
        println!("packed_0: {:?}", &packed_raw.get_poly(0, 0)[..10]);

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed);
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            rescaled.data[i] = rescale(dec_raw.data[i], params.modulus, params.pt_modulus);
        }

        println!("rescaled: {:?}", &rescaled.as_slice()[..50]);
        assert_eq!(rescaled.as_slice()[0], sentinel_val);
    }

    #[test]
    fn test_automorph_tables() {
        let params = params_for_scenario(1 << 30, 1);

        let now = Instant::now();
        let tables = generate_automorph_tables_brute_force(&params);
        println!("Generating tables took {} us", now.elapsed().as_micros());

        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);

        for t in ([3, 9, 17, 33, 65, 129, 257, 513, 1025, 2049])
            .into_iter()
            .rev()
        {
            let poly = PolyMatrixRaw::random_rng(&params, 1, 1, &mut rng);
            let poly_ntt = poly.ntt();

            let poly_auto = automorph_alloc(&poly, t);
            let poly_auto_ntt = poly_auto.ntt();
            let mut poly_auto_ntt_using_tables = PolyMatrixNTT::zero(&params, 1, 1);
            apply_automorph_ntt(
                &params,
                &tables,
                &poly_ntt,
                &mut poly_auto_ntt_using_tables,
                t,
            );

            println!("poly_ntt: {:?}", &poly_ntt.as_slice()[..30]);
            println!(
                "poly_auto_ntt_using_tables: {:?}",
                &poly_auto_ntt_using_tables.as_slice()[..30]
            );

            assert_eq!(
                &poly_auto_ntt.as_slice(),
                &poly_auto_ntt_using_tables.as_slice(),
                "t: {}",
                t
            );
        }
    }
}
