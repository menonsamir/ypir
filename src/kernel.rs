use std::arch::x86_64::*;
use std::ops::Range;

use log::debug;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{arith::*, params::*};

use super::server::{ToM512, ToU64};

pub fn fast_dot_product_avx512(
    params: &Params,
    a: &[u64],
    a_elems: usize,
    b_t: &[u16], // transposed
    b_rows: usize,
    b_cols: usize,
) -> AlignedMemory64 {
    assert_eq!(a_elems, b_rows);

    debug!("Multiplying {}x{} by {}x{}", 1, a_elems, b_rows, b_cols);

    let simd_width = 8;

    let chunk_size = 8192.min(a_elems / simd_width);
    let num_chunks = (a_elems / simd_width) / chunk_size;
    debug!("k_chunk_size: {}, k_num_chunks: {}", chunk_size, num_chunks);

    let j_chunk_size = 1;
    let j_num_chunks = b_cols / j_chunk_size;
    debug!(
        "j_chunk_size: {}, j_num_chunks: {}",
        j_chunk_size, j_num_chunks
    );

    let mut result = AlignedMemory64::new(b_cols);
    let res_mut_slc = result.as_mut_slice();

    unsafe {
        let a_ptr = a.as_ptr();
        let b_ptr = b_t.as_ptr();
        // let res_ptr = res_mut_slc.as_mut_ptr();

        for k_outer in 0..num_chunks {
            for j_outer in 0..j_num_chunks {
                for j_inner in 0..j_chunk_size {
                    let j = j_outer * j_chunk_size + j_inner;

                    let mut total_sum_lo = _mm512_setzero_si512();
                    let mut total_sum_hi = _mm512_setzero_si512();

                    for k_inner in 0..chunk_size {
                        let k = simd_width * (k_outer * chunk_size + k_inner);

                        let a_idx = k;
                        let b_idx = j * b_rows + k;

                        let a_val_simd = _mm512_load_si512(a_ptr.add(a_idx) as *const _);
                        let a_val_lo = a_val_simd;
                        let a_val_hi = _mm512_srli_epi64(a_val_simd, 32);

                        let b_val_simd = b_ptr.add(b_idx).to_m512();

                        total_sum_lo =
                            _mm512_add_epi64(total_sum_lo, _mm512_mul_epu32(a_val_lo, b_val_simd));
                        total_sum_hi =
                            _mm512_add_epi64(total_sum_hi, _mm512_mul_epu32(a_val_hi, b_val_simd));
                    }

                    let mut values_lo = [0u64; 8];
                    let mut values_hi = [0u64; 8];
                    _mm512_store_si512((&mut values_lo).as_mut_ptr() as *mut _, total_sum_lo);
                    _mm512_store_si512((&mut values_hi).as_mut_ptr() as *mut _, total_sum_hi);

                    let res_lo = values_lo[0]
                        + values_lo[1]
                        + values_lo[2]
                        + values_lo[3]
                        + values_lo[4]
                        + values_lo[5]
                        + values_lo[6]
                        + values_lo[7];

                    let res_hi = values_hi[0]
                        + values_hi[1]
                        + values_hi[2]
                        + values_hi[3]
                        + values_hi[4]
                        + values_hi[5]
                        + values_hi[6]
                        + values_hi[7];

                    let (lo, hi) = (
                        barrett_coeff_u64(params, res_lo as u64, 0),
                        barrett_coeff_u64(params, res_hi as u64, 1),
                    );

                    let res = params.crt_compose_2(lo, hi);
                    res_mut_slc[j] = barrett_u64(params, res_mut_slc[j] + res);
                }
            }
        }
    }

    result
}

pub fn fast_batched_dot_product_avx512<const K: usize, T: Copy>(
    params: &Params,
    c: &mut [u64],
    a: &[u64],
    a_elems: usize,
    b_t: &[T], // transposed
    b_rows: usize,
    b_cols: usize,
) where
    *const T: ToM512,
{
    assert_eq!(a_elems, b_rows);

    // debug!("Multiplying {}x{} by {}x{}", K, a_elems, b_rows, b_cols);

    let simd_width = 8;

    let chunk_size = (8192 / K.next_power_of_two()).min(a_elems / simd_width);
    let num_chunks = (a_elems / simd_width) / chunk_size;
    // debug!("k_chunk_size: {}, k_num_chunks: {}", chunk_size, num_chunks);

    let j_chunk_size = 1;
    let j_num_chunks = b_cols / j_chunk_size;
    // debug!(
    //     "j_chunk_size: {}, j_num_chunks: {}",
    //     j_chunk_size, j_num_chunks
    // );

    // let mut result = AlignedMemory64::new(b_cols);
    // let res_mut_slc = result.as_mut_slice();
    let res_mut_slc = c;

    unsafe {
        let mut a_slcs: [&[u64]; K] = [&[]; K];
        for (slc_mut, chunk) in a_slcs.iter_mut().zip(a.chunks_exact(a.len() / K)) {
            *slc_mut = chunk;
        }

        // let a_ptr = a.as_ptr();
        let b_ptr = b_t.as_ptr();

        for k_outer in 0..num_chunks {
            for j_outer in 0..j_num_chunks {
                for j_inner in 0..j_chunk_size {
                    let j = j_outer * j_chunk_size + j_inner;

                    let mut total_sum_lo = [_mm512_setzero_si512(); K];
                    let mut total_sum_hi = [_mm512_setzero_si512(); K];
                    let mut tmp = [_mm512_setzero_si512(); K];

                    for k_inner in 0..chunk_size {
                        let k = simd_width * (k_outer * chunk_size + k_inner);

                        let a_idx = k;
                        let b_idx = j * b_rows + k;
                        let b_val_simd = b_ptr.add(b_idx).to_m512();

                        for batch in 0..K {
                            tmp[batch] =
                                _mm512_load_si512(a_slcs[batch].as_ptr().add(a_idx) as *const _);
                        }

                        for batch in 0..K {
                            let a_val_lo = tmp[batch];
                            let a_val_hi = _mm512_srli_epi64(tmp[batch], 32);

                            total_sum_lo[batch] = _mm512_add_epi64(
                                total_sum_lo[batch],
                                _mm512_mul_epu32(a_val_lo, b_val_simd),
                            );
                            total_sum_hi[batch] = _mm512_add_epi64(
                                total_sum_hi[batch],
                                _mm512_mul_epu32(a_val_hi, b_val_simd),
                            );
                        }
                    }

                    let res_mut_slcs = res_mut_slc.chunks_exact_mut(res_mut_slc.len() / K);
                    for (batch, res_mut_slc) in (0..K).zip(res_mut_slcs) {
                        let mut values_lo = [0u64; 8];
                        let mut values_hi = [0u64; 8];
                        _mm512_store_si512(
                            (&mut values_lo).as_mut_ptr() as *mut _,
                            total_sum_lo[batch],
                        );
                        _mm512_store_si512(
                            (&mut values_hi).as_mut_ptr() as *mut _,
                            total_sum_hi[batch],
                        );

                        let res_lo = values_lo[0]
                            + values_lo[1]
                            + values_lo[2]
                            + values_lo[3]
                            + values_lo[4]
                            + values_lo[5]
                            + values_lo[6]
                            + values_lo[7];

                        let res_hi = values_hi[0]
                            + values_hi[1]
                            + values_hi[2]
                            + values_hi[3]
                            + values_hi[4]
                            + values_hi[5]
                            + values_hi[6]
                            + values_hi[7];

                        // res_mut_slc[j] = res_lo + res_hi;

                        let (lo, hi) = (
                            barrett_coeff_u64(params, res_lo as u64, 0),
                            barrett_coeff_u64(params, res_hi as u64, 1),
                        );

                        // res_mut_slc[j] = lo | (hi << 32);

                        let res = params.crt_compose_2(lo, hi);
                        res_mut_slc[j] = barrett_u64(params, res_mut_slc[j] + res);
                    }
                }
            }
        }
    }
}

pub fn lwe_fast_batched_dot_product_avx512<const K: usize, T: Copy>(
    _params: &Params,
    c: &mut [u64],
    a: &[u64],
    a_elems: usize,
    b_t: &[T], // transposed
    b_rows: usize,
    b_cols: usize,
) where
    *const T: ToM512,
{
    assert_eq!(a_elems, b_rows);

    debug!("Multiplying {}x{} by {}x{}", K, a_elems, b_rows, b_cols);

    let simd_width = 8;

    let chunk_size = (8192 / K.next_power_of_two()).min(a_elems / simd_width);
    let num_chunks = (a_elems / simd_width) / chunk_size;
    assert_eq!(num_chunks, 1);
    // debug!("k_chunk_size: {}, k_num_chunks: {}", chunk_size, num_chunks);

    let j_chunk_size = 1;
    let j_num_chunks = b_cols / j_chunk_size;
    // debug!(
    //     "j_chunk_size: {}, j_num_chunks: {}",
    //     j_chunk_size, j_num_chunks
    // );

    // let mut result = AlignedMemory64::new(b_cols);
    // let res_mut_slc = result.as_mut_slice();
    let res_mut_slc = c;

    unsafe {
        let mut a_slcs: [&[u64]; K] = [&[]; K];
        for (slc_mut, chunk) in a_slcs.iter_mut().zip(a.chunks_exact(a.len() / K)) {
            *slc_mut = chunk;
        }

        // let a_ptr = a.as_ptr();
        let b_ptr = b_t.as_ptr();

        for k_outer in 0..num_chunks {
            for j_outer in 0..j_num_chunks {
                for j_inner in 0..j_chunk_size {
                    let j = j_outer * j_chunk_size + j_inner;

                    let mut total_sum_lo = [_mm512_setzero_si512(); K];

                    for k_inner in 0..chunk_size {
                        let k = simd_width * (k_outer * chunk_size + k_inner);

                        let a_idx = k;
                        let b_idx = j * b_rows + k;
                        let b_val_simd = b_ptr.add(b_idx).to_m512();

                        for batch in 0..K {
                            let a_val_lo =
                                _mm512_load_si512(a_slcs[batch].as_ptr().add(a_idx) as *const _);

                            total_sum_lo[batch] = _mm512_add_epi64(
                                total_sum_lo[batch],
                                _mm512_mul_epu32(a_val_lo, b_val_simd),
                            );
                        }
                    }

                    let res_mut_slcs = res_mut_slc.chunks_exact_mut(res_mut_slc.len() / K);
                    for (batch, res_mut_slc) in (0..K).zip(res_mut_slcs) {
                        let mut values_lo = [0u64; 8];
                        _mm512_store_si512(
                            (&mut values_lo).as_mut_ptr() as *mut _,
                            total_sum_lo[batch],
                        );

                        let res_lo = values_lo[0]
                            + values_lo[1]
                            + values_lo[2]
                            + values_lo[3]
                            + values_lo[4]
                            + values_lo[5]
                            + values_lo[6]
                            + values_lo[7];

                        res_mut_slc[j] = (res_mut_slc[j] + res_lo) % (1u64 << 32);
                    }
                }
            }
        }
    }
}

const MASK: u32 = 0x0000_00FF;
const BASIS: u32 = 8;
const BASIS2: u32 = 16;

pub unsafe fn lwe_u32_fast_batched_dot_product_avx512<const K: usize>(
    out: *mut u32,
    a: *const u32, // reversed a and b
    a_rows: usize,
    a_cols: usize,
    b: *const u32,
) {
    let mut index = 0;
    for i in 0..a_rows {
        let mut tmp = 0;
        let mut tmp2 = 0;
        let mut tmp3 = 0;
        let mut tmp4 = 0;
        let mut tmp5 = 0;
        let mut tmp6 = 0;
        let mut tmp7 = 0;
        let mut tmp8 = 0;

        let mut index2 = 0;
        for _j in 0..a_cols {
            let db = *a.add(index);
            let db2 = *a.add(index + 1 * a_cols);
            let db3 = *a.add(index + 2 * a_cols);
            let db4 = *a.add(index + 3 * a_cols);
            let db5 = *a.add(index + 4 * a_cols);
            let db6 = *a.add(index + 5 * a_cols);
            let db7 = *a.add(index + 6 * a_cols);
            let db8 = *a.add(index + 7 * a_cols);

            let mut val = db & MASK;
            let mut val2 = db2 & MASK;
            let mut val3 = db3 & MASK;
            let mut val4 = db4 & MASK;
            let mut val5 = db5 & MASK;
            let mut val6 = db6 & MASK;
            let mut val7 = db7 & MASK;
            let mut val8 = db8 & MASK;
            tmp += val * (*b.add(index2));
            tmp2 += val2 * (*b.add(index2));
            tmp3 += val3 * (*b.add(index2));
            tmp4 += val4 * (*b.add(index2));
            tmp5 += val5 * (*b.add(index2));
            tmp6 += val6 * (*b.add(index2));
            tmp7 += val7 * (*b.add(index2));
            tmp8 += val8 * (*b.add(index2));
            index2 += 1;

            val = (db >> BASIS) & MASK;
            val2 = (db2 >> BASIS) & MASK;
            val3 = (db3 >> BASIS) & MASK;
            val4 = (db4 >> BASIS) & MASK;
            val5 = (db5 >> BASIS) & MASK;
            val6 = (db6 >> BASIS) & MASK;
            val7 = (db7 >> BASIS) & MASK;
            val8 = (db8 >> BASIS) & MASK;
            tmp += val * (*b.add(index2));
            tmp2 += val2 * (*b.add(index2));
            tmp3 += val3 * (*b.add(index2));
            tmp4 += val4 * (*b.add(index2));
            tmp5 += val5 * (*b.add(index2));
            tmp6 += val6 * (*b.add(index2));
            tmp7 += val7 * (*b.add(index2));
            tmp8 += val8 * (*b.add(index2));
            index2 += 1;

            val = (db >> BASIS2) & MASK;
            val2 = (db2 >> BASIS2) & MASK;
            val3 = (db3 >> BASIS2) & MASK;
            val4 = (db4 >> BASIS2) & MASK;
            val5 = (db5 >> BASIS2) & MASK;
            val6 = (db6 >> BASIS2) & MASK;
            val7 = (db7 >> BASIS2) & MASK;
            val8 = (db8 >> BASIS2) & MASK;
            tmp += val * (*b.add(index2));
            tmp2 += val2 * (*b.add(index2));
            tmp3 += val3 * (*b.add(index2));
            tmp4 += val4 * (*b.add(index2));
            tmp5 += val5 * (*b.add(index2));
            tmp6 += val6 * (*b.add(index2));
            tmp7 += val7 * (*b.add(index2));
            tmp8 += val8 * (*b.add(index2));
            index2 += 1;
            index += 1;
        }
        *out.add(i) += tmp;
        *out.add(i + 1) += tmp2;
        *out.add(i + 2) += tmp3;
        *out.add(i + 3) += tmp4;
        *out.add(i + 4) += tmp5;
        *out.add(i + 5) += tmp6;
        *out.add(i + 6) += tmp7;
        *out.add(i + 7) += tmp8;
        index += a_cols * 7;
    }
}

// pub fn lwe_u32_fast_batched_dot_product_avx512<const K: usize>(
//     _params: &Params,
//     out: &mut [u32],
//     a: &[u32], // reversed a and b
//     a_rows: usize,
//     a_cols: usize,
//     b: &[u32],
// ) {
//     let mut index = 0;
//     for i in 0..a_rows {
//         let mut tmp = 0;
//         let mut tmp2 = 0;
//         let mut tmp3 = 0;
//         let mut tmp4 = 0;
//         let mut tmp5 = 0;
//         let mut tmp6 = 0;
//         let mut tmp7 = 0;
//         let mut tmp8 = 0;

//         let mut index2 = 0;
//         for j in 0..a_cols {
//             let db = a[index];
//             let db2 = a[index + 1 * a_cols];
//             let db3 = a[index + 2 * a_cols];
//             let db4 = a[index + 3 * a_cols];
//             let db5 = a[index + 4 * a_cols];
//             let db6 = a[index + 5 * a_cols];
//             let db7 = a[index + 6 * a_cols];
//             let db8 = a[index + 7 * a_cols];

//             let mut val = db & MASK;
//             let mut val2 = db2 & MASK;
//             let mut val3 = db3 & MASK;
//             let mut val4 = db4 & MASK;
//             let mut val5 = db5 & MASK;
//             let mut val6 = db6 & MASK;
//             let mut val7 = db7 & MASK;
//             let mut val8 = db8 & MASK;
//             tmp += val * b[index2];
//             tmp2 += val2 * b[index2];
//             tmp3 += val3 * b[index2];
//             tmp4 += val4 * b[index2];
//             tmp5 += val5 * b[index2];
//             tmp6 += val6 * b[index2];
//             tmp7 += val7 * b[index2];
//             tmp8 += val8 * b[index2];
//             index2 += 1;

//             val = (db >> BASIS) & MASK;
//             val2 = (db2 >> BASIS) & MASK;
//             val3 = (db3 >> BASIS) & MASK;
//             val4 = (db4 >> BASIS) & MASK;
//             val5 = (db5 >> BASIS) & MASK;
//             val6 = (db6 >> BASIS) & MASK;
//             val7 = (db7 >> BASIS) & MASK;
//             val8 = (db8 >> BASIS) & MASK;
//             tmp += val * b[index2];
//             tmp2 += val2 * b[index2];
//             tmp3 += val3 * b[index2];
//             tmp4 += val4 * b[index2];
//             tmp5 += val5 * b[index2];
//             tmp6 += val6 * b[index2];
//             tmp7 += val7 * b[index2];
//             tmp8 += val8 * b[index2];
//             index2 += 1;

//             val = (db >> BASIS2) & MASK;
//             val2 = (db2 >> BASIS2) & MASK;
//             val3 = (db3 >> BASIS2) & MASK;
//             val4 = (db4 >> BASIS2) & MASK;
//             val5 = (db5 >> BASIS2) & MASK;
//             val6 = (db6 >> BASIS2) & MASK;
//             val7 = (db7 >> BASIS2) & MASK;
//             val8 = (db8 >> BASIS2) & MASK;
//             tmp += val * b[index2];
//             tmp2 += val2 * b[index2];
//             tmp3 += val3 * b[index2];
//             tmp4 += val4 * b[index2];
//             tmp5 += val5 * b[index2];
//             tmp6 += val6 * b[index2];
//             tmp7 += val7 * b[index2];
//             tmp8 += val8 * b[index2];
//             index2 += 1;
//             index += 1;
//         }
//         out[i] += tmp;
//         out[i + 1] += tmp2;
//         out[i + 2] += tmp3;
//         out[i + 3] += tmp4;
//         out[i + 4] += tmp5;
//         out[i + 5] += tmp6;
//         out[i + 6] += tmp7;
//         out[i + 7] += tmp8;
//         index += a_cols * 7;
//     }
// }

// pub fn lwe_u32_fast_batched_dot_product_avx512<const K: usize>(
//     _params: &Params,
//     c: &mut [u32],
//     a: &[u32],
//     a_elems: usize,
//     b_t: &[u32], // transposed
//     b_rows: usize,
//     b_cols: usize,
// ) {
//     assert_eq!(a_elems, b_rows);
//     assert_eq!(a_elems % 4, 0);
//     assert_eq!(b_cols % 8, 0);

//     debug!("Multiplying {}x{} by {}x{}", K, a_elems, b_rows, b_cols);

//     for col in (0..b_cols).step_by(8) {
//         let (mut sum0, mut sum1, mut sum2, mut sum3, mut sum4, mut sum5, mut sum6, mut sum7) =
//             (0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
//         for row in 0..b_rows {
//             let db_0 = b_t[(col * b_rows) + row + 0 * b_rows];
//             let db_1 = b_t[(col * b_rows) + row + 1 * b_rows];
//             let db_2 = b_t[(col * b_rows) + row + 2 * b_rows];
//             let db_3 = b_t[(col * b_rows) + row + 3 * b_rows];
//             let db_4 = b_t[(col * b_rows) + row + 4 * b_rows];
//             let db_5 = b_t[(col * b_rows) + row + 5 * b_rows];
//             let db_6 = b_t[(col * b_rows) + row + 6 * b_rows];
//             let db_7 = b_t[(col * b_rows) + row + 7 * b_rows];

//             sum0 = sum0.wrapping_add(((db_0 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 0]));
//             sum1 = sum1.wrapping_add(((db_1 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 1]));
//             sum2 = sum2.wrapping_add(((db_2 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 2]));
//             sum3 = sum3.wrapping_add(((db_3 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 3]));
//             sum4 = sum4.wrapping_add(((db_4 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 4]));
//             sum5 = sum5.wrapping_add(((db_5 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 5]));
//             sum6 = sum6.wrapping_add(((db_6 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 6]));
//             sum7 = sum7.wrapping_add(((db_7 >> (0 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 7]));

//             sum0 = sum0.wrapping_add(((db_0 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 0]));
//             sum1 = sum1.wrapping_add(((db_1 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 1]));
//             sum2 = sum2.wrapping_add(((db_2 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 2]));
//             sum3 = sum3.wrapping_add(((db_3 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 3]));
//             sum4 = sum4.wrapping_add(((db_4 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 4]));
//             sum5 = sum5.wrapping_add(((db_5 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 5]));
//             sum6 = sum6.wrapping_add(((db_6 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 6]));
//             sum7 = sum7.wrapping_add(((db_7 >> (1 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 7]));

//             sum0 = sum0.wrapping_add(((db_0 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 0]));
//             sum1 = sum1.wrapping_add(((db_1 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 1]));
//             sum2 = sum2.wrapping_add(((db_2 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 2]));
//             sum3 = sum3.wrapping_add(((db_3 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 3]));
//             sum4 = sum4.wrapping_add(((db_4 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 4]));
//             sum5 = sum5.wrapping_add(((db_5 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 5]));
//             sum6 = sum6.wrapping_add(((db_6 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 6]));
//             sum7 = sum7.wrapping_add(((db_7 >> (2 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 7]));

//             sum0 = sum0.wrapping_add(((db_0 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 0]));
//             sum1 = sum1.wrapping_add(((db_1 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 1]));
//             sum2 = sum2.wrapping_add(((db_2 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 2]));
//             sum3 = sum3.wrapping_add(((db_3 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 3]));
//             sum4 = sum4.wrapping_add(((db_4 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 4]));
//             sum5 = sum5.wrapping_add(((db_5 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 5]));
//             sum6 = sum6.wrapping_add(((db_6 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 6]));
//             sum7 = sum7.wrapping_add(((db_7 >> (3 * 8)) & MASK_LO).wrapping_mul(a[4 * row + 7]));
//         }

//         c[col + 0] = c[col + 0].wrapping_add(sum0);
//         c[col + 1] = c[col + 1].wrapping_add(sum1);
//         c[col + 2] = c[col + 2].wrapping_add(sum2);
//         c[col + 3] = c[col + 3].wrapping_add(sum3);
//         c[col + 4] = c[col + 4].wrapping_add(sum4);
//         c[col + 5] = c[col + 5].wrapping_add(sum5);
//         c[col + 6] = c[col + 6].wrapping_add(sum6);
//         c[col + 7] = c[col + 7].wrapping_add(sum7);
//     }
// }

pub fn multiply_matrices_raw<T>(
    params: &Params,
    a: &[u64],
    a_rows: usize,
    a_cols: usize,
    b_t: &[T], // transposed
    b_rows: usize,
    _b_cols: usize,
    col_range: Range<usize>,
) -> AlignedMemory64
where
    T: ToU64 + Copy,
{
    assert_eq!(a_cols, b_rows);

    debug!(
        "Multiplying {}x{} by {}x[{:?}]",
        a_rows, a_cols, b_rows, col_range
    );

    let mut result = vec![0u128; a_rows * col_range.clone().len()];

    for i in 0..a_rows {
        for k in 0..a_cols {
            for j in col_range.clone() {
                let a_idx = i * a_cols + k;
                let b_idx = j * b_rows + k; // on purpose, since transposed
                let res_idx = i * col_range.len() + (j - col_range.start);

                unsafe {
                    let a_val = *a.get_unchecked(a_idx);
                    let b_val = (*b_t.get_unchecked(b_idx)).to_u64();

                    // result[res_idx] += a_val as u128 * b_val as u128;

                    // result[res_idx] = (result[res_idx] + (a_val * b_val)) % params.modulus;

                    // *result.get_unchecked_mut(res_idx) += a_val * b_val;
                    // }

                    // result[res_idx] += barrett_u64(params, a_val * b_val);
                    let prod = a_val as u128 * b_val as u128;
                    result[res_idx] += prod;
                }

                // if j % 128 == 0 {
                //     result[res_idx] = barrett_u64(params, result[res_idx]); // do every _ accumulations
                // }
            }
        }
    }

    // for i in 0..result.len() {
    //     result[i] %= params.modulus;
    // }

    // result

    let mut result_u64 = AlignedMemory64::new(a_rows * col_range.len());
    for i in 0..result.len() {
        result_u64[i] = barrett_reduction_u128(params, result[i]);
    }

    result_u64
}

pub fn as_u32_slc(a: &[u64]) -> &[u32] {
    unsafe { std::slice::from_raw_parts(a.as_ptr() as *const u32, a.len() * 2) }
}

pub fn as_u32_slc_mut(a: &mut [u64]) -> &mut [u32] {
    unsafe { std::slice::from_raw_parts_mut(a.as_mut_ptr() as *mut u32, a.len() * 2) }
}

pub fn multiply_matrices<T: ToU64 + Copy>(
    a: &[u32],
    a_rows: usize,
    a_cols: usize,
    b_t: &[T], // transposed
    b_rows: usize,
    b_cols: usize,
    is_b_transposd: bool,
) -> Vec<u32> {
    // performs wrapping arithmetic

    assert_eq!(a_cols, b_rows);

    // debug!("Multiplying {}x{} by {}x{}", a_rows, a_cols, b_rows, b_cols);

    let mut result = vec![0u32; a_rows * b_cols];
    for i in 0..a_rows {
        for j in 0..b_cols {
            for k in 0..a_cols {
                let a_idx = i * a_cols + k;
                let b_idx = if is_b_transposd {
                    j * b_rows + k // on purpose, since transposed
                } else {
                    k * b_cols + j
                };
                let res_idx = i * b_cols + j;

                unsafe {
                    let a_val = *a.get_unchecked(a_idx);
                    let b_val = (b_t.get_unchecked(b_idx)).to_u64() as u32;

                    result[res_idx] = result[res_idx].wrapping_add(a_val.wrapping_mul(b_val));
                }
            }
        }
    }

    result
}

#[cfg(test)]
mod test {
    use std::time::Instant;

    use spiral_rs::poly::*;

    use super::super::util::test_params;
    use super::*;

    #[test]
    fn test_fast_batched_dot_product_avx512() {
        let params = test_params();

        const A_ROWS: usize = 8;
        let a_cols = 2048;
        let b_rows = a_cols;
        let b_cols = 4096;

        let a = PolyMatrixRaw::random(&params, A_ROWS, a_cols);
        let mut b = AlignedMemory64::new(b_rows * b_cols);
        let mut c = AlignedMemory64::new(A_ROWS * b_cols);
        let trials = 128;
        let mut sum = 0u64;
        let mut sum_time = 0;
        for _ in 0..trials {
            for i in 0..b.len() {
                b[i] = fastrand::u64(..);
            }
            let b_u16_slc =
                unsafe { std::slice::from_raw_parts(b.as_ptr() as *const u16, b.len() * 4) };

            let now = Instant::now();
            // fast_dot_product_avx512(&params, a.as_slice(), rows, b_as_t_slice, rows, cols);
            // fast_dot_product_avx512_fastest(&params, a.as_slice(), rows, b_as_t_slice, rows, cols);
            fast_batched_dot_product_avx512::<A_ROWS, _>(
                &params,
                c.as_mut_slice(),
                a.as_slice(),
                a_cols,
                b_u16_slc,
                b_rows,
                b_cols,
            );
            sum_time += now.elapsed().as_micros();
            sum += c.as_slice()[fastrand::usize(..c.len())];
        }
        debug!(
            "fast_matmul_avx512 in {} us ({}: {} x {})",
            sum_time, trials, a_cols, b_cols
        );
        debug!("");
        debug!("{}", sum);
    }

    #[test]
    fn test_lwe_fast_batched_dot_product_avx512() {
        let params = test_params();

        const A_ROWS: usize = 1;
        let a_cols = 32768;
        let b_rows = a_cols;
        let b_cols = 32768;

        let a = PolyMatrixRaw::random(&params, A_ROWS, a_cols);
        let mut b = AlignedMemory64::new((b_rows + 16) * b_cols);
        let mut c = AlignedMemory64::new(A_ROWS * b_cols + 16);
        let trials = 1;
        let mut sum = 0u64;
        let mut sum_time = 0;
        for _ in 0..trials {
            for i in 0..b.len() {
                b[i] = fastrand::u64(..);
            }
            let b_u32_slc = as_u32_slc(b.as_slice());
            // let b_u8_slc =
            //     unsafe { std::slice::from_raw_parts(b.as_ptr() as *const u8, b.len() * 2) };

            let now = Instant::now();
            unsafe {
                lwe_u32_fast_batched_dot_product_avx512::<A_ROWS>(
                    as_u32_slc_mut(c.as_mut_slice()).as_mut_ptr(),
                    b_u32_slc.as_ptr(),
                    b_rows,
                    b_cols / 4,
                    as_u32_slc(a.as_slice()).as_ptr(),
                );
            }
            sum_time += now.elapsed().as_micros();
            sum += c.as_slice()[fastrand::usize(..c.len())];
        }
        debug!(
            "fast_matmul_avx512 in {} us ({}: {} x {})",
            sum_time, trials, a_cols, b_cols
        );
        debug!("");
        debug!("{}", sum);
    }
}
