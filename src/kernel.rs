use std::arch::x86_64::*;

use spiral_rs::{arith::*, params::*};

use super::server::ToM512;

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

#[cfg(test)]
mod test {
    use std::time::Instant;

    use log::debug;
    use spiral_rs::aligned_memory::AlignedMemory64;
    use spiral_rs::poly::*;

    use super::super::util::test_params;
    use super::*;
    use test_log::test;

    #[test]
    fn test_fast_batched_dot_product_avx512() {
        let params = test_params();

        const A_ROWS: usize = 1;
        let a_cols = 131072;
        let b_rows = a_cols;
        let b_cols = 32768;

        let a = PolyMatrixRaw::random(&params, A_ROWS, a_cols);
        let mut b = AlignedMemory64::new(b_rows * b_cols);
        let mut c = AlignedMemory64::new(A_ROWS * b_cols);
        let trials = 1;
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
}
