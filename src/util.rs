use spiral_rs::{arith::*, params::*, poly::*, util};

pub fn negacyclic_perm(a: &[u64], shift: usize, modulus: u64) -> Vec<u64> {
    let n = a.len();
    let mut out = vec![0u64; n];

    for i in 0..shift + 1 {
        out[i] = a[shift - i];
    }

    for i in shift + 1..n {
        out[i] = modulus - (a[n - (i - shift)] % modulus);
        if out[i] == modulus {
            out[i] = 0;
        }
    }

    out
}

pub fn negacyclic_matrix(a: &[u64], modulus: u64) -> Vec<u64> {
    let n = a.len();
    let mut out = vec![0u64; n * n];

    for i in 0..n {
        let perm = negacyclic_perm(a, i, modulus);
        for j in 0..n {
            out[j * n + i] = perm[j];
        }
    }

    out
}

pub fn get_negacylic<'a>(poly: &PolyMatrixRaw<'a>) -> PolyMatrixRaw<'a> {
    let mut out = poly.clone();

    (&mut out.as_mut_slice()[1..]).reverse();

    for z in 1..out.data.len() {
        out.data[z] = out.params.modulus - out.data[z];
    }

    out
}

pub fn reduce_copy(params: &Params, out: &mut [u64], inp: &[u64]) {
    for n in 0..params.crt_count {
        for z in 0..params.poly_len {
            out[n * params.poly_len + z] = barrett_coeff_u64(params, inp[z], n);
        }
    }
}

pub fn add_into_no_reduce(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT) {
    assert!(res.rows == a.rows);
    assert!(res.cols == a.cols);

    for i in 0..res.rows {
        for j in 0..res.cols {
            let res_poly = res.get_poly_mut(i, j);
            let pol2 = a.get_poly(i, j);
            for (res_poly, pol2) in res_poly.iter_mut().zip(pol2.iter()) {
                *res_poly += *pol2;
            }
        }
    }
}

pub fn add_into_at_no_reduce(
    res: &mut PolyMatrixNTT,
    a: &PolyMatrixNTT,
    t_row: usize,
    t_col: usize,
) {
    let params = res.params;
    for i in 0..a.rows {
        for j in 0..a.cols {
            let res_poly = res.get_poly_mut(t_row + i, t_col + j);
            let pol2 = a.get_poly(i, j);
            for z in 0..params.crt_count * params.poly_len {
                res_poly[z] += pol2[z];
            }
        }
    }
}

pub fn modular_reduce_poly<'a>(a: &mut PolyMatrixNTT<'a>) {
    let params = a.params;
    for i in 0..a.rows {
        for j in 0..a.cols {
            let pol = a.get_poly_mut(i, j);
            for z in 0..params.crt_count * params.poly_len {
                pol[z] = barrett_coeff_u64(params, pol[z], z / params.poly_len);
            }
        }
    }
}

pub fn concat_horizontal(v_a: &[Vec<u64>], a_rows: usize, a_cols: usize) -> Vec<u64> {
    let mut out = vec![0u64; a_rows * a_cols * v_a.len()];

    for i in 0..a_rows {
        for j in 0..a_cols {
            for k in 0..v_a.len() {
                let idx = i * a_cols + j;
                let out_idx = i * a_cols * v_a.len() + k * a_cols + j;
                out[out_idx] = v_a[k][idx];
            }
        }
    }

    out
}

#[cfg(target_feature = "avx2")]
pub fn is_avx() -> bool {
    true
}

#[cfg(not(target_feature = "avx2"))]
pub fn is_avx() -> bool {
    false
}

pub fn test_params() -> Params {
    let params_str = r#"{
        "n": 1,
        "nu_1": 0,
        "nu_2": 0,
        "p": 256,
        "q2_bits": 22,
        "t_gsw": 3,
        "t_conv": 2,
        "t_exp_left": 2,
        "t_exp_right": 2,
        "instances": 1,
        "db_item_size": 0,
        "version": 2
    }"#;
    let params = util::params_from_json(&params_str);
    params
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_negacyclic_perm() {
        let a = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let modulus = 11;
        let shift = 3;

        let result = negacyclic_perm(&a, shift, modulus);
        assert_eq!(result, vec![4, 3, 2, 1, 11 - 8, 11 - 7, 11 - 6, 11 - 5]);
    }

    #[test]
    fn test_crt() {
        let params = test_params();

        let a = (1u64 << 45) + 77;

        let a0 = a % params.moduli[0];
        let a1 = a % params.moduli[1];

        let a_goal = params.crt_compose_2(a0, a1) % params.modulus;
        assert_eq!(a_goal, a);

        let a0_mod = (a0 * 5123) % params.moduli[0];
        let a1_mod = (a1 * 5123) % params.moduli[1];

        let a_modified = params.crt_compose_2(a0_mod, a1_mod) % params.modulus;
        assert_eq!(a_modified, (a * 5123) % params.modulus);
    }
}
