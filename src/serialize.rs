use spiral_rs::{aligned_memory::AlignedMemory64, params::*, poly::*};

use crate::client::{YPIRQuery, YPIRSimpleQuery};

pub trait ToBytes {
    fn to_bytes(&self) -> Vec<u8>;
}

pub trait AsBytes {
    fn as_bytes(&self) -> &[u8];
}

pub trait FromBytes {
    fn from_bytes(data: &[u8]) -> Self;
}

impl ToBytes for &[u32] {
    fn to_bytes(&self) -> Vec<u8> {
        // fast
        unsafe {
            let ptr = self.as_ptr() as *const u8;
            std::slice::from_raw_parts(ptr, self.len() * 4).to_vec()
        }
    }
}

impl AsBytes for &[u32] {
    fn as_bytes(&self) -> &[u8] {
        // fast
        unsafe {
            let ptr = self.as_ptr() as *const u8;
            std::slice::from_raw_parts(ptr, self.len() * 4)
        }
    }
}

impl FromBytes for Vec<u32> {
    fn from_bytes(data: &[u8]) -> Self {
        // fast
        unsafe {
            let mut out = Vec::with_capacity(data.len() / 4);
            let u8_mut = std::slice::from_raw_parts(data.as_ptr(), data.len());
            out.set_len(data.len() / 4);
            let ptr = out.as_mut_ptr() as *mut u8;
            std::ptr::copy_nonoverlapping(u8_mut.as_ptr(), ptr, data.len());
            out
        }
    }
}

impl ToBytes for &[u64] {
    fn to_bytes(&self) -> Vec<u8> {
        // fast
        unsafe {
            let ptr = self.as_ptr() as *const u8;
            std::slice::from_raw_parts(ptr, self.len() * 8).to_vec()
        }
    }
}

impl AsBytes for &[u64] {
    fn as_bytes(&self) -> &[u8] {
        // fast
        unsafe {
            let ptr = self.as_ptr() as *const u8;
            std::slice::from_raw_parts(ptr, self.len() * 8)
        }
    }
}

impl FromBytes for AlignedMemory64 {
    fn from_bytes(data: &[u8]) -> Self {
        // fast
        unsafe {
            let mut out = AlignedMemory64::new(data.len() / 8);
            let u8_mut = std::slice::from_raw_parts(data.as_ptr(), data.len());
            let ptr = out.as_mut_ptr() as *mut u8;
            std::ptr::copy_nonoverlapping(u8_mut.as_ptr(), ptr, data.len());
            out
        }
    }
}

impl ToBytes for YPIRQuery {
    fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(self.0.as_slice().as_bytes());
        out.extend_from_slice(self.1.as_slice().as_bytes());
        out.extend_from_slice(self.2.as_slice().as_bytes());
        out
    }
}

impl ToBytes for YPIRSimpleQuery {
    fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(self.0.as_slice().as_bytes());
        out.extend_from_slice(self.1.as_slice().as_bytes());
        out
    }
}

pub fn pack_vec_pm(
    params: &Params,
    rows: usize,
    cols: usize,
    v_cts: &[PolyMatrixNTT],
) -> AlignedMemory64 {
    assert_eq!(v_cts[0].rows, rows);
    assert_eq!(v_cts[0].cols, cols);
    assert_eq!(params.crt_count, 2);
    let mut aligned_out = AlignedMemory64::new(v_cts.len() * rows * cols * params.poly_len);
    let mut iter = aligned_out
        .as_mut_slice()
        .chunks_exact_mut(rows * cols * params.poly_len);
    for ct in v_cts {
        let out = iter.next().unwrap();
        for row in 0..rows {
            for col in 0..cols {
                let out_offs = (row * cols + col) * params.poly_len;
                let inp_offs = (row * cols + col) * 2 * params.poly_len;
                for z in 0..params.poly_len {
                    out[out_offs + z] =
                        ct.data[inp_offs + z] | (ct.data[inp_offs + z + params.poly_len] << 32);
                }
            }
        }
    }
    aligned_out
}

pub fn unpack_vec_pm<'a>(
    params: &'a Params,
    rows: usize,
    cols: usize,
    data: &[u64],
) -> Vec<PolyMatrixNTT<'a>> {
    assert_eq!(params.crt_count, 2);
    let mut v_cts = Vec::with_capacity(data.len() / (rows * cols * params.poly_len));
    let mut iter = data.chunks_exact(rows * cols * params.poly_len);
    for _ in 0..v_cts.capacity() {
        let in_data = iter.next().unwrap();
        let mut ct = PolyMatrixNTT::zero(params, rows, cols);
        for row in 0..rows {
            for col in 0..cols {
                // this is (on purpose) not the inverse of pack_vec_pm
                let in_offs = (row * cols + col) * params.poly_len;
                let out_offs = (row * cols + col) * 2 * params.poly_len;
                for z in 0..params.poly_len {
                    ct.data[out_offs + z] = in_data[in_offs + z];
                }
            }
        }
        v_cts.push(ct);
    }
    v_cts
}

#[cfg(test)]
mod test {
    use crate::{packing::uncondense_matrix, params::params_for_scenario};

    use super::*;

    #[test]
    fn test_pack_unpack_vec_pm() {
        let params = params_for_scenario(1 << 10, 1);
        let rows = 5;
        let cols = 3;
        let len = 7;
        let mut v_cts = Vec::new();
        for _ in 0..len {
            v_cts.push(PolyMatrixRaw::random(&params, rows, cols).ntt());
        }
        let data = pack_vec_pm(&params, rows, cols, &v_cts);
        let v_cts2_weird = unpack_vec_pm(&params, rows, cols, data.as_slice());
        let v_cts2 = v_cts2_weird
            .into_iter()
            .map(|ct| uncondense_matrix(&params, &ct))
            .collect::<Vec<_>>();
        for (ct1, ct2) in v_cts.iter().zip(v_cts2.iter()) {
            assert_eq!(ct1.raw().as_slice(), ct2.raw().as_slice());
        }
    }
}
