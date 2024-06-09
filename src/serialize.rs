use std::{
    fs::File,
    io::{BufReader, Read, Seek, Write},
};

use spiral_rs::{aligned_memory::AlignedMemory64, params::*, poly::*, util::read_arbitrary_bits};

use crate::client::{YPIRQuery, YPIRSimpleQuery};

pub type Precomp<'a> = Vec<(PolyMatrixNTT<'a>, Vec<PolyMatrixNTT<'a>>, Vec<Vec<usize>>)>;

#[cfg(feature = "server")]
use crate::server::YServer;

#[cfg(not(feature = "server"))]
type YServer<'a, T> = T;

#[derive(Clone)]
pub struct OfflinePrecomputedValues<'a> {
    pub hint_0: Vec<u64>,
    pub hint_1: Vec<u64>,
    pub pseudorandom_query_1: Vec<PolyMatrixNTT<'a>>,
    pub y_constants: (Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
    pub smaller_server: Option<YServer<'a, u16>>,
    pub prepacked_lwe: Vec<Vec<PolyMatrixNTT<'a>>>,
    pub fake_pack_pub_params: Vec<PolyMatrixNTT<'a>>,
    pub precomp: Precomp<'a>,
}

pub trait ToBytes {
    fn to_bytes(&self) -> Vec<u8>;
}

pub trait AsBytes {
    fn as_bytes(&self) -> &[u8];
}

pub trait AsBytesMut {
    fn as_bytes_mut(&mut self) -> &mut [u8];
}

pub trait FromBytes {
    fn from_bytes(data: &[u8]) -> Self;
}

pub trait FromBytesParams<'a> {
    fn from_bytes(data: &[u8], params: &'a Params) -> Self;
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

impl AsBytesMut for &mut [u64] {
    fn as_bytes_mut(&mut self) -> &mut [u8] {
        // fast
        unsafe {
            let ptr = self.as_mut_ptr() as *mut u8;
            std::slice::from_raw_parts_mut(ptr, self.len() * 8)
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

impl FromBytes for Vec<u64> {
    fn from_bytes(data: &[u8]) -> Self {
        // fast
        unsafe {
            let mut out = Vec::with_capacity(data.len() / 8);
            let u8_mut = std::slice::from_raw_parts(data.as_ptr(), data.len());
            out.set_len(data.len() / 8);
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

impl ToBytes for Vec<usize> {
    fn to_bytes(&self) -> Vec<u8> {
        // fast
        unsafe {
            let ptr = self.as_ptr() as *const u8;
            std::slice::from_raw_parts(ptr, self.len() * 8).to_vec()
        }
    }
}

impl FromBytesParams<'_> for Vec<usize> {
    fn from_bytes(data: &[u8], _params: &Params) -> Vec<usize> {
        // fast
        unsafe {
            let mut out = Vec::with_capacity(data.len() / std::mem::size_of::<usize>());
            let u8_mut = std::slice::from_raw_parts(data.as_ptr(), data.len());
            out.set_len(data.len() / std::mem::size_of::<usize>());
            let ptr = out.as_mut_ptr() as *mut u8;
            std::ptr::copy_nonoverlapping(u8_mut.as_ptr(), ptr, data.len());
            out
        }
    }
}

impl<'a> ToBytes for PolyMatrixNTT<'a> {
    fn to_bytes(&self) -> Vec<u8> {
        // write rows, cols, and data
        let mut out = Vec::new();
        out.extend_from_slice(&self.rows.to_be_bytes());
        out.extend_from_slice(&self.cols.to_be_bytes());
        out.extend_from_slice(self.as_slice().as_bytes());
        out
    }
}

impl<'a> FromBytesParams<'a> for PolyMatrixNTT<'a> {
    fn from_bytes(data: &[u8], params: &'a Params) -> PolyMatrixNTT<'a> {
        let rows = u64::from_be_bytes(data[0..8].try_into().unwrap()) as usize;
        let cols = u64::from_be_bytes(data[8..16].try_into().unwrap()) as usize;
        let data = &data[16..];
        let mut out = PolyMatrixNTT::zero(params, rows, cols);
        out.as_mut_slice().as_bytes_mut().copy_from_slice(data);
        out
    }
}

impl<T: ToBytes> ToBytes for Vec<T> {
    fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&(self.len() as u64).to_be_bytes());
        for item in self {
            let item_bytes = item.to_bytes();
            out.extend_from_slice(&(item_bytes.len() as u64).to_be_bytes());
            out.extend_from_slice(&item_bytes);
        }
        out
    }
}

impl<'a, T: FromBytesParams<'a>> FromBytesParams<'a> for Vec<T> {
    fn from_bytes(data: &[u8], params: &'a Params) -> Vec<T> {
        let mut out = Vec::new();
        let mut data = data;
        let len = u64::from_be_bytes(data[0..8].try_into().unwrap()) as usize;
        data = &data[8..];
        for _ in 0..len {
            let item_len = u64::from_be_bytes(data[0..8].try_into().unwrap()) as usize;
            let item = T::from_bytes(&data[8..8 + item_len], params);
            out.push(item);
            data = &data[8 + item_len..];
        }
        out
    }
}

// length 2 tuple
impl<T1: ToBytes, T2: ToBytes> ToBytes for (T1, T2) {
    fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let item1 = self.0.to_bytes();
        let item2 = self.1.to_bytes();
        out.extend_from_slice(&(item1.len() as u64).to_be_bytes());
        out.extend_from_slice(&item1);
        out.extend_from_slice(&item2);
        out
    }
}

impl<'a, T1: FromBytesParams<'a>, T2: FromBytesParams<'a>> FromBytesParams<'a> for (T1, T2) {
    fn from_bytes(data: &[u8], params: &'a Params) -> (T1, T2) {
        let len1 = u64::from_be_bytes(data[0..8].try_into().unwrap()) as usize;
        let item1 = T1::from_bytes(&data[8..8 + len1], params);
        let item2 = T2::from_bytes(&data[8 + len1..], params);
        (item1, item2)
    }
}

// length 3 tuple
impl<T1: ToBytes, T2: ToBytes, T3: ToBytes> ToBytes for (T1, T2, T3) {
    fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let item1 = self.0.to_bytes();
        let item2 = self.1.to_bytes();
        let item3 = self.2.to_bytes();
        out.extend_from_slice(&(item1.len() as u64).to_be_bytes());
        out.extend_from_slice(&item1);
        out.extend_from_slice(&(item2.len() as u64).to_be_bytes());
        out.extend_from_slice(&item2);
        out.extend_from_slice(&item3);
        out
    }
}

impl<'a, T1: FromBytesParams<'a>, T2: FromBytesParams<'a>, T3: FromBytesParams<'a>>
    FromBytesParams<'a> for (T1, T2, T3)
{
    fn from_bytes(data: &[u8], params: &'a Params) -> (T1, T2, T3) {
        let len1 = u64::from_be_bytes(data[0..8].try_into().unwrap()) as usize;
        let item1 = T1::from_bytes(&data[8..8 + len1], params);
        let len2 = u64::from_be_bytes(data[8 + len1..16 + len1].try_into().unwrap()) as usize;
        let item2 = T2::from_bytes(&data[16 + len1..16 + len1 + len2], params);
        let item3 = T3::from_bytes(&data[16 + len1 + len2..], params);
        (item1, item2, item3)
    }
}

// pub hint_0: Vec<u64>,
// pub hint_1: Vec<u64>,
// pub pseudorandom_query_1: Vec<PolyMatrixNTT<'a>>,
// pub y_constants: (Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
// pub smaller_server: Option<YServer<'a, u16>>,
// pub prepacked_lwe: Vec<Vec<PolyMatrixNTT<'a>>>,
// pub fake_pack_pub_params: Vec<PolyMatrixNTT<'a>>,
// pub precomp: Precomp<'a>,

impl<'a> ToBytes for OfflinePrecomputedValues<'a> {
    fn to_bytes(&self) -> Vec<u8> {
        assert!(self.smaller_server.is_none());

        let mut out = Vec::new();
        out.extend_from_slice(&self.hint_0.as_slice().to_bytes());
        out.extend_from_slice(&self.hint_1.as_slice().to_bytes());
        out.extend_from_slice(&self.pseudorandom_query_1.to_bytes());
        out.extend_from_slice(&self.y_constants.0.to_bytes());
        out.extend_from_slice(&self.y_constants.1.to_bytes());
        out.extend_from_slice(&self.prepacked_lwe.to_bytes());
        out.extend_from_slice(&self.fake_pack_pub_params.to_bytes());
        out.extend_from_slice(&self.precomp.to_bytes());
        out
    }
}

impl<'a> FromBytesParams<'a> for OfflinePrecomputedValues<'a> {
    fn from_bytes(data: &[u8], params: &'a Params) -> OfflinePrecomputedValues<'a> {
        let mut data = data;
        let hint_0 = Vec::<u64>::from_bytes(data);
        data = &data[hint_0.len() * 8..];
        let hint_1 = Vec::<u64>::from_bytes(data);
        data = &data[hint_1.len() * 8..];
        let pseudorandom_query_1 = Vec::<PolyMatrixNTT>::from_bytes(data, params);
        data = &data[pseudorandom_query_1.to_bytes().len()..];
        let y_constants = (
            Vec::<PolyMatrixNTT>::from_bytes(data, params),
            Vec::<PolyMatrixNTT>::from_bytes(data, params),
        );
        data = &data[y_constants.0.to_bytes().len()..];
        let prepacked_lwe = Vec::<Vec<PolyMatrixNTT>>::from_bytes(data, params);
        data = &data[prepacked_lwe.to_bytes().len()..];
        let fake_pack_pub_params = Vec::<PolyMatrixNTT>::from_bytes(data, params);
        data = &data[fake_pack_pub_params.to_bytes().len()..];
        let precomp = Precomp::from_bytes(data, params);
        OfflinePrecomputedValues {
            hint_0,
            hint_1,
            pseudorandom_query_1,
            y_constants,
            smaller_server: None,
            prepacked_lwe,
            fake_pack_pub_params,
            precomp,
        }
    }
}

pub fn read_file_to_vec_u64(filename: &str) -> Vec<u64> {
    let mut file = File::open(filename).unwrap();
    let mut data = Vec::new();
    file.read_to_end(&mut data).unwrap();
    let mut out = Vec::with_capacity(data.len() / 8);
    let mut iter = data.chunks_exact(8);
    for _ in 0..iter.len() {
        out.push(u64::from_le_bytes(iter.next().unwrap().try_into().unwrap()));
    }
    out
}

pub fn write_vec_u64_to_file(filename: &str, data: &[u64]) {
    let mut file = File::create(filename).unwrap();
    for &x in data {
        file.write_all(&x.to_le_bytes()).unwrap();
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

pub struct FilePtIter<R: Read + Seek> {
    file: R,
    pt_bits: usize,
    bytes_per_row: usize,
    db_cols: usize,
    row_idx: usize,
    col_idx: usize,
    buf_pos: usize,
    buf: Vec<u8>,
    buf_vals: Vec<u16>,
}

impl<R: Read + Seek> FilePtIter<R> {
    pub fn new(file: R, bytes_per_row: usize, db_cols: usize, pt_bits: usize) -> Self {
        let max_filled_col = bytes_per_row / pt_bits;
        assert!(max_filled_col <= db_cols);

        Self {
            file,
            pt_bits,
            bytes_per_row,
            db_cols,
            col_idx: 0,
            row_idx: 0,
            buf_pos: 8,
            buf: vec![0; pt_bits * 16],
            buf_vals: vec![0; pt_bits * 8],
        }
    }
}

impl FilePtIter<BufReader<File>> {
    pub fn from_file(filename: &str, bytes_per_row: usize, db_cols: usize, pt_bits: usize) -> Self {
        println!("bytes_per_row: {}, pt_bits: {}", bytes_per_row, pt_bits);
        Self::new(
            BufReader::new(File::open(filename).unwrap()),
            bytes_per_row,
            db_cols,
            pt_bits,
        )
    }
}

impl<R: Read + Seek> Iterator for FilePtIter<R> {
    type Item = u16;

    fn next(&mut self) -> Option<Self::Item> {
        // reads file, pt_bits at a time

        // max_filled_col pt-bits sized words contain data in each row (rest are zeros)
        let max_filled_col = self.bytes_per_row * 8 / self.pt_bits;
        assert!(max_filled_col <= self.db_cols);
        if self.col_idx >= self.db_cols {
            self.col_idx = 0;
            self.row_idx += 1;
            let seeked_to = self
                .file
                .seek(std::io::SeekFrom::Start(
                    self.row_idx as u64 * self.bytes_per_row as u64,
                ))
                .unwrap();
            assert_eq!(seeked_to, self.row_idx as u64 * self.bytes_per_row as u64);
            self.buf_vals.fill(0);
            self.buf_pos = 8;
        } else if self.col_idx >= max_filled_col {
            self.col_idx += 1;
            return Some(0);
        }

        // if at end of buffer, read pt_bits * 8 BITS (pt_bits bytes) of data
        if self.buf_pos == 8 {
            self.buf_pos = 0;
            self.buf.fill(0);

            let read = self.file.read_exact(&mut self.buf[..self.pt_bits]);
            if read.is_err() {
                // TODO: behavior at end-of-file
                // For now, just produce zeros infinitely now
                self.buf_pos = 8;
                return Some(0);
            }

            // now, populate buf_vals with the (up to) 8 pt_bits-sized words
            self.buf_vals.fill(0);
            let mut bit_offs = 0;
            for i in 0..8 {
                let val = read_arbitrary_bits(&self.buf, bit_offs, self.pt_bits);
                self.buf_vals[i] = val as u16;
                bit_offs += self.pt_bits;
            }
        }

        // return the next value in buf_vals
        let val = self.buf_vals[self.buf_pos];
        self.buf_pos += 1;
        self.col_idx += 1;
        Some(val)
    }
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use crate::{
        bits::u64s_to_contiguous_bytes,
        params::{params_for_scenario, params_for_scenario_simplepir, DbRowsCols, PtModulusBits},
        server::{ToU64, YServer},
    };

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

    #[test]
    fn test_pt_iter() {
        let num_items = 1 << 14;
        let item_size_bytes = 16384;

        let params = params_for_scenario_simplepir(num_items, item_size_bytes as u64 * 8);
        let inp_data = (0..params.db_rows())
            .flat_map(|i| {
                let mut out = vec![0u8; item_size_bytes];
                out[0] = i as u8;
                (&mut out[1..5]).copy_from_slice(&(i as u32).to_be_bytes());
                out
            })
            .collect::<Vec<_>>();
        let cursor = Cursor::new(inp_data);
        let pt_iter = FilePtIter::new(
            cursor,
            item_size_bytes,
            params.db_cols_simplepir(),
            params.pt_modulus_bits(),
        );
        let y_server = YServer::<u16>::new(&params, pt_iter, true, false, true);

        for i in 0..params.db_rows() {
            let row = y_server
                .get_row(i)
                .iter()
                .map(|x| x.to_u64())
                .collect::<Vec<_>>();
            let ci_bytes = u64s_to_contiguous_bytes(&row, params.pt_modulus_bits());
            assert_eq!(ci_bytes[0], i as u8);
            assert_eq!(
                i as u32,
                u32::from_be_bytes(ci_bytes[1..5].try_into().unwrap())
            );
        }
    }
}
