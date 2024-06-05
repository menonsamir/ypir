use js_sys::wasm_bindgen;
use spiral_rs::{arith::rescale, client::*, params::*, poly::*};
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen_test::wasm_bindgen_test;
use ypir::{
    bits::u64s_to_contiguous_bytes, client::*, modulus_switch::ModulusSwitch, params::*,
    util::qt_hash,
};

#[wasm_bindgen]
extern "C" {
    // Use `js_namespace` here to bind `console.log(..)` instead of just
    // `log(..)`
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}
macro_rules! console_log {
    ($($t:tt)*) => (log(&format_args!($($t)*).to_string()))
}

// wasm_bindgen_test_configure!(run_in_browser);

fn custom_decode_response_simplepir_yclient(
    params: &Params,
    y_client: &YClient,
    response_data: &[u8],
) -> Vec<u64> {
    let db_cols = params.instances * params.poly_len;
    let num_rlwe_outputs = db_cols / params.poly_len;

    assert_eq!(response_data.len() % num_rlwe_outputs, 0);
    let response_vecs = response_data
        .chunks_exact(response_data.len() / num_rlwe_outputs)
        .map(|chunk| chunk.to_vec())
        .collect::<Vec<_>>();

    // rescale
    let rlwe_q_prime_1 = params.get_q_prime_1();
    let rlwe_q_prime_2 = params.get_q_prime_2();
    let mut response = Vec::new();
    for ct_bytes in response_vecs.iter() {
        let ct = PolyMatrixRaw::recover(&params, rlwe_q_prime_1, rlwe_q_prime_2, ct_bytes);
        response.push(ct);
    }

    // debug!("decrypting outer cts...");
    let outer_ct = response
        .iter()
        .flat_map(|ct| {
            custom_decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len)
                .as_slice()
                .to_vec()
        })
        .collect::<Vec<_>>();
    assert_eq!(outer_ct.len(), num_rlwe_outputs * params.poly_len);
    // debug!("outer_ct: {:?}", &outer_ct[..]);
    // let outer_ct_t_u8 = u64s_to_contiguous_bytes(&outer_ct, pt_bits);
    outer_ct
}

fn custom_decrypt_ct_reg_measured<'a>(
    client: &Client<'a>,
    params: &'a Params,
    ct: &PolyMatrixNTT<'a>,
    coeffs_to_measure: usize,
) -> PolyMatrixRaw<'a> {
    // console_log!("dec inp: {} x {}", ct.rows, ct.cols);
    // console_log!("dec inp len: {:?}", &ct.as_slice().len());
    // console_log!("dec inp: {:?}", &ct.raw().ntt().as_slice());
    // console_log!("dec inp: {:?}", &ct.raw().as_slice()[..16]);
    // console_log!(
    //     "self.sk_reg_full: {:?}",
    //     &client.sk_reg_full.as_slice()[..16]
    // );
    // console_log!(
    //     "self.sk_reg_full hash: {:?}",
    //     qt_hash(client.sk_reg_full.as_slice())
    // );
    // console_log!("ct.as_slice() hash: {:?}", qt_hash(ct.as_slice()));
    // console_log!("ct.as_slice() hash: {:?}", qt_hash(ct.as_slice()));
    // assert_eq!(pol.raw().as_slice(), pol.raw().ntt().raw().as_slice());
    let dec_result = client.decrypt_matrix_reg(ct).raw();
    // println!("dec_result: {:?}", &dec_result.as_slice()[..16]);
    // console_log!("dec_result: {:?}", &dec_result.as_slice()[..16]);

    let mut dec_rescaled = PolyMatrixRaw::zero(&params, dec_result.rows, dec_result.cols);
    for z in 0..dec_rescaled.data.len() {
        dec_rescaled.data[z] = rescale(dec_result.data[z], params.modulus, params.pt_modulus);
    }

    // measure noise width
    // let s_2 = measure_noise_width_squared(params, client, ct, &dec_rescaled, coeffs_to_measure);
    // debug!("log2(measured noise): {}", s_2.log2());

    dec_rescaled
}

// This runs a unit test in native Rust, so it can only use Rust APIs.
#[wasm_bindgen_test]
fn rust_test() {
    assert_eq!(1, 1);
    use ypir::data::{RESP1, SEED1};

    let client = YPIRClient::from_db_sz(1 << 14, 16384 * 8, true);
    let decoded = client.decode_response_simplepir(SEED1, RESP1);
    assert_eq!(decoded[0], 17);
    assert_eq!(decoded[1], 0);
    assert_eq!(decoded[2], 0);
    // console_log!("raw:  {:?}", &decoded[..32]);
    // let bytes = u64s_to_contiguous_bytes(&decoded, client.params().pt_modulus_bits());
    // console_log!("as u8: {:?}", &bytes[..32]);

    // let params = params_for_scenario_simplepir(1 << 14, 16384 * 8);
    // let mut client = Client::init(&params);
    // client.generate_secret_keys_from_seed(SEED1);
    // let y_client = YClient::from_seed(&mut client, &params, SEED1);
    // let decoded = custom_decode_response_simplepir_yclient(&params, &y_client, RESP1);
    // let decoded_bytes = u64s_to_contiguous_bytes(&decoded, params.pt_modulus_bits());
    // assert_eq!(decoded_bytes[0], 17);
    // assert_eq!(decoded_bytes[1], 0);
    // assert_eq!(decoded_bytes[2], 0);
}

// // This runs a unit test in the browser, so it can use browser APIs.
// #[wasm_bindgen_test]
// fn web_test() {
//     assert_eq!(1, 1);
// }

// #[wasm_bindgen_test]
// fn decode_test() {
//     use ypir::data::{RESP1, SEED1};

//     let client = YPIRClient::from_db_sz(1 << 14, 16384 * 8, true);
//     let decoded = client.decode_response_simplepir(SEED1, RESP1);
//     println!("raw:  {:?}", &decoded[..32]);
//     let bytes = u64s_to_contiguous_bytes(&decoded, client.params().pt_modulus_bits());
//     println!("as u8: {:?}", &bytes[..32]);
// }

// // This runs a unit test in the browser, and in addition it supports asynchronous Future APIs.
// #[wasm_bindgen_test(async)]
// fn async_test() -> impl Future<Item = (), Error = JsValue> {
//     // Creates a JavaScript Promise which will asynchronously resolve with the value 42.
//     let promise = js_sys::Promise::resolve(&JsValue::from(42));

//     // Converts that Promise into a Future.
//     // The unit test will wait for the Future to resolve.
//     JsFuture::from(promise).map(|x| {
//         assert_eq!(x, 42);
//     })
// }
