use std::convert::TryInto;

use serde::{Deserialize, Serialize};
use spiral_rs::{arith::rescale, client::Client, params::*, poly::*};
use wasm_bindgen::prelude::*;
use web_sys::console;
use ypir::{
    bits::{read_bits, u64s_to_contiguous_bytes, write_bits},
    client::*,
    modulus_switch::ModulusSwitch,
    params::*,
    serialize::ToBytes,
};

// #[wasm_bindgen]
// pub struct State {
//     client: YPIRClient,
// }

fn get_db_conf() -> (u64, u64, bool) {
    // TODO: don't hardcode this
    let num_items = 1u64 << 14;
    let bits_per_item = 16384u64 * 8;
    let is_simplepir = true;
    (num_items, bits_per_item, is_simplepir)
}

fn test_stuff() {
    let mut buffer = [0u8; 4];

    // Test write_bits
    write_bits(&mut buffer, 0b11010101, 1, 6);
    assert_eq!(buffer, [0b00101010, 0b00000000, 0b00000000, 0b00000000]);

    // Test read_bits
    let value = read_bits(&buffer, 1, 6);
    assert_eq!(value, 0b010101);

    // Additional tests
    let mut buffer2 = [0u8; 4];
    write_bits(&mut buffer2, 0b11111111, 0, 8);
    assert_eq!(buffer2, [0b11111111, 0b00000000, 0b00000000, 0b00000000]);
    let value2 = read_bits(&buffer2, 0, 8);
    assert_eq!(value2, 0b11111111);

    let mut buffer3 = [0u8; 4];
    write_bits(&mut buffer3, 0b10101010, 4, 4);
    assert_eq!(buffer3, [0b10100000, 0b00000000, 0b00000000, 0b00000000]);
    let value3 = read_bits(&buffer3, 4, 4);
    assert_eq!(value3, 0b1010);

    let num = 42;
    let bits_per = 14;
    let total_sy_bytes = ((num * bits_per) as f64 / 8.0).ceil() as usize;
    let vals = (0..num)
        .map(|i| (17 * i as u64 ^ 3 >> 1) % (1 << bits_per))
        .collect::<Vec<_>>();
    let mut buffer4 = vec![0u8; total_sy_bytes];
    let mut bit_offs = 0;
    for i in 0..num {
        write_bits(&mut buffer4, vals[i], bit_offs, bits_per);
        bit_offs += bits_per;
    }

    for i in 0..num {
        let val = read_bits(&buffer4, i * bits_per, bits_per);
        assert_eq!(val, vals[i]);
    }
}

// #[wasm_bindgen]
// pub fn client_init() -> State {
//     let (num_items, bits_per_item, is_simplepir) = get_db_conf();
//     let client = YPIRClient::from_db_sz(num_items, bits_per_item, is_simplepir);
//     let state = State { client };
//     state
// }

// #[wasm_bindgen]
// pub struct QueryAndSeed(pub Vec<u8>, pub Vec<u8>);

#[derive(Serialize, Deserialize)]
pub struct QueryAndSeed {
    pub query: Vec<u8>,
    pub seed: Vec<u8>,
}

#[wasm_bindgen(js_name = generateQuery)]
pub fn generate_query(desired_index: usize) -> Result<JsValue, JsValue> {
    test_stuff();

    let (num_items, bits_per_item, is_simplepir) = get_db_conf();
    let client = YPIRClient::from_db_sz(num_items, bits_per_item, is_simplepir);
    let (query, seed) = client.generate_query_simplepir(desired_index as usize);
    web_sys::console::log_1(&JsValue::from_str(&format!("seed {:?}", seed)));
    let query_bytes = query.to_bytes();
    assert_eq!(query_bytes.len(), 671744);
    Ok(serde_wasm_bindgen::to_value(&QueryAndSeed {
        query: query_bytes,
        seed: seed.to_vec(),
    })?)
}

#[wasm_bindgen(js_name = decodeResponse)]
pub fn decode_response(response: &[u8], seed: &[u8]) -> Vec<u8> {
    web_sys::console::log_1(&JsValue::from_str(&format!("seed {:?}", seed)));
    web_sys::console::log_1(&JsValue::from_str(&format!(
        "response.len() {:?}",
        response.len()
    )));
    web_sys::console::log_1(&JsValue::from_str(&format!(
        "response {:?}...",
        &response[..32]
    )));
    let (num_items, bits_per_item, is_simplepir) = get_db_conf();
    let client = YPIRClient::from_db_sz(num_items, bits_per_item, is_simplepir);

    // let mut client = Client::init(&self.params);
    // client.generate_secret_keys_from_seed(client_seed);
    // let y_client = YClient::from_seed(&mut client, &self.params, client_seed);
    // y_client

    let decoded = client.decode_response_simplepir(seed.try_into().unwrap(), &response);
    web_sys::console::log_1(&JsValue::from_str(&format!(
        "decoded: {:?}...",
        &decoded[..64]
    )));
    decoded
}

#[wasm_bindgen(start)]
pub fn main_js() -> Result<(), JsValue> {
    // This provides better error messages in debug mode.
    // It's disabled in release mode so it doesn't bloat up the file size.
    #[cfg(debug_assertions)]
    console_error_panic_hook::set_once();

    // Your code goes here!
    console::log_1(&JsValue::from_str("Hello world!"));

    Ok(())
}
