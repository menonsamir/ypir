use wasm_bindgen_test::wasm_bindgen_test;
use ypir_client::*;

#[wasm_bindgen_test]
fn rust_test() {
    let _ = generate_query_to_check_item_inclusion("test");
}
