// import("../pkg/ypir_client.js").catch(console.error);
// export * from '../pkg/ypir_client.js';
import init, { decodeResponse, generateQuery } from "../pkg/ypir_client.js";
export { decodeResponse, generateQuery, init };

// export function postRequest(url, data) {
//   // data is Uint8Array, returns a Uint8Array
//   return fetch(url, {
//     method: 'POST',
//     body: data,
//     headers: {
//       'Content-Type': 'application/octet-stream',
//     },
//   }).then(response => {
//     if (!response.ok) {
//       throw new Error('Network response was not ok');
//     }
//     return response.arrayBuffer();
//   }).then(buffer => new Uint8Array(buffer));
// }

export function postRequest(url, data) {
  data = new Uint8Array(data);
  console.log("sending", data)
  return fetch(url, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/octet-stream'
    },
    body: data
  })
  .then(response => {
    if (!response.ok) {
      return Promise.reject(new Error(`HTTP error! Status: ${response.status}`));
    }
    return response.arrayBuffer();
  })
  .then(arrayBuffer => {
    return new Uint8Array(arrayBuffer);
  });
}
