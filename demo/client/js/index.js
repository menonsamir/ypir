import init, { decodeResponse, decodeResponseToInclusionQuery, generateQuery, generateQueryToCheckItemInclusion } from "../pkg/ypir_client.js";
export { decodeResponse, decodeResponseToInclusionQuery, generateQuery, generateQueryToCheckItemInclusion, init };

export function postRequest(url, data) {
  data = new Uint8Array(data);
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
