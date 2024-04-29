"""
# HintlessPIR query + response size calculation
```
Upload = (k + l) (n log q) + db_cols (log Q)
Download = 2k(db_rows + n) log q + db_rows (log Q)

k     = 2
l     = 6
n     = 4096
log q = 90
log Q = 32
```

## DB sizes
```
 1 GB: 32768 x 32768
 2 GB: 32768 x 65536
 4 GB: 32768 x 131072
 8 GB: 32768 x 262144
16 GB: 65536 x 262144
32 GB: 65536 x 524288
```
"""

k = 2
l = 6
n = 4096
log_q = 90
log_Q = 32


def upload_bits(db_rows, db_cols):
    return (k + l) * n * log_q + db_cols * log_Q


def download_bits(db_rows, db_cols):
    return 2 * k * (db_rows + n) * log_q + db_rows * log_Q


db_sizes = {
    "1GB": (32768, 32768),
    "2GB": (32768, 65536),
    "4GB": (32768, 131072),
    "8GB": (32768, 262144),
    "16GB": (65536, 262144),
    "32GB": (65536, 524288),
}

for db_size, (db_rows, db_cols) in db_sizes.items():
    print(f"{db_size} DB ({db_rows} x {db_cols}):")
    print(f"Upload: {upload_bits(db_rows, db_cols)//8} bytes")
    print(f"Download: {download_bits(db_rows, db_cols)//8} bytes")
    print()
