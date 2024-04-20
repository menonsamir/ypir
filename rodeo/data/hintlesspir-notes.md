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

## Query
for 16384x16384:
= ((2 + 6) * (4096 * 90) + 2^14 * 32)/8 
= 434176
= 424 KB

for 32768x32768:
= ((2 + 6) * (4096 * 90) + 2^15 * 32)/8 
= 499712
= 488 KB

for 32768x262144:
= ((2 + 6) * (4096 * 90) + 2^18 * 32)/8 
= 1417216
= 1384 KB

for 65536x524288:
= ((2 + 6) * (4096 * 90) + 2^19 * 32)/8 
= 2465792
= 2408 KB


## Response
for 16384x16384:
= ((2 * 2) * (2^14 + 4096) * 90 + 2^14 * 32)/8 
= 987136
= 964 KB

for 32768x32768:
= ((2 * 2) * (2^15 + 4096) * 90 + 2^15 * 32)/8 
= 1789952
= 1748 KB

for 32768x262144:
= ((2 * 2) * (2^15 + 4096) * 90 + 2^15 * 32)/8 
= 1789952
= 1748 KB

for 65536x524288:
= ((2 * 2) * (2^16 + 4096) * 90 + 2^16 * 32)/8 
= 3395584
= 3316 KB

## Runtime
for 32768x32768:
1425 ms

for 32768x262144:
3506 ms

for 65536x524288:
?

```
taskset -c 0 sh -c "CC=clang++-14 bazel run -c opt --copt='-march=native' --cxxopt='-std=c++17' hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1 --benchmark_context=db_dims='1024x1024' --num_rows=65536 --num_cols=524288"
```

```
ubuntu@ip-172-31-4-28:~/hintless_pir$ taskset -c 0 sh -c "CC=clang++-14 bazel run -c opt --copt='-march=native' --cxxopt='-std=c++17' hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1 --num_rows=65536 --num_cols=524288"
INFO: Analyzed target //hintless_simplepir:hintless_simplepir_benchmarks (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //hintless_simplepir:hintless_simplepir_benchmarks up-to-date:
  bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks
INFO: Elapsed time: 3.979s, Critical Path: 3.68s
INFO: 3 processes: 1 internal, 2 linux-sandbox.
INFO: Build completed successfully, 3 total actions
INFO: Running command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' '--num_rows=65536' '--num_cols=524288'
exec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-04-18T22:54:36+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (32 X 3500.04 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x16)
  L1 Instruction 32 KiB (x16)
  L2 Unified 1280 KiB (x16)
  L3 Unified 55296 KiB (x1)
Load Average: 0.40, 0.22, 0.09
Preprocess: 5732511ms
SimplePIR: 10835ms
NTTlessPIR: 1949ms
SimplePIR: 10833ms
NTTlessPIR: 1938ms
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64 1.2785e+10 ns   1.2761e+10 ns            1
```

for 32768x32768:
```
LinPIR preprocessing: 16096ms
```