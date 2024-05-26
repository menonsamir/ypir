# Notes

For commit `4be2ae8` in `google/hintless_pir`.

`time bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1 --num_rows 32768 --num_cols 262144`

`bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1  --benchmark_out_format=json --num_rows 32768 --num_cols 65536 --benchmark_out=2GB.json && bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1  --benchmark_out_format=json --num_rows 32768 --num_cols 131072 --benchmark_out=4GB.json && bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1  --benchmark_out_format=json --num_rows 65536 --num_cols 262144 --benchmark_out=16GB.json`

## DB sizes
```
 1 GB: 32768 x 32768
 2 GB: 32768 x 65536
 4 GB: 32768 x 131072
 8 GB: 32768 x 262144
16 GB: 65536 x 262144
32 GB: 65536 x 524288

1GB DB (32768 x 32768):
Upload: 499712 bytes
Download: 1789952 bytes

2GB DB (32768 x 65536):
Upload: 630784 bytes
Download: 1789952 bytes

4GB DB (32768 x 131072):
Upload: 892928 bytes
Download: 1789952 bytes

8GB DB (32768 x 262144):
Upload: 1417216 bytes
Download: 1789952 bytes

16GB DB (65536 x 262144):
Upload: 1417216 bytes
Download: 3395584 bytes

32GB DB (65536 x 524288):
Upload: 2465792 bytes
Download: 3395584 bytes
```

## Runs

### 1 GB
```
ubuntu@ip-172-31-15-66:~/hintless_pir$ taskset -c 0 bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repe
titions=1 --benchmark_min_time=0.1 --num_rows 32768 --num_cols 32768
INFO: Analyzed target //hintless_simplepir:hintless_simplepir_benchmarks (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //hintless_simplepir:hintless_simplepir_benchmarks up-to-date:
  bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks
INFO: Elapsed time: 0.152s, Critical Path: 0.00s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
INFO: Running command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' --num_rows 32768 --num_cols 32768
exec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-05-24T03:19:58+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (16 X 3500.53 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1280 KiB (x8)
  L3 Unified 55296 KiB (x1)
Load Average: 7.37, 2.86, 1.04
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64  742863894 ns    742844681 ns            1
```

## 2 GB

```
INFO: Running command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' '--benchmark_out_format=json' --num_rows 32768 --num_cols 65536 '--benchmark_out=2GB.json'
exec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-05-26T00:44:30+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (16 X 3519.87 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1280 KiB (x8)
  L3 Unified 55296 KiB (x1)
Load Average: 9.98, 3.02, 1.04
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64  860082865 ns    860056388 ns            1
```

## 4 GB

```
32mINFO: 0mRunning command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' --num_rows 32768 --num_cols 131072
0mexec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-05-24T03:38:59+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (16 X 3500.07 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1280 KiB (x8)
  L3 Unified 55296 KiB (x1)
Load Average: 1.00, 1.39, 1.34
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64 1121242523 ns   1121213617 ns            1
```


## 8 GB
```
ubuntu@ip-172-31-15-66:~/hintless_pir$ time bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1 --num_rows 32768 --num_cols 262144
INFO: Analyzed target //hintless_simplepir:hintless_simplepir_benchmarks (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //hintless_simplepir:hintless_simplepir_benchmarks up-to-date:
  bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks
INFO: Elapsed time: 0.087s, Critical Path: 0.00s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
INFO: Running command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' --num_rows 32768 --num_cols 262144
exec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-05-23T21:55:15+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (16 X 3500.07 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1280 KiB (x8)
  L3 Unified 55296 KiB (x1)
Load Average: 0.07, 0.27, 0.51
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64 1619558811 ns   1619518229 ns            1

real    27m31.718s
user    27m27.183s
sys     0m5.445s
```
## 16 GB

```
32mINFO: 0mRunning command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' --num_rows 65536 --num_cols 262144
0mexec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-05-24T03:53:03+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (16 X 3501.39 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1280 KiB (x8)
  L3 Unified 55296 KiB (x1)
Load Average: 1.00, 1.01, 1.11
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64 2954728842 ns   2954720931 ns            1
```

## 32 GB
```
ubuntu@ip-172-31-15-66:~/hintless_pir$ time bazel run -c opt hintless_simplepir/hintless_simplepir_benchmarks --  --benchmark_filter=BM_HintlessPirRlwe64 --benchmark_repetitions=1 --benchmark_min_time=0.1 --num_rows 65536 --num_cols 524288
INFO: Analyzed target //hintless_simplepir:hintless_simplepir_benchmarks (0 packages loaded, 0 targets configured).
INFO: Found 1 target...
Target //hintless_simplepir:hintless_simplepir_benchmarks up-to-date:
  bazel-bin/hintless_simplepir/hintless_simplepir_benchmarks
INFO: Elapsed time: 0.069s, Critical Path: 0.00s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
INFO: Running command line: external/bazel_tools/tools/test/test-setup.sh hintless_simplepir/hintless_simplepir_benchmarks '--benchmark_filter=BM_HintlessPirRlwe64' '--benchmark_repetitions=1' '--benchmark_min_time=0.1' --num_rows 65536 --num_cols 524288
exec ${PAGER:-/usr/bin/less} "$0" || exit 1
Executing tests from //hintless_simplepir:hintless_simplepir_benchmarks
-----------------------------------------------------------------------------
2024-05-23T22:46:30+00:00
Running /home/ubuntu/.cache/bazel/_bazel_ubuntu/2aa86797b13019953d76ea99599532c8/execroot/_main/bazel-out/k8-opt/bin/hintless_simplepir/hintless_simplepir_benchmarks.runfiles/_main/hintless_simplepir/hintless_simplepir_benchmarks
Run on (16 X 3483.54 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x8)
  L1 Instruction 32 KiB (x8)
  L2 Unified 1280 KiB (x8)
  L3 Unified 55296 KiB (x1)
Load Average: 0.00, 0.00, 0.17
---------------------------------------------------------------
Benchmark                     Time             CPU   Iterations
---------------------------------------------------------------
BM_HintlessPirRlwe64 4997175694 ns   4997114429 ns            1

real    107m44.998s
user    107m31.511s
sys     0m17.099s
```