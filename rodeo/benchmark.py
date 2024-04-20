import json
import math
import subprocess
import sys

# ex run: python rodeo/benchmark.py ypir rodeo/workload.json rodeo/output.json


def run_benchmark(
    scheme: str,
    num_items: int,
    item_size_bits: int,
    num_clients: int,
    trials: int,
    is_simplepir: bool,
):
    out_json_filename = "report.json"
    sp_flag = "--is-simplepir" if is_simplepir else ""

    cmd = f"./target/release-with-debug/run {num_items} {item_size_bits} {num_clients} {trials} {out_json_filename} {sp_flag}"
    print(cmd)

    # run commmand and get output
    subprocess.run(cmd, shell=True, cwd=".")

    # read output json file
    json_str = None
    with open(out_json_filename, "r") as fh:
        json_str = fh.read()

    # print json
    data = json.loads(json_str)

    # delete json file
    subprocess.check_output(f"rm {out_json_filename}", shell=True)

    return data


def run_benchmarks(scheme: str, trials: int, workload_file: str, output_json_file: str):
    workload_json = None
    with open(workload_file, "r") as fh:
        workload_json = json.load(fh)

    template_json_file = "rodeo/rodeo.json"
    merged_json = None
    with open(template_json_file, "r") as fh:
        merged_json = json.load(fh)

    results = []
    for scenario in workload_json["workloads"]:
        num_items = int(scenario["db"]["numItems"])
        item_size_bits = int(scenario["db"]["itemSizeBits"])
        num_clients = 1
        if "clients" in scenario and "numClients" in scenario["clients"]:
            num_clients = int(scenario["clients"]["numClients"])
        is_simplepir = scheme == "ypir-simplepir"
        measurement = run_benchmark(
            scheme, num_items, item_size_bits, num_clients, trials, is_simplepir
        )
        result = {"scenario": scenario, "measurement": measurement}
        results.append(result)

        # merge and write
        merged_json["results"] = results
        with open(output_json_file, "w") as fh:
            json.dump(merged_json, fh, indent=2)


if __name__ == "__main__":
    scheme = sys.argv[1]
    workload_file = sys.argv[2]
    output_json_file = sys.argv[3]
    trials = 5
    if len(sys.argv) > 4:
        trials = int(sys.argv[4])
    print(f"Running benchmark for {scheme} with {trials} trials", file=sys.stderr)
    run_benchmarks(scheme, trials, workload_file, output_json_file)
