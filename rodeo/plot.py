import argparse
import csv
import json
import math
import sys
from calendar import c

import matplotlib.pyplot as plt
import numpy as np
import tikzplotlib

scheme_styles = {
    "ypir": {"color": "red"},
    "simplepir": {"color": "green"},
    "doublepir": {"color": "blue"},
    "simplepir*": {"color": "forestgreen"},
    "doublepir*": {"color": "royalblue"},
    "zpir": {"color": "yellow"},
    "tiptoe": {"color": "purple"},
    "hintlesspir": {"color": "orange"},
    "hintlesspir*": {"color": "teal"},
}

replace = {
    "semithick, red": "ypir",
    "semithick, green": "simplepir",
    "semithick, blue": "doublepir",
    "semithick, forestgreen": "simplepirstar",
    "semithick, royalblue": "doublepirstar",
    "semithick, yellow": "zpir",
    "semithick, purple": "tiptoe",
    "semithick, orange": "hintlesspir",
    "semithick, teal": "hintlesspirplus",
    "\\addplot [ypir]": "\\addlegendentry{\ypir}\n\\addplot [ypir]",
    "\\addplot [zpir]": "\\addlegendentry{\zpir}\n\\addplot [zpir]",
    "\\addplot [simplepir]": "\\addlegendentry{\simplepir}\n\\addplot [simplepir]",
    "\\addplot [doublepir]": "\\addlegendentry{\doublepir}\n\\addplot [doublepir]",
    "\\addplot [simplepirstar]": "\\addlegendentry{\simplepirstar}\n\\addplot [simplepirstar]",
    "\\addplot [doublepirstar]": "\\addlegendentry{\doublepirstar}\n\\addplot [doublepirstar]",
    "\\addplot [tiptoe]": "\\addlegendentry{\\tiptoe}\n\\addplot [tiptoe]",
    "\\addplot [hintlesspir]": "\\addlegendentry{\\hintlesspir}\n\\addplot [hintlesspir]",
    "\\addplot [hintlesspirplus]": "\\addlegendentry{\\hintlesspirplus}\n\\addplot [hintlesspirplus]",
}

# can change this size to 1000 if desired
BASE_SZ = 1024

# these are really constants (don't change)
KB_SZ = BASE_SZ
MB_SZ = BASE_SZ * KB_SZ
GB_SZ = BASE_SZ * MB_SZ
TB_SZ = BASE_SZ * GB_SZ
TRUE_GB_SZ = 1024 * 1024 * 1024
MS_PER_S = 1000


def gather_1_bit_retrieval_data(data_files_json: list[str]):
    # each file has something like 'rodeo.json'

    # make sure any files containing "prep" are last
    data_files_json = sorted(data_files_json, key=lambda x: 1 if "prep" in x else 0)

    # load data, each file has a list of results for a scheme
    data = []
    for data_file in data_files_json:
        with open(data_file, "r") as fh:
            data.append(json.load(fh))

    # get all scenarios
    scenarios = []
    for datum in data:
        for result in datum["results"]:
            scenario = result["scenario"]
            scenarios.append(scenario)

    # dedup scenarios by combination of numItems and itemSizeBits
    scenarios = list(
        {
            (scenario["db"]["numItems"], scenario["db"]["itemSizeBits"]): scenario
            for scenario in scenarios
        }.values()
    )

    # get result of each scheme for each scenario
    scheme_results = {}
    for datum in data:
        scheme = datum["scheme"]["variant"]
        if scheme not in scheme_results:
            scheme_results[scheme] = {}
        for result in datum["results"]:
            scenario = result["scenario"]
            if "clients" in scenario and scenario["clients"]["numClients"] != 1:
                continue
            scenario_key = (scenario["db"]["numItems"], scenario["db"]["itemSizeBits"])
            measurement = result["measurement"]
            if scenario_key in scheme_results[scheme]:
                cur_data = scheme_results[scheme][scenario_key]
                # if we're missing the offline time, use the one from this one
                if (
                    cur_data["offline"]["serverTimeMs"] == 0
                    and measurement["offline"]["serverTimeMs"] > 0
                ):
                    cur_data["offline"]["serverTimeMs"] = measurement["offline"][
                        "serverTimeMs"
                    ]
                    # don't use any other data
                    continue
            scheme_results[scheme][scenario_key] = measurement

    # get all scenarios
    all_scenarios = None
    for scheme, results in scheme_results.items():
        if all_scenarios is None:
            all_scenarios = set(results.keys())
        else:
            all_scenarios = all_scenarios.union(results.keys())
    all_scenarios = sorted(list(all_scenarios))

    # remove 64 GB
    if (549755813888, 1) in all_scenarios:
        all_scenarios.remove((549755813888, 1))

    # calculate the percent std. dev. for each scheme
    percent_devs = []
    for scheme, results in scheme_results.items():
        for scenario, measurement in results.items():
            if "stdDevServerTimeMs" not in measurement["online"]:
                print("No std dev for", scheme, scenario)
                continue
            std_dev = measurement["online"]["stdDevServerTimeMs"]
            percent_dev = std_dev / measurement["online"]["serverTimeMs"] * 100
            percent_devs.append((percent_dev, scheme, scenario))
            print(f"{scheme} {scenario} {percent_dev:.2f}%")
    largest_percent_dev = max(percent_devs, key=lambda x: x[0])
    print(f"Largest percent dev: {largest_percent_dev}")
    if largest_percent_dev[0] > 10:
        print(f"WARNING: Largest percent dev: {largest_percent_dev}")

    return scheme_results, all_scenarios


def get_underlying_data_scheme(scheme):
    if scheme == "simplepir*":
        return "ypir"
    elif scheme == "doublepir*":
        return "ypir"
    elif scheme == "hintlesspir*":
        return "hintlesspir"
    else:
        return scheme


def plot_1_bit_retrieval(data_files_json: list[str], output_type: str):
    filtered_data_files_json = [x for x in data_files_json if "large-items" not in x]
    scheme_results, all_scenarios = gather_1_bit_retrieval_data(
        filtered_data_files_json
    )
    print("all_scenarios", all_scenarios)

    # schemes = ["simplepir", "doublepir", "simplepir*", "doublepir*", "ypir"]
    schemes = [
        "simplepir*",
        "doublepir*",
        "tiptoe",
        "hintlesspir",
        "ypir",
    ]

    # plot
    fig, ax = plt.subplots()
    for scheme in schemes:
        results = None
        if scheme in scheme_results:
            results = scheme_results[scheme]
        else:
            assert "*" in scheme, f"Unknown scheme: {scheme}"
            results = scheme_results[get_underlying_data_scheme(scheme)]
        xs = []
        ys = []
        for scenario in all_scenarios:
            if scenario not in results:
                continue
            measurement = results[scenario]
            db_size_gb = scenario[0] * scenario[1] / (8 * GB_SZ)
            db_size_gb_nice = scenario[0] * scenario[1] / (8 * TRUE_GB_SZ)
            xs.append(int(db_size_gb_nice))
            # throughput = db_size_gb / (measurement["online"]["serverTimeMs"] / MS_PER_S)
            cur_throughput = val_throughput(scheme, scenario, measurement)
            ys.append(cur_throughput)
        if scheme in scheme_styles:
            ax.plot(
                xs,
                ys,
                label=scheme,
                **scheme_styles[scheme],
            )
        else:
            ax.plot(xs, ys, label=scheme, marker="o")

    ax.set_xlabel("Database size (GB)")
    ax.set_ylabel("Throughput (GB/s)")

    # save
    if output_type == "tex":
        tikzplotlib.clean_figure()
        code = tikzplotlib.get_tikz_code()
        for old, new in replace.items():
            code = code.replace(old, new)
        with open("rodeo/plot.tex", "w") as fh:
            fh.write(code)
    elif output_type == "pdf":
        plt.savefig("rodeo/plot.pdf")
    else:
        print(f"Unknown output type: {output_type}", file=sys.stderr)
        sys.exit(1)


"""
Model table:
      \textbf{Database} & \textbf{Metric} & \textbf{Best Previous} & {\bf \scshape Spiral} & {\bf \scshape SpiralStream} & {\bf \scshape SpiralPack} & {\bf \scshape SpiralStreamPack} \\ \midrule
                                          & {\bf Param. Size}   & 1 MB            & 14 MB            & \emphcell{344 KB} & 14 MB             & 16 MB                \\
                                          & {\bf Query Size}    & 34 MB           & \emphcell{14 KB} & 8 MB              & \emphcell{14 KB}  & 15 MB                \\
 $\mathbf{2^{20} \times 256} \textbf{B}$  & {\bf Response Size} & 66 KB           & 21 KB            & \emphcell{20 KB}  & \emphcell{20 KB}  & 71 KB                \\
 {\bf (268 MB)}                           & {\bf Computation}   & 1.44 s          & 1.68 s           & 0.86 s            & 1.37 s            & \emphcell{0.42 s}    \\ \cmidrule{2-7}
                                          & {\bf Rate}          & 0.0039          & 0.0122           & \emphcell{0.0125} & \emphcell{0.0125} & 0.0036               \\
                                          & {\bf Throughput}    & 186 MB/s        & 159 MB/s         & 312 MB/s          & 196 MB/s          & \emphcell{635 MB/s}  \\ \midrule
                                          & {\bf Param. Size}   & 5 MB            & 18 MB            & \emphcell{3 MB}   & 18 MB             & 16 MB                \\
                                          & {\bf Query Size}    & 63 KB           & \emphcell{14 KB} & 15 MB             & \emphcell{14 KB}  & 30 MB                \\
 $\mathbf{2^{18} \times 30} \textbf{KB}$  & {\bf Response Size} & 127 KB          & 84 KB            & \emphcell{62 KB}  & 86 KB             & 96 KB                \\
 {\bf (7.9 GB)}                           & {\bf Computation}   & 52.99 s         & 24.52 s          & 9.00 s            & 17.69 s           & \emphcell{5.33 s}    \\ \cmidrule{2-7}
                                          & {\bf Rate}          & 0.2363          & 0.3573           & \emphcell{0.4803} & 0.3488            & 0.3117               \\
                                          & {\bf Throughput}    & 148 MB/s        & 321 MB/s         & 874 MB/s          & 444 MB/s          & \emphcell{1.48 GB/s} \\ \midrule
                                          & {\bf Param. Size}   & 5 MB            & 17 MB            & \emphcell{1 MB}   & 47 MB             & 24 MB                \\
                                          & {\bf Query Size}    & 63 KB           & \emphcell{14 KB} & 8 MB              & \emphcell{14 KB}  & 30 MB                \\
 $\mathbf{2^{14} \times 100} \textbf{KB}$ & {\bf Response Size} & 508 KB          & 242 KB           & 208 KB            & 188 KB            & \emphcell{150 KB}    \\
 {\bf (1.6 GB)}                           & {\bf Computation}   & 14.35 s         & 4.92 s           & 2.40 s            & 4.58 s            & \emphcell{1.21 s}    \\ \cmidrule{2-7}
                                          & {\bf Rate}          & 0.1969          & 0.4129           & 0.4811            & 0.5307            & \emphcell{0.6677}    \\
                                          & {\bf Throughput}    & 114 MB/s        & 333 MB/s         & 683 MB/s          & 358 MB/s          & \emphcell{1.35 GB/s} \\

"""


DASH = "---"


def format_bytes(n: int):
    """Takes a number of bytes and returns a nicely formatted string. Round."""
    if n == 0:
        return DASH
    if n < KB_SZ:
        return f"{n} B"
    elif n < MB_SZ:
        return f"{n / KB_SZ:.0f} KB"
    elif n < GB_SZ:
        if n / MB_SZ < 10:
            return f"{n / (MB_SZ):.1f} MB"
        else:
            return f"{n / (MB_SZ):.0f} MB"
    elif n < TB_SZ:
        if n / GB_SZ < 10:
            return f"{n / (GB_SZ):.1f} GB"
        else:
            return f"{n / (GB_SZ):.0f} GB"
    else:
        return f"{n / (TB_SZ):.1f} TB"


def format_time(scenario: tuple, server_time_ms: int) -> str:
    if server_time_ms < MS_PER_S:
        return f"{server_time_ms:.0f} ms"
    elif server_time_ms <= MS_PER_S * 60 * 60 * 5:
        return f"{server_time_ms / MS_PER_S:.2f} s"
    else:
        return f"{server_time_ms / (MS_PER_S * 60 * 60):.0f} h"


def alt_format_time(scenario: tuple, server_time_ms: int) -> str:
    if server_time_ms <= MS_PER_S * 60 * 60 * 5:
        return f"{server_time_ms / MS_PER_S:.2f} s"
    else:
        return f"{server_time_ms / (MS_PER_S * 60 * 60):.0f} h"


def calc_tput(scenario: tuple, server_time_ms: int):
    if server_time_ms == 0:
        return DASH
    db_size_gb = scenario[0] * scenario[1] / (8 * GB_SZ)
    throughput = db_size_gb / (server_time_ms / MS_PER_S)
    return throughput


def format_tput(scenario: tuple, server_time_ms: int) -> str:
    if server_time_ms == 0:
        return DASH
    db_size_gb = scenario[0] * scenario[1] / (8 * GB_SZ)
    throughput = db_size_gb / (server_time_ms / MS_PER_S)

    # nicely format throughput, down to mb
    if throughput < 1:
        if throughput * BASE_SZ < 10:
            return f"{throughput * BASE_SZ:.1f} MB/s"
        else:
            return f"{throughput * BASE_SZ:.0f} MB/s"
    else:
        return f"{throughput:.1f} GB/s"


def pad(s: str, width: int) -> str:
    assert len(s) < width
    return s + " " * (width - len(s))


def calc_download(scheme, scenario, x):
    if "simplepir" in scheme or "doublepir" in scheme:
        # fixes minor bug in collection of stats
        return x["online"]["downloadBytes"] - x["online"]["uploadBytes"]
    else:
        return x["online"]["downloadBytes"]


def prep_tput(scheme, scenario, x):
    if scheme == "simplepir*":
        return format_tput(scenario, x["offline"]["simplepirPrepTimeMs"])
    else:
        return format_tput(scenario, x["offline"]["serverTimeMs"])


def off_download(scheme, scenario, x):
    if scheme == "simplepir*":
        return format_bytes(x["offline"]["simplepirHintBytes"])
    elif scheme == "doublepir*":
        return format_bytes(x["offline"]["doublepirHintBytes"])
    else:
        return format_bytes(x["offline"]["downloadBytes"])


def val_off_download(scheme, scenario, x):
    if scheme == "simplepir*":
        return x["offline"]["simplepirHintBytes"]
    elif scheme == "doublepir*":
        return x["offline"]["doublepirHintBytes"]
    else:
        return x["offline"]["downloadBytes"]


def upload(scheme, scenario, x):
    if scheme == "simplepir*":
        return format_bytes(x["online"]["simplepirQueryBytes"])
    elif scheme == "doublepir*":
        return format_bytes(
            x["online"]["simplepirQueryBytes"] + x["online"]["doublepirQueryBytes"]
        )
    else:
        return format_bytes(x["online"]["uploadBytes"])


def download(scheme, scenario, x):
    if scheme == "simplepir*":
        return format_bytes(x["online"]["simplepirRespBytes"])
    elif scheme == "doublepir*":
        return format_bytes(x["online"]["doublepirRespBytes"])
    else:
        return format_bytes(calc_download(scheme, scenario, x))


def val_download(scheme, scenario, x):
    if scheme == "simplepir*":
        return x["online"]["simplepirRespBytes"]
    elif scheme == "doublepir*":
        return x["online"]["doublepirRespBytes"]
    else:
        return calc_download(scheme, scenario, x)


def val_upload(scheme, scenario, x):
    if scheme == "simplepir*":
        return x["online"]["simplepirQueryBytes"]
    elif scheme == "doublepir*":
        return x["online"]["simplepirQueryBytes"] + x["online"]["doublepirQueryBytes"]
    else:
        return x["online"]["uploadBytes"]


def server_time(scheme, scenario, x):
    if scheme == "hintlesspir*":
        return format_time(
            scenario, x["online"]["firstPassTimeMs"] + x["online"]["ringPackingTimeMs"]
        )
    elif scheme == "simplepir*":
        return format_time(scenario, x["online"]["firstPassTimeMs"])
    elif scheme == "doublepir*":
        return format_time(
            scenario, x["online"]["firstPassTimeMs"] + x["online"]["secondPassTimeMs"]
        )
    else:
        return format_time(scenario, x["online"]["serverTimeMs"])


def alt_server_time(scheme, scenario, x):
    return alt_format_time(scenario, x["online"]["serverTimeMs"])


def val_throughput(scheme, scenario, x):
    if scheme == "hintlesspir*":
        return calc_tput(
            scenario, x["online"]["firstPassTimeMs"] + x["online"]["ringPackingTimeMs"]
        )
    elif scheme == "simplepir*":
        return calc_tput(scenario, x["online"]["firstPassTimeMs"])
    elif scheme == "doublepir*":
        return calc_tput(
            scenario, x["online"]["firstPassTimeMs"] + x["online"]["secondPassTimeMs"]
        )
    else:
        return calc_tput(scenario, x["online"]["serverTimeMs"])


def throughput(scheme, scenario, x):
    print("throughput", scheme)
    if scheme == "hintlesspir*":
        return format_tput(
            scenario, x["online"]["firstPassTimeMs"] + x["online"]["ringPackingTimeMs"]
        )
    elif scheme == "simplepir*":
        return format_tput(scenario, x["online"]["firstPassTimeMs"])
    elif scheme == "doublepir*":
        return format_tput(
            scenario, x["online"]["firstPassTimeMs"] + x["online"]["secondPassTimeMs"]
        )
    else:
        return format_tput(scenario, x["online"]["serverTimeMs"])


def rate(scheme, scenario, x):
    # ratio of plaintext item size to ciphertext item size
    plaintext_item_size = 256  # scenario[1] / 8  # bits to bytes
    ciphertext_item_size = x["online"]["downloadBytes"]
    rate_val = plaintext_item_size / ciphertext_item_size
    return f"{rate_val:.4f}"


def table_1_bit_retrieval(args, data_files_json: list[str], output_type: str):
    disp_width = 30
    padw = lambda s: pad(s, disp_width)
    scheme_results, all_scenarios = gather_1_bit_retrieval_data(data_files_json)
    schemes = ["simplepir", "doublepir", "tiptoe", "hintlesspir", "ypir"]
    nice_schemes = ["SimplePIR", "DoublePIR", "Tiptoe", "HintlessPIR", "YPIR"]
    if args.star_variants:
        schemes = [
            "simplepir",
            "simplepir*",
            "doublepir",
            "doublepir*",
            "hintlesspir",
            "ypir",
        ]
        nice_schemes = [
            "SimplePIR",
            "SimplePIR*",
            "DoublePIR",
            "DoublePIR*",
            "\\hintlesspir",
            "YPIR",
        ]
        # schemes = ["simplepir*", "doublepir*", "ypir"]
        # nice_schemes = ["SimplePIR*", "DoublePIR*", "YPIR"]
    if args.ypir_only:
        schemes = ["ypir"]
        nice_schemes = ["YPIR"]

    db_scenarios = {
        "1 GB": (8589934592, 1),
        "8 GB": (68719476736, 1),
        "16 GB": (137438953472, 1),
        "32 GB": (274877906944, 1),
    }
    db_scenarios_keys = ["1 GB", "8 GB", "32 GB"]
    # db_scenarios_keys = ["1 GB", "8 GB", "16 GB", "32 GB"]
    metrics = {
        "Prep. Throughput": prep_tput,
        "Off. Download": off_download,
        "Upload": upload,
        "Download": download,
        "Server Time": server_time,
        "Throughput": throughput,
    }
    metric_keys = [
        "Prep. Throughput",
        "Off. Download",
        "Upload",
        "Download",
        "Server Time",
        "Throughput",
    ]
    if args.star_variants:
        metrics = {
        "Prep. Speed": prep_tput,
        "Off. Comm.": off_download,
        "Upload": upload,
        "Download": download,
        "Server Time": server_time,
        "Throughput": throughput,
        }
        metric_keys = [
            "Prep. Speed",
            "Off. Comm.",
            "Upload",
            "Download",
            "Server Time",
            "Throughput",
        ]
    output = (
        padw("\\textbf{Database}")
        + "& "
        + padw("\\textbf{Metric}")
        + "& "
        + "& ".join([padw("\\textbf{" + scheme + "}") for scheme in nice_schemes])
        + "\\\\ \\midrule\n"
    )

    for db_sz in db_scenarios_keys:
        scenario = db_scenarios[db_sz]
        for metric in metric_keys:
            row = "& " + padw("{\\bf " + metric + "}")
            if metric == "Upload":
                row = padw("{\\bf " + db_sz + "}") + row
            else:
                row = padw("") + row
            for scheme in schemes:
                if "hintlesspir" in scheme and metric == "Off. Download":
                    row += "& " + padw(DASH)
                    continue

                # if "hintlesspir" in scheme and scenario not in scheme_results[scheme]:
                #     char = "-" if metric == "Off. Download" else "$\\ast$"
                #     row += "& " + padw(char)
                #     continue

                res = None
                if scheme not in scheme_results:
                    if "*" in scheme:
                        # use stripped name for *-variants
                        stripped_scheme = scheme[:-1]
                        if stripped_scheme != "hintlesspir":
                            # use YPIR measurements for simplepir* and doublepir*
                            stripped_scheme = "ypir"
                        res = scheme_results[stripped_scheme][scenario]
                    else:
                        assert False, f"Unknown scheme: {scheme}"
                elif scenario not in scheme_results[scheme]:
                    row += "& " + padw(DASH)
                    continue
                else:
                    res = scheme_results[scheme][scenario]

                row += "& " + padw(metrics[metric](scheme, scenario, res))
            row += "\\\\"
            if metric == "Off. Download":  # or metric == "Server Time":
                row += f" \\cmidrule{{2-{len(schemes) + 2}}}"
            if metric == "Throughput":
                if db_sz != db_scenarios_keys[-1]:
                    row += " \\midrule"
                else:
                    row += " \\bottomrule"
            output += row + "\n"

    print(output)


def plot_large_items(args, data_files_json: list[str], output_type: str):
    scheme_results, all_scenarios = gather_1_bit_retrieval_data(data_files_json)
    print(scheme_results)
    # assert output_type == "tex"

    schemes = ["simplepir", "hintlesspir", "ypir-sp"]
    nice_schemes = ["SimplePIR", "HintlessPIR", "\\YPIRSP"]
    if args.star_variants:
        schemes = ["simplepir", "hintlesspir", "hintlesspir*", "ypir-sp"]
        nice_schemes = ["SimplePIR", "HintlessPIR", "HintlessPIR+", "\\YPIRSP"]

    if args.respire:
        schemes = ["ypir-sp"]
        nice_schemes = ["YPIR+SimplePIR"]

    # generate table
    disp_width = 45
    padw = lambda s: pad(s, disp_width)
    output = (
        padw("\\textbf{Database}")
        + "& "
        + padw("\\textbf{Metric}")
        + "& "
        + "& ".join([padw("\\textbf{" + scheme + "}") for scheme in nice_schemes])
        + "\\\\ \\midrule\n"
    )

    db_keys = [
        "$\\boldsymbol{2^{15} \\times 32}$ {\\bf KB}",
        "$\\boldsymbol{2^{18} \\times 32}$ {\\bf KB}",
        "$\\boldsymbol{2^{19} \\times 64}$ {\\bf KB}",
    ]

    db_scenarios = {
        db_keys[0]: (32768, 262144),
        db_keys[1]: (131072, 524288),
        db_keys[2]: (262144, 1048576),
    }

    db_caption = {
        db_keys[0]: "{\\bf (1 GB)}",
        db_keys[1]: "{\\bf (8 GB)}",
        db_keys[2]: "{\\bf (32 GB)}",
    }

    alt_scenarios = {
        (32768, 262144): (8589934592, 1),
        (131072, 524288): (68719476736, 1),
        (262144, 1048576): (274877906944, 1),
    }

    metrics = {
        "Prep. Throughput": prep_tput,
        "Off. Download": off_download,
        "Upload": upload,
        "Download": download,
        "Server Time": server_time,
        "Throughput": throughput,
    }

    metric_keys = [
        "Prep. Throughput",
        "Off. Download",
        "Upload",
        "Download",
        "Server Time",
        "Throughput",
    ]

    if args.respire:
        db_keys = [
            "$2^{14} \\times 16\ \\text{KB}$",
            "$2^{15} \\times 32\ \\text{KB}$",
            "$2^{17} \\times 64\ \\text{KB}$",
        ]
        db_scenarios = {
            db_keys[0]: (16384, 131072),
            db_keys[1]: (32768, 262144),
            db_keys[2]: (131072, 524288),
        }
        db_caption = {
            db_keys[0]: "(256 MB)",
            db_keys[1]: "(1 GB)",
            db_keys[2]: "(8 GB)",
        }
        metrics["Rate"] = rate
        metrics["Server Time"] = alt_server_time
        metric_keys = [
            "Off. Download",
            "Upload",
            "Download",
            "Server Time",
            "Rate",
            "Throughput",
        ]

    for db_sz in db_keys:
        scenario = db_scenarios[db_sz]
        for metric in metric_keys:
            row = "& " + padw("{\\bf " + metric + "}")
            if metric == "Upload":
                row = padw("{ " + db_sz + "}") + row
            elif metric == "Download":
                row = padw(db_caption[db_sz]) + row
            else:
                row = padw("") + row
            for scheme in schemes:
                if scheme == "hintlesspir" and metric == "Off. Download":
                    row += "& " + padw(DASH)
                    continue

                if scheme == "hintlesspir" and scenario not in scheme_results[scheme]:
                    row += "& " + padw(DASH)
                    continue

                # if scheme == "hintlesspir" and (
                #     metric == "Throughput"
                #     or metric == "Prep. Throughput"
                #     or metric == "Server Time"
                # ):
                #     row += "& " + padw(DASH)
                #     continue

                res = None

                # if scheme not in scheme_results:
                #     assert False, f"Unknown scheme: {scheme}"

                stripped_scheme = scheme
                if "*" in scheme:
                    # use stripped name for *-variants
                    stripped_scheme = scheme[:-1]

                if (
                    scenario not in scheme_results[stripped_scheme]
                    and scenario not in alt_scenarios
                ):
                    row += "& " + padw(DASH)
                    continue

                maybe_scenario = scenario
                if scenario not in scheme_results[stripped_scheme]:
                    maybe_scenario = alt_scenarios[scenario]

                res = scheme_results[stripped_scheme][maybe_scenario]

                row += "& " + padw(metrics[metric](scheme, maybe_scenario, res))
            row += "\\\\"
            if metric == "Off. Download":
                row += f" \\cmidrule{{2-{len(schemes) + 2}}}"

            if metric == "Throughput":
                if db_sz != list(db_scenarios.keys())[-1]:
                    row += " \\midrule"
                else:
                    row += " \\bottomrule"
            output += row + "\n"

    print(output)


def plot_comm_comp_tradeoff(args, data_files_json: list[str], output_type: str):
    # first get all the data
    scheme_results, all_scenarios = gather_1_bit_retrieval_data(data_files_json)

    schemes = ["simplepir*", "doublepir*", "hintlesspir", "ypir"]

    # we have 1, 2, 4, 8, 16, and 32 GB
    # we'll consider 32 x 1 GB, 16 x 2 GB, etc
    # we want to gather the (throughput, total comm.) pairs for every instantiation
    # then plot in a scatter
    db_sz_to_scenario = {
        0.125: (1073741824, 1),
        0.25: (2147483648, 1),
        0.5: (4294967296, 1),
        1.0: (8589934592, 1),
        2.0: (17179869184, 1),
        4.0: (34359738368, 1),
        8.0: (68719476736, 1),
        16.0: (137438953472, 1),
        32.0: (274877906944, 1),
    }
    db_configs = [(32 / x, x) for x in db_sz_to_scenario.keys()]
    print(db_configs)
    off_comms = {}
    scheme_points = {}
    for scheme in schemes:
        scheme_points[scheme] = []
        off_comms[scheme] = []
    for instances, db_sz in db_configs:
        instances = int(instances)
        scenario = db_sz_to_scenario[db_sz]
        for scheme in schemes:
            scheme_to_use = get_underlying_data_scheme(scheme)
            if scenario not in scheme_results[scheme_to_use]:
                continue
            res = scheme_results[scheme_to_use][scenario]
            throughput = val_throughput(scheme, scenario, res)
            up = val_upload(scheme, scenario, res)
            down = val_download(scheme, scenario, res)
            total_comm = up + instances * down
            scheme_points[scheme].append((throughput, total_comm))

            off_comm = instances * val_off_download(scheme, scenario, res)
            off_comms[scheme].append(off_comm)

    for scheme in schemes:
        scheme_points[scheme] = list(set(scheme_points[scheme]))
    print(scheme_points)
    # interesting bit: remove 'dominated' points
    # for a given scheme, if there's a point
    # that has higher throughput *and* lower
    # comm. than another point, remove the latter
    remove_dominated = True
    if remove_dominated:
        dominates = lambda p1, p2: p1[0] >= p2[0] and p1[1] <= p2[1]
        for scheme in schemes:
            while True:
                points = scheme_points[scheme]
                to_remove_idxs = set()
                for i in range(len(points)):
                    for j in range(len(points)):
                        if i == j:
                            continue
                        if dominates(points[i], points[j]):
                            print("Removing " + scheme + "  " + str(points[j]))
                            print("(dominated by " + str(points[i]) + ")")
                            to_remove_idxs.add(j)
                if len(to_remove_idxs) == 0:
                    break
                to_remove_idxs = list(set(to_remove_idxs))
                to_remove_idxs.sort(reverse=True)
                for idx in to_remove_idxs:
                    del scheme_points[scheme][idx]
                    del off_comms[scheme][idx]
            scheme_points[scheme].sort()

    for scheme in schemes:
        print("Max off. comm. for " + scheme + ": " + str(max(off_comms[scheme])))
        print("Min off. comm. for " + scheme + ": " + str(min(off_comms[scheme])))

    # matplotlib lines
    fig, ax = plt.subplots()
    for scheme in schemes:
        xs = [x[1] for x in scheme_points[scheme]]
        ys = [x[0] for x in scheme_points[scheme]]
        # ax.scatter(xs, ys, label=scheme, **scheme_styles[scheme])
        # plot lines with dots
        ax.plot(xs, ys, label=scheme, **scheme_styles[scheme], marker="o")

    ax.set_xlabel("Total Communication (B)")
    ax.set_ylabel("Throughput (GB/s)")
    # add legend in south west
    ax.legend(loc="lower left")
    # set y axis to start at 0
    ax.set_ylim(bottom=0)
    # invert the x axis (comm)
    ax.invert_xaxis()
    plt.savefig("rodeo/plot.pdf")

    for scheme in schemes:
        print()
        print("\\addplot [" + scheme + "]\ntable {%")
        for point in scheme_points[scheme]:
            print(f"{point[1]/(1<<20)} {point[0]}")
        print("")


def plot_ccb(data_files_json: list[str], output_type: str):
    # load data, each file has a list of results for a scheme
    data = []
    for data_file in data_files_json:
        with open(data_file, "r") as fh:
            data.append(json.load(fh))

    # get all scenarios
    # scenarios = []
    # for datum in data:
    #     for result in datum["results"]:
    #         scenario = result["scenario"]
    #         scenarios.append(scenario)
    # print("scenarios", scenarios)

    # only use the 32 GB scenario
    base_scenario = (274877906944, 1)
    scenarios = []
    for num_clients in [1, 2, 4, 8]:
        scenario = {
            "db": {
                "numItems": base_scenario[0],
                "itemSizeBits": base_scenario[1],
            },
            "clients": {"numClients": num_clients},
        }
        scenarios.append(scenario)

    db_size_gb = (
        scenarios[0]["db"]["numItems"]
        * scenarios[0]["db"]["itemSizeBits"]
        / (8 * GB_SZ)
    )
    print("db_size_gb", db_size_gb)

    # dedup scenarios by numClients
    get_clients = lambda s: s["clients"]["numClients"] if "clients" in s else 1
    scenarios = list(
        {get_clients(scenario): scenario for scenario in scenarios}.values()
    )

    # get result of each scheme for each scenario
    scheme_results = {}
    for datum in data:
        scheme = datum["scheme"]["variant"]
        if scheme not in scheme_results:
            scheme_results[scheme] = {}
        for result in datum["results"]:
            scenario = result["scenario"]
            measurement = result["measurement"]
            scheme_results[scheme][get_clients(scenario)] = measurement

    # get all scenarios
    all_scenarios = None
    for scheme, results in scheme_results.items():
        if all_scenarios is None:
            all_scenarios = set(results.keys())
        else:
            all_scenarios = all_scenarios.union(results.keys())
    all_scenarios = sorted(list(all_scenarios))

    schemes = ["simplepir*", "doublepir*", "hintlesspir*", "ypir"]

    # plot
    fig, ax = plt.subplots()
    for scheme in schemes:
        results = None
        if scheme in scheme_results:
            results = scheme_results[scheme]
        else:
            assert "*" in scheme, f"Unknown scheme: {scheme}"
            results = scheme_results["ypir"]  # on purpose
        xs = []
        ys = []
        for scenario in all_scenarios:
            if scenario not in results:
                continue
            measurement = results[scenario]
            xs.append(scenario)
            scenario_traditional = (scenario * db_size_gb * GB_SZ * 8, 1)
            print(scenario, scenario_traditional)
            throughput = val_throughput(scheme, scenario_traditional, measurement)
            if scheme == "hintlesspir*":
                print(scheme_results["hintlesspir"])
                total_packing_time = (
                    scenario
                    * scheme_results["hintlesspir"][1]["online"]["ringPackingTimeMs"]
                )
                print("total_packing_time", total_packing_time)
                matmul_time = scheme_results["ypir"][scenario]["online"]["serverTimeMs"]
                total_time = total_packing_time + matmul_time
                throughput = scenario * db_size_gb / (total_time / MS_PER_S)
            print(
                "!",
                scheme,
                throughput,
                scenario,
                db_size_gb,
                measurement["online"]["serverTimeMs"],
            )
            print(
                scheme,
                scenario,
                throughput,
                db_size_gb,
                measurement["online"]["serverTimeMs"],
            )
            ys.append(throughput)
        print(scheme, list(xs), list(ys))
        if scheme in scheme_styles:
            ax.plot(
                xs,
                ys,
                label=scheme,
                **scheme_styles[scheme],
            )
        else:
            ax.plot(xs, ys, label=scheme, marker="o")

    ax.set_xlabel("Number of clients to batch across")
    ax.set_ylabel("Effective throughput per query (GB/s)")

    # save
    if output_type == "tex":
        tikzplotlib.clean_figure()
        code = tikzplotlib.get_tikz_code()
        for old, new in replace.items():
            code = code.replace(old, new)
        with open("rodeo/plot.tex", "w") as fh:
            fh.write(code)
    elif output_type == "pdf":
        plt.savefig("rodeo/plot.pdf")
    else:
        print(f"Unknown output type: {output_type}", file=sys.stderr)
        sys.exit(1)


def plot_ccb_rlwe(data_files_json: list[str], output_type: str):
    # load data, each file has a list of results for a scheme
    data = []
    for data_file in data_files_json:
        with open(data_file, "r") as fh:
            data.append(json.load(fh))

    # get all scenarios
    scenarios = []
    for datum in data:
        for result in datum["results"]:
            scenario = result["scenario"]
            scenarios.append(scenario)

    db_size_gb = (
        scenarios[0]["db"]["numItems"]
        * scenarios[0]["db"]["itemSizeBits"]
        / (8 * GB_SZ)
    )

    # dedup scenarios by numClients
    scenarios = list(
        {scenario["clients"]["numClients"]: scenario for scenario in scenarios}.values()
    )

    # get result of each scheme for each scenario
    scheme_results = {}
    for datum in data:
        scheme = datum["scheme"]["variant"]
        if scheme not in scheme_results:
            scheme_results[scheme] = {}
        for result in datum["results"]:
            scenario = result["scenario"]
            measurement = result["measurement"]
            scheme_results[scheme][scenario["clients"]["numClients"]] = measurement

    # get all scenarios
    all_scenarios = None
    for scheme, results in scheme_results.items():
        if all_scenarios is None:
            all_scenarios = set(results.keys())
        else:
            all_scenarios = all_scenarios.union(results.keys())
    all_scenarios = sorted(list(all_scenarios))

    # plot
    fig, ax = plt.subplots()
    for scheme, results in scheme_results.items():
        xs = []
        ys = []
        for scenario in all_scenarios:
            if scenario not in results:
                continue
            measurement = results[scenario]
            xs.append(scenario)
            throughput = (
                scenario
                * db_size_gb
                / (measurement["online"]["serverTimeMs"] / MS_PER_S)
            )
            print(
                scheme,
                scenario,
                throughput,
                db_size_gb,
                measurement["online"]["serverTimeMs"],
            )
            ys.append(throughput)
        if scheme in scheme_styles:
            ax.plot(
                xs,
                ys,
                label=scheme,
                **scheme_styles[scheme],
            )
        else:
            ax.plot(xs, ys, label=scheme, marker="o")

    ax.set_xlabel("Number of clients to batch across")
    ax.set_ylabel("Effective throughput per query (GB/s)")

    # save
    if output_type == "tex":
        tikzplotlib.clean_figure()
        code = tikzplotlib.get_tikz_code()
        for old, new in replace.items():
            code = code.replace(old, new)
        with open("rodeo/plot.tex", "w") as fh:
            fh.write(code)
    elif output_type == "pdf":
        plt.savefig("rodeo/plot.pdf")
    else:
        print(f"Unknown output type: {output_type}", file=sys.stderr)
        sys.exit(1)


def custom_disp_seconds(t: float, ms_cutoff: float) -> str:
    if t < ms_cutoff:
        # use ms
        return f"{t * 1000:.0f} ms"
    else:
        return f"{t:.2f} s"


def ypir_breakdown(data_files_json: list[str], output_type: str):
    scheme_results, all_scenarios = gather_1_bit_retrieval_data(data_files_json)

    # plot stacked bar chart
    fig, ax = plt.subplots()
    scheme = "ypir"
    results = scheme_results[scheme]

    xs = []
    ys_first = []
    ys_second = []
    ys_ring = []
    for scenario in all_scenarios:
        if scenario not in results:
            continue
        measurement = results[scenario]
        db_size_gb = scenario[0] * scenario[1] / (8 * GB_SZ)
        xs.append(db_size_gb)

        # throughput = db_size_gb / (measurement["online"]["serverTimeMs"] / MS_PER_S)
        first_pass = measurement["online"]["firstPassTimeMs"]
        second_pass = measurement["online"]["secondPassTimeMs"]
        ring_packing = measurement["online"]["ringPackingTimeMs"]

        ys_first.append(first_pass)
        ys_second.append(second_pass)
        ys_ring.append(ring_packing)

    for i in range(len(xs)):
        xs[i] = math.log2(xs[i])

    # generate table like this:
    """
    {\bf Size} & {\bf SimplePIR} & {\bf DoublePIR} & {\bf Packing} & {\bf Total} \\ \midrule
    1~GB       & 0.08 s (18\%)   & 0.02 s (4\%)   & 0.35 s (78\%) & 0.45 s       \\
    4~GB       & 0.34 s (47\%)   & 0.03 s (4\%)   & 0.35 s (49\%) & 0.72 s       \\
    8~GB       & 0.73 s (66\%)   & 0.03 s (3\%)   & 0.34 s (31\%) & 1.1  s       \\
    16~GB      & 1.46 s (78\%)   & 0.06 s (3\%)   & 0.34 s (18\%) & 1.9  s       \\
    32~GB      & 3.04 s (88\%)   & 0.06 s (2\%)   & 0.34 s (10\%) & 3.4  s       \\
    """

    disp_width = 16
    padw = lambda s: pad(s, disp_width)

    output = (
        padw("{\\bf Size}")
        + " & "
        + padw("{\\bf SimplePIR}")
        + " & "
        + padw("{\\bf DoublePIR}")
        + " & "
        + padw("{\\bf Packing}")
        + " & "
        + padw("{\\bf Total}")
        + " \\\\ \\midrule\n"
    )

    for i in range(len(xs)):
        db_size_gb = round(all_scenarios[i][0] * all_scenarios[i][1] / (8 * GB_SZ))
        first_pass = ys_first[i] / MS_PER_S
        second_pass = ys_second[i] / MS_PER_S
        ring_packing = ys_ring[i] / MS_PER_S
        total = first_pass + second_pass + ring_packing
        output += (
            padw(f"{db_size_gb}~GB")
            + " & "
            + padw(
                f"{custom_disp_seconds(first_pass, 0.0)} ({first_pass / total * 100:.0f}\\%)"
            )
            + " & "
            + padw(
                f"{custom_disp_seconds(second_pass, 1.0)} ({second_pass / total * 100:.0f}\\%)"
            )
            + " & "
            + padw(
                f"{custom_disp_seconds(ring_packing, 1.0)} ({ring_packing / total * 100:.0f}\\%)"
            )
            + " & "
            + padw(f"{total:.2f} s")
            + " \\\\ \n"
        )

    print(output)


def query_breakdown(data_files_json: list[str], output_type: str):
    scheme_results, all_scenarios = gather_1_bit_retrieval_data(data_files_json)

    # plot stacked bar chart
    fig, ax = plt.subplots()
    scheme = "ypir"
    results = scheme_results[scheme]

    xs = []
    ys_first = []
    ys_second = []
    ys_ring = []
    for scenario in all_scenarios:
        if scenario not in results:
            continue
        measurement = results[scenario]
        db_size_gb = scenario[0] * scenario[1] / (8 * GB_SZ)
        xs.append(db_size_gb)

        # throughput = db_size_gb / (measurement["online"]["serverTimeMs"] / MS_PER_S)
        # total query size
        query_size = measurement["online"]["uploadBytes"]
        simplepir_query_size = measurement["online"]["simplepirQueryBytes"]

        # correction for 7/8 difference in query size computation (conservative):
        doublepir_query_size = measurement["online"]["doublepirQueryBytes"] * (8 / 7)
        key_switching_matrices_size = query_size - (
            simplepir_query_size + doublepir_query_size
        )

        ys_first.append(simplepir_query_size)
        ys_second.append(doublepir_query_size)
        ys_ring.append(key_switching_matrices_size)

    for i in range(len(xs)):
        xs[i] = math.log2(xs[i])

    disp_width = 20
    padw = lambda s: pad(s, disp_width)

    output = (
        padw("{\\bf Database Size}")
        + " & "
        + padw("{\\bf SimplePIR}")
        + " & "
        + padw("{\\bf DoublePIR}")
        + " & "
        + padw("{\\bf Packing}")
        + " & "
        + padw("{\\bf Total}")
        + " \\\\ \\midrule\n"
    )

    for i in range(len(xs)):
        db_size_gb = round(all_scenarios[i][0] * all_scenarios[i][1] / (8 * GB_SZ))
        first_pass = ys_first[i]
        second_pass = ys_second[i]
        ring_packing = ys_ring[i]
        total = first_pass + second_pass + ring_packing
        output += (
            padw(f"{db_size_gb}~GB")
            + " & "
            + padw(f"{format_bytes(first_pass)} ({first_pass / total * 100:.0f}\\%)")
            + " & "
            + padw(f"{format_bytes(second_pass)} ({second_pass / total * 100:.0f}\\%)")
            + " & "
            + padw(
                f"{format_bytes(ring_packing)} ({ring_packing / total * 100:.0f}\\%)"
            )
            + " & "
            + padw(f"{format_bytes(total)}")
            + " \\\\ \n"
        )

    print(output)


queries_per_week = 10**4 / 500
assert queries_per_week == 20.0
cost_per_byte_down = 0.09 / 1e9
# cost_per_cpu_s = 1.829 / (60 * 60 * 64)  # r6i.16xlarge 3y Reserved Instance price
cost_per_cpu_s = 1.5e-5  # cost quoted in original SimplePIR paper
num_clients = 1


def get_variant_multiple(variant) -> float:
    if variant is None:
        return 1.0

    if variant == "m":  # monthly
        assert False, "Monthly not supported"
    elif variant == "w":  # weekly
        return 1.0
    elif variant == "d":  # daily
        return 27.95 / 4.0
    elif variant == "h":  # hourly
        raise Exception(f"Unknown variant: {variant}")
    else:
        raise Exception(f"Unknown variant: {variant}")


def off_download_cost(scheme, scenario, x) -> float:
    cost = (
        get_variant_multiple(get_variant(scheme)) * num_clients * x * cost_per_byte_down
    )
    return cost


def off_download_bytes(scheme, scenario, x) -> float:
    cost = get_variant_multiple(get_variant(scheme)) * num_clients * x
    return cost


def download_cost(scheme, scenario, x) -> float:
    cost = num_clients * x * queries_per_week * cost_per_byte_down
    return cost


def cpu_cost(x_ms):
    cost = num_clients * x_ms / MS_PER_S * queries_per_week * cost_per_cpu_s
    return cost


def cost_str(cost):
    return f"\\${cost:.6f}"


def get_variant(scheme):
    if "-" in scheme:
        scheme, variant = scheme.split("-")
        return variant
    else:
        return None


def plot_sct(data_files_json: list[str], output_type: str):
    disp_width = 32
    padw = lambda s: pad(s, disp_width)

    scheme_results, all_scenarios = gather_1_bit_retrieval_data(data_files_json)

    schemes = [
        "doublepir-w",
        "doublepir-d",
        "tiptoe",
        "hintlesspir",
        "ypir",
    ]
    nice_schemes = [
        "DoublePIR",
        "DoublePIR",
        "Tiptoe",
        "HintlessPIR",
        "YPIR",
    ]
    nice_freqs = [
        "Weekly",
        "Daily",
        DASH,
        DASH,
        DASH,
    ]

    scenario = (68719476736, 1)  # 8 GB

    # print table
    output = (
        padw(" ")
        + "".join(
            [padw("& \\textbf{" + scheme + "}") + "  " for scheme in nice_schemes]
        )
        + "\\\\ \n"
    )
    output = (
        output
        + padw("\\textbf{Update Frequency}")
        + "".join([padw("& \emph{" + freq + "}") + "  " for freq in nice_freqs])
        + "\\\\ \\midrule\n"
    )
    metrics = {
        "Offline Download": lambda sch, s, x: format_bytes(
            off_download_bytes(sch, s, x["offline"]["downloadBytes"])
        ),
        "Upload": lambda sch, s, x: format_bytes(
            queries_per_week * x["online"]["uploadBytes"]
        ),
        "Download": lambda sch, s, x: format_bytes(
            queries_per_week * calc_download(sch, s, x)
        ),
        "Computation": lambda sch, s, x: format_time(
            s, queries_per_week * x["online"]["serverTimeMs"]
        ),
        "Communication Cost": lambda sch, s, x: cost_str(
            download_cost(sch, s, calc_download(sch, s, x))
            + off_download_cost(sch, s, x["offline"]["downloadBytes"])
        ),
        "Computation Cost": lambda sch, s, x: cost_str(
            cpu_cost(x["online"]["serverTimeMs"])
        ),
        "Total Cost": lambda sch, s, x: 0,
    }

    metric_keys = [
        "Offline Download",
        "Upload",
        "Download",
        "Computation",
        "Communication Cost",
        "Computation Cost",
        "Total Cost",
    ]

    for metric in metric_keys:
        row = padw("{\\bf " + metric + "}")
        for scheme in schemes:
            full_scheme = scheme
            if "-" in scheme:
                scheme, _ = scheme.split("-")

            # if scheme == "hintlesspir" and (
            #     metric == "Computation" or metric == "Computation Cost"
            # ):
            #     row += "& " + padw("$*$")
            #     continue

            if metric == "Offline Download" and scheme != "doublepir":
                row += "& " + padw(DASH)
                continue

            val = None
            if metric == "Total Cost":
                x = scheme_results[scheme][scenario]
                cost = (
                    download_cost(
                        full_scheme, scenario, calc_download(scheme, scenario, x)
                    )
                    + cpu_cost(x["online"]["serverTimeMs"])
                    + off_download_cost(
                        full_scheme, scenario, x["offline"]["downloadBytes"]
                    )
                )

                val = cost_str(cost)
            else:
                val = metrics[metric](
                    full_scheme, scenario, scheme_results[scheme][scenario]
                )

            row += "& " + padw(val)
        row += "\\\\"
        if (
            metric == "Computation"
            or metric == "Computation Cost"
            # or metric == "Offline Download"
        ):
            row += " \\midrule"
        output += row + "\n"

    print(output)


figures = {
    "bit-retrieval": plot_1_bit_retrieval,
    "table-bit-retrieval": table_1_bit_retrieval,
    "ccb": plot_ccb,
    "ypir-breakdown": ypir_breakdown,
    "query-breakdown": query_breakdown,
    "ccb-rlwe": plot_ccb_rlwe,
    "sct": plot_sct,
    "large-items": plot_large_items,
    "comm-comp-tradeoff": plot_comm_comp_tradeoff,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot results from PIR Rodeo experiments"
    )

    parser.add_argument(
        "figure",
        choices=figures.keys(),
        help="Figure to plot",
    )
    parser.add_argument(
        "data_files_json",
        nargs="+",
        help="JSON files containing results from RODEO experiments",
    )
    parser.add_argument(
        "--output-type",
        choices=["pdf", "tex"],
        default="tex",
        help="Output type for plot",
    )
    # add flag for simplepir* and doublepir* variants
    parser.add_argument(
        "--star-variants",
        action="store_true",
        help="Use simplepir* and doublepir* variants instead of simplepir and doublepir",
    )
    # flag for respire db sizes (smaller)
    parser.add_argument(
        "--respire",
        action="store_true",
        help="Use respire database sizes",
    )
    # add flag for ypir only
    parser.add_argument(
        "--ypir-only",
        action="store_true",
        help="Plot only the results for YPIR",
    )
    return parser.parse_args()


# ex:
# python rodeo/plot.py --output-type tex table-bit-retrieval rodeo/data/*.json rodeo/latest_data/v3-*
if __name__ == "__main__":
    args = parse_args()

    if args.figure == "bit-retrieval":
        plot_1_bit_retrieval(args.data_files_json, args.output_type)
    elif args.figure == "table-bit-retrieval":
        table_1_bit_retrieval(args, args.data_files_json, args.output_type)
    elif args.figure == "ccb":
        plot_ccb(args.data_files_json, args.output_type)
    elif args.figure == "ccb-rlwe":
        plot_ccb_rlwe(args.data_files_json, args.output_type)
    elif args.figure == "ypir-breakdown":
        ypir_breakdown(args.data_files_json, args.output_type)
    elif args.figure == "query-breakdown":
        query_breakdown(args.data_files_json, args.output_type)
    elif args.figure == "sct":
        plot_sct(args.data_files_json, args.output_type)
    elif args.figure == "large-items":
        plot_large_items(args, args.data_files_json, args.output_type)
    elif args.figure == "comm-comp-tradeoff":
        plot_comm_comp_tradeoff(args, args.data_files_json, args.output_type)
