import argparse
import json
import math
import random
from itertools import repeat

import numpy as np
from matplotlib import pyplot as plt

audits_per_month = 80

SCALE = 1

new_certs_per_day = 10e6 / SCALE

# PIR database
db_sz_bytes = 8 * (2**30)
word_sz_bits = 9
db_sz_words = db_sz_bytes * 8 / word_sz_bits
db_dim = int(math.ceil(math.sqrt(db_sz_words)) / math.sqrt(SCALE))
print(f"db_dim = {db_dim}")
hint_row_sz_bytes = 4 * 1024
hint_sz_bytes = 4 * 4 * 1024 * 1024

# poisson
seconds_per_day = 24 * 60 * 60
seconds_per_month = seconds_per_day * 30
lambda_new_certs = new_certs_per_day / seconds_per_day
lambda_audits = audits_per_month / seconds_per_month


def arrivals(lambda_: float, total_time: float):
    """Generate a list of arrivals according to a Poisson process."""
    t = 0
    arrivals = []
    while t < total_time:
        t += np.random.exponential(1 / lambda_)
        arrivals.append(t)
    return arrivals


def num_empty_bins(balls: int, bins: int) -> float:
    """Generate a random variable representing the number of empty bins."""
    # don't simulate, use probability (poisson approximation)
    lamb = balls / bins
    # for each bin, the probability that it is empty is exp(-lambda)
    # so we have the sum of bins independent Bernoulli trials, each with
    # probability exp(-lambda) of success.
    # This is normally distributed with mean bins * exp(-lambda) and variance
    # bins * exp(-lambda) * (1 - exp(-lambda))
    # We can approximate this with a Poisson distribution with the same mean
    # and variance.
    return np.random.normal(
        bins * np.exp(-lamb), bins * np.exp(-lamb) * (1 - np.exp(-lamb))
    )


def simulate_updates(
    db_updates_per_month: float,
    num_clients: int = 100,
    sim_total_times_s: float = 12 * seconds_per_month,
) -> float:
    """Simulate the average number of hint bytes downloaded per month."""

    full_set = set(range(db_dim))

    client_has = {}
    for c in range(num_clients):
        client_has[c] = set()

    # Do not simulate the database updates as a Poisson process.
    # lambda_db_updates = db_updates_per_month / (seconds_per_day * 30)
    # db_update_times = arrivals(lambda_db_updates, sim_total_times_s)

    # Instead, space db updates evenly throughout the month.
    total_db_updates = int(
        math.ceil(db_updates_per_month * sim_total_times_s / seconds_per_month)
    )
    db_update_times = list(
        np.linspace(0, sim_total_times_s, total_db_updates, endpoint=False)
    )

    # audit_times = arrivals(num_clients * lambda_audits, sim_total_times_s)

    # Select *exactly* X audits per month per client.
    audits_per_client = int(
        math.ceil(audits_per_month * sim_total_times_s / seconds_per_month)
    )
    print(f"audits_per_client = {audits_per_client}")
    assert (
        abs(
            audits_per_client
            - (audits_per_month * sim_total_times_s / seconds_per_month)
        )
        <= 1e-6
    )
    audits = []
    for c in range(num_clients):
        audit_times = [
            (random.uniform(0, sim_total_times_s), c) for _ in range(audits_per_client)
        ]
        audits.extend(audit_times)

    time_steps = sorted(list(zip(db_update_times, repeat("U"))) + audits)
    # replace any string of updates with the last update:
    new_time_steps = []
    for t, action in time_steps:
        if action == "U":
            if len(new_time_steps) > 0 and new_time_steps[-1][1] == "U":
                new_time_steps.pop()
        new_time_steps.append((t, action))
    time_steps = new_time_steps

    print("Simulating...")
    progress_every = sim_total_times_s / 100
    cur = 0

    percent_of_hint_downloaded = []

    total_downloaded = 0.0
    last_update = 0.0
    for t, action in time_steps:
        # if t > cur:
        #     print(f"t = {t / seconds_per_month} months")
        #     cur += progress_every
        if action == "U":
            # database update
            # print("database update at t = ", t / seconds_per_month)

            # new cert events are poisson
            seconds_since_last_update = t - last_update
            num_new_certs_since_last_update = np.random.poisson(
                lambda_new_certs * seconds_since_last_update
            )
            # print(
            #     f"num_new_certs_since_last_update = {num_new_certs_since_last_update}"
            # )

            if num_new_certs_since_last_update >= 5 * db_dim:
                # don't need to simulate, since w.v.h.p. all dims will be updated

                for c in range(num_clients):
                    client_has[c] = set()
            else:
                # print("updating", num_new_certs_since_last_update, "certs")
                dims_that_updated = set()
                for _ in range(num_new_certs_since_last_update):
                    dims_that_updated.add(random.randint(0, db_dim - 1))
                # num_to_update = min(
                #     db_dim
                #     - int(num_empty_bins(num_new_certs_since_last_update, db_dim)),
                #     db_dim,
                # )
                # if num_to_update < 0 or num_to_update > db_dim:
                #     print(num_to_update)
                # dims_that_updated = set(random.sample(range(db_dim), num_to_update))
                # print(len(dims_that_updated))

                for d in dims_that_updated:
                    for c in range(num_clients):
                        client_has[c].discard(d)

            last_update = t
        else:
            # audit
            client_performing_audit = action
            has_hint_for_rows = len(client_has[client_performing_audit])
            num_hint_rows_to_download = db_dim - has_hint_for_rows
            num_bytes_to_download = num_hint_rows_to_download * hint_row_sz_bytes
            # if num_bytes_to_download < hint_sz_bytes:
            #     print(
            #         f"downloading partial hint {num_bytes_to_download}/{hint_sz_bytes}"
            #     )
            num_bytes_to_download = min(num_bytes_to_download, hint_sz_bytes)
            # print(
            #     f"client {client_performing_audit} downloading {num_bytes_to_download}"
            # )
            percent_of_hint_downloaded.append(num_bytes_to_download / hint_sz_bytes)
            total_downloaded += num_bytes_to_download
            client_has[client_performing_audit] = full_set.copy()

    print(
        "percent hint downloaded on average: ",
        sum(percent_of_hint_downloaded) / len(percent_of_hint_downloaded),
    )
    # save histogram
    plt.hist(percent_of_hint_downloaded, bins=100)
    plt.xlabel("Percent of hint downloaded")
    plt.ylabel("Frequency")
    plt.savefig("percent_of_hint_downloaded.png")

    return total_downloaded / num_clients


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db-updates-per-month",
        type=float,
        default=30.0,
        help="Number of database updates per month",
    )
    parser.add_argument(
        "--num-clients",
        type=int,
        default=10,
        help="Number of clients",
    )
    parser.add_argument(
        "--sim-total-times-s",
        type=float,
        default=seconds_per_month * 12.0,
        help="Total simulation time in seconds",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    hint_bytes = simulate_updates(
        args.db_updates_per_month,
        args.num_clients,
        args.sim_total_times_s,
    )
    num_months = args.sim_total_times_s / seconds_per_month
    print(f"Average hint bytes / month: {hint_bytes / num_months}")


if __name__ == "__main__":
    main()
