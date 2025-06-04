import time
from qs_impl import qs
from sympy import randprime
import numpy as np
import multiprocessing

def generate_n_p_q(n_size):
    j = n_size // 2
    k, p, q = 0, 0, 0
    while k != n_size:
        p = randprime(2**(j-1), 2**j)
        q = randprime(2**(j-1), 2**j)
        if p == q:
            continue
        n = p*q
        k = len(bin(n)[2:])
    return n, p, q

def test_qs(bit_length=20):
    # Find primes p and q of bit_length // 2
    n,p,q = generate_n_p_q(bit_length)
    print(f"Generated primes p={p}, q={q}, N=p*q={n}")

    # Factor with quadratic sieve
    time_start = time.time()
    p,q = qs(n, debug=True)
    time_end = time.time()
    print(f"Quadratic Sieve completed in {time_end - time_start:.2f} seconds.")

    # 6. Verify
    if p is None or q is None:
        print("Failed to factor N.")
    else:
        print(f"Factors found: p={p}, q={q}")
        if p * q == n:
            print("Factors verified successfully.")

def test_single_n_time(args):
    n_size, o1 = args
    n, p, q = generate_n_p_q(n_size)
    time_start = time.time()
    p, q = qs(n, o1=o1, debug=False)
    time_end = time.time()
    if p is None or q is None:
        return (0, False)
    if p and q:
        return (time_end - time_start, True)
    else:
        return (0, False)

def test_different_sizes(lower_limit, upper_limit, runs):
    print("Testing quadratic sieve with different sizes of n...")
    times = []
    while lower_limit <= upper_limit:
        o1 = 0.3557 * np.exp(-0.02377 * lower_limit)
        print(f"Testing size: {lower_limit}, o1: {o1:.10f}")
        with multiprocessing.Pool(processes=8) as pool:
            results = pool.map(test_single_n_time, [(lower_limit, o1)] * runs)
        failure_count = 0
        success_times = []

        for r in results:
            if r == 0:
                failure_count += 1
            else:
                elapsed, did_factor = r
                if not did_factor:
                    failure_count += 1
                else:
                    success_times.append(elapsed)
        failure_rate = failure_count / runs
        if success_times:
            avg_time = np.mean(success_times)
            print(f"  • Failure rate: {(failure_rate*100):.2f}% ({failure_count}/{runs})")
            print(f"  • Average time (successful runs only): {avg_time:.10f} seconds")
            print(f"({lower_limit}, {avg_time:.10f})\n")
            times.append(avg_time)
        else:
            print(f"  • Failure rate: {failure_rate:.4f} ({failure_count}/{runs})")
            print("  • All runs failed, so no average time.\n")
            times.append(None)
        lower_limit += 2
    for i, t in enumerate(times):
        print(f"({20 + 2 * i}, {t:.10f})")

if __name__ == '__main__':
    runs = 200
    lower_limit = 20
    upper_limit = 40
    test_different_sizes(lower_limit, upper_limit, runs)