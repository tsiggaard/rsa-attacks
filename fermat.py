import math
import time
import numpy as np
import matplotlib.pyplot as plt
import gmpy2
from sympy import randprime, nextprime
import multiprocessing

def fermat(n):
    """
    Fermat's factorization method.
    """
    counter = 0
    x = math.ceil(math.sqrt(n))
    y2 = x*x - n
    while not math.sqrt(y2).is_integer():
        counter += 1
        x += 1
        y2 = x*x - n
    y = int(math.sqrt(y2))
    p = x - y
    q = x + y
    return p, q, counter

def fermat_isqrt(n):
    """
    Fermat's factorization method using integer square root.
    """
    counter = 0
    x = math.isqrt(n) + 1
    y2 = x*x - n
    while not math.isqrt(y2)**2 == y2:
        counter += 1
        x += 1
        y2 = x*x - n
    y = int(math.isqrt(y2))
    p = x - y
    q = x + y
    return p, q, counter

def fermat_GMPY2(n):
    """
    Fermat's factorization method using gmpy2 for high precision.
    """
    counter = 0
    x = gmpy2.isqrt(gmpy2.mpz(n)) + 1
    y2 = x*x - n
    y = gmpy2.isqrt(gmpy2.mpz(y2))
    while not y*y == y2:
        counter += 1
        x += 1
        y2 = x*x - n 
        y = gmpy2.isqrt(gmpy2.mpz(y2))
    y = int(gmpy2.sqrt(y2))
    p = x - y
    q = x + y
    return p, q, counter

def generate_n_p_q(n_size):
    """
    Generate a random number n of size n_size, and two primes p and q such that n = p * q.
    Returns the number n, and the primes p and q.
    """
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

def fermat_iterations(n, p, q):
    """
    Calculate the number of iterations for Fermat's factorization method.
    """
    return (p + q) / 2 - math.ceil(math.sqrt(n))

def fermat_iterations_upper_bound(p, q):
    """
    Calculate the upper bound for the number of iterations for Fermat's factorization method.
    """
    return abs((q - p) / 2)

def test_single_n(fermat_function, n):
    """
    Test a single number n with the given Fermat function.
    Returns the factors p and q, and the time taken for the factorization.
    """
    start_time = time.time()
    p, q, counter = fermat_function(n)
    end_time = time.time()
    k = len(bin(n)[2:])
    print(f"n: {n}, k: {k}, p: {p}, q: {q}, counter: {counter}, duration: {end_time - start_time} seconds using {fermat_function.__name__}")
    return p, q

def get_key_lists(start, end, runs=1000):
    """
    Generate a list of lists of random numbers n, where each n is a product of two primes p and q.
    """
    n_lists = []
    for i in range(start, end, 2):
        n_list = []
        for _ in range(runs):
            n, p, q = generate_n_p_q(i)
            n_list.append([n, p, q])
        n_lists.append(n_list)
    return n_lists

def process_single_n(args):
    """
    Process a single number n with the given Fermat function in a separate process.
    Returns the time taken for the factorization and the number of iterations (counter).
    """
    fermat_function, n = args
    start_time = time.time_ns()
    p, q, counter = fermat_function(n)
    if p * q != n:
        return 0, 0  # Return 0 duration if the factorization is incorrect
    end_time = time.time_ns()
    return end_time - start_time, counter

def test_many_n(fermat_function, start, end, key_lists, show_table=False):
    """
    Test multiple numbers n with the given Fermat function in parallel.
    Returns the average time taken for the factorization and the number of iterations (counter).
    """
    time_list = []
    counter_list = []
    total_start_time = time.time()
    for key_list in key_lists:
        n_list = [n[0] for n in key_list]
        with multiprocessing.Pool(processes=8) as pool:
            res = pool.map(process_single_n, [(fermat_function, n) for n in n_list])
            print(f"Processed {len(n_list)} numbers of size {len(bin(n_list[0])[2:])} with {fermat_function.__name__} in {time.time() - total_start_time} seconds")
        durations, counters = zip(*res)  # Unzip the results
        # Filter out incorrect factorizations
        durations = [d for d in durations if d != 0]  
        counters = [c for c in counters if c != 0]
        time_list.append(sum(durations) / len(durations) if durations else 0)
        counters = sum(counters) / len(counters) if counters else 0
        counter_list.append(counters)
    total_end_time = time.time()
    print(f"Total time: {total_end_time - total_start_time} seconds using {fermat_function.__name__}")

    if show_table:
        print("Bitlen".ljust(10) + "Counter".ljust(20) + "Average Time (s)".ljust(30))
        for i in range(len(counter_list)):
            print(str(2*i+start).ljust(10) + str(counter_list[i]).ljust(20) + str(time_list[i]/1e9).ljust(30))   

    x = np.arange(start, end, 2)
    y = np.array(time_list) / 1e9  # Convert to seconds
    print(np.exp(np.polyfit(x, np.log(y), 1)))
    z = np.array(counter_list)
    return x, y, z

def check_many_n(start, end, runs=1000):
    key_lists = get_key_lists(start, end, runs)
    print("Bitlen, Average Counter, Average Upper Bound")
    average_counter = []
    average_upper_bound = []
    x = np.arange(start, end, 2)
    for i in range(len(key_lists)):
        key_list = key_lists[i]
        counter_list = [fermat_iterations(n, p, q) for n, p, q in key_list] 
        upper_bound_list = [fermat_iterations_upper_bound(p, q) for _, p, q in key_list]
        print("(" + str(i*2 + start) + ", " + str(sum(counter_list)/len(counter_list)) +
              ", " + str(sum(upper_bound_list)/len(upper_bound_list)) + ")")
        average_counter.append(sum(counter_list) / len(counter_list))
        average_upper_bound.append(sum(upper_bound_list) / len(upper_bound_list))
    print(np.exp(np.polyfit(x, np.log(average_counter), 1)))
    print(np.exp(np.polyfit(x, np.log(average_upper_bound), 1)))

def show_time_plot(x, y, runs, low, high):
    methods = ["Fermat", "Isqrt", "GMPY2"]
    for i in range(len(x)):
        plt.plot(x[i], y[i], label=methods[i])
    plt.legend()
    plt.xticks(np.arange(low, high, 2))
    plt.title(f"Fermat's factorization method - {runs} runs")
    plt.xlabel("Bit length")
    plt.ylabel("Average time (s)")
    plt.show()

def benchmark_fermat_methods(low, high, runs=1000):
    """
    Benchmark the Fermat methods and print/plot the results.
    """
    key_lists = get_key_lists(low, high, runs)
    x1, y1, z1 = test_many_n(fermat, low, high, key_lists, show_table=True)
    x2, y2, z2 = test_many_n(fermat_isqrt, low, high, key_lists, show_table=True)
    x3, y3, z3 = test_many_n(fermat_GMPY2, low, high, key_lists, show_table=True)
    x = [x1, x2, x3]
    y = [y1, y2, y3]
    z = [z1, z2, z3]
    show_time_plot(x, y, runs, low, high)

def test_fermat_method(low, high, runs=1000):
    key_lists = get_key_lists(low, high, runs)
    x, y, z = test_many_n(fermat_isqrt, low, high, key_lists, show_table=True)
    show_time_plot([x], [y], runs, low, high)

def show_upper_bound_plot(iterations=1000):
    """
    Show the upper bound plot for Fermat's factorization method.
    """
    p = 3
    q = nextprime(p)
    ns = []
    iterations = []
    upper_bound = []
    for i in range(1, iterations):
        n = p * q
        ns.append(n)
        iterations.append(fermat_iterations(n, p, q))
        upper_bound.append(fermat_iterations_upper_bound(p, q))
        q = nextprime(q)
    plt.plot(ns, iterations, label="Iterations")
    plt.plot(ns, upper_bound, label="Upper bound")
    plt.xlabel("n")
    plt.ylabel("Iterations")
    plt.title("Fermat's factorization method iterations")
    plt.legend()
    plt.show()
    print(np.polyfit(ns, iterations, 1))

if __name__ == "__main__":
    # benchmark_fermat_methods(10, 51, 1000)
    # test_fermat_method(20, 61, 1000)
    check_many_n(10, 101, 10000)
    # show_upper_bound_plot(1000)