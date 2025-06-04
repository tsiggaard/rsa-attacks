import math
from sympy import randprime
import numpy as np
import time

def qs(N, o1, debug = False):
    if not o1:
        o1 = 0.3557 * np.exp(-0.02377 * N.bit_length())
    if debug:
        print("Starting Quadratic Sieve with N =", N)
    if N < 2:
        raise ValueError("n must be greater than 1.")
    if N % 2 == 0:
        return 2, N // 2
    sqrt_n = math.isqrt(N)
    if sqrt_n * sqrt_n == N:
        return sqrt_n, sqrt_n

    # Step 1: Choose a smoothness bound B (p. 78)
    if debug:
        print("Step 1: Computing smoothness bound B...")
    B = compute_smoothness_bound(N, o1)
    if debug:
        print(f"Step 1 finished, smoothness bound B = {B}")
    
    # Make relations
    start_number = sqrt_n + 1

    # Step 2: Find all primes up to the smoothness bound
    if debug:
        print("Step 2: Computing factor base...")
    factor_base = compute_factor_base(N,B)
    if debug:
        print("Step 2 finished")
    max_steps = compute_sieve_steps(N, o1)

    # Step 3: Find smooth relations
    if debug:
        print("Step 3: Finding smooth relations...")
    smooth_relations, relation_indexes = find_smooth_relations(N, B, factor_base, start_number, max_steps, debug)
    if smooth_relations is None or relation_indexes is None:
        if debug:
            print("Failed to find smooth relations.")
        return None, None
    if debug:
        print(f"Step 3 finished, found {len(smooth_relations)} smooth relations.")

    # Step 4: Build the exponent matrix
    if debug:
        print("Step 4: Building exponent matrix...")
    M = np.array([vec for _, vec in smooth_relations], dtype=int)
    if debug:
        print("Step 4 finished")
    
    # Step 5: Find the null space of the exponent matrix
    if debug:
        print("Step 5: Finding null space of the exponent matrix...")
    basis = nullspace_mod2(M)
    if debug:
        print("Step 5 finished")

    if basis == []:
        print("No nullspace found.")
        return None, None
    
    # Step 6: Find non-trivial factors
    if debug:
        print("Step 6: Finding non-trivial factors...")
    factor1, factor2 = find_nontrivial_factors(basis, factor_base, N, relation_indexes, start_number, M)
    if debug:
        print("Step 6 finished")
    return factor1, factor2

def find_nontrivial_factors(basis, factor_base, N, relation_indexes, start_number, M):
    for c in basis:
        # Compute a
        a = 1
        indexes = [i for i, v in enumerate(c) if int(v) % 2 == 1]
        for i in indexes:
            index = relation_indexes[i]
            x = start_number + index
            a = (a * x) % N

        # Compute b
        exponents_sum = np.zeros(len(factor_base), dtype=int)
        for i in indexes:
            v = M[i]
            for j in range(len(v)):
                exponents_sum[j] += v[j]

        b = 1
        for p, e in zip(factor_base, exponents_sum):
            b = (b * pow(int(p), int(e) // 2, N)) % N

        # Verify a and b
        if a % N == b % N or -a % N == b % N:
            continue
        
        # Compute the GCD
        factor1 = math.gcd(a - b, N)
        factor2 = math.gcd(a + b, N)
        if factor1 == 1 or factor1 == N or factor2 == 1 or factor2 == N or factor1 == factor2:
            continue
        return factor1, factor2
    return None, None

def nullspace_mod2(A):
    A = (A % 2).T
    A = A.copy() % 2
    m, n = A.shape
    pivots, row = [], 0
    for col in range(n):
        pivot = next((r for r in range(row, m) if A[r, col]), None)
        if pivot is None:
            continue
        A[[row, pivot]] = A[[pivot, row]]
        for r in range(m):
            if r != row and A[r, col]:
                A[r] ^= A[row]
        pivots.append(col)
        row += 1
    free_cols = [c for c in range(n) if c not in pivots]
    basis = []
    for f in free_cols:
        v = np.zeros(n, dtype=int)
        v[f] = 1
        for i, col in enumerate(pivots):
            if A[i, f]:
                v[col] = 1
        basis.append(v % 2)
    return basis

def compute_factor_base(N, B):
    sieve = [True] * (B + 1)
    sieve[0:2] = [False, False]
    for i in range(2, int(B**0.5) + 1):
        if sieve[i]:
            sieve[i*i:B+1:i] = [False] * len(range(i*i, B + 1, i))

    factor_base = [2] if N % 2 != 0 else []

    for p in range(3, B + 1, 2):
        if not sieve[p]:
            continue
        if pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)

    return factor_base

def factorize(n, primes):
    factors = []
    for p in primes:
        if p*p > n:
            break
        while n % p == 0:
            factors.append(p)
            n //= p
    if n > 1:
        factors.append(n)
    return factors

def find_smooth_relations(N, B, factor_base, start_number, max_steps, debug = False):
    relations = []
    current_xi = start_number
    for x in range(start_number, start_number + len(factor_base) + 10):
        number = x * x - N
        relations.append(number)
        current_xi = x
    smooth_relations = []
    relation_indexes = []
    index = 0
    steps = 0
    while len(smooth_relations) <= len(factor_base):
        num = relations[index]
        exponent_vector = [0] * len(factor_base)
        factors = factorize(num, factor_base)
        if all(p <= B for p in factors):
            for factor in factors:
                if factor in factor_base:
                    factor_index = factor_base.index(factor)
                    exponent_vector[factor_index] += 1

            smooth_relations.append([num, exponent_vector])
            relation_indexes.append(index)
        else:
            x = current_xi + 1
            number = x * x - N
            relations.append(number)
            current_xi = x
        index += 1
        steps += 1
        if steps > max_steps:
            if debug:
                print("Max steps reached, stopping search for smooth relations.")
            return None, None
    return smooth_relations, relation_indexes

def compute_sieve_steps(n, o1):
    log_n = math.log(n)
    log_log_n = math.log(log_n)
    exp = (1 + o1) * math.sqrt(log_n * log_log_n)
    return int(math.exp(exp))

def compute_smoothness_bound(n, o1):
    log_n = math.log(n)
    log_log_n = math.log(log_n)
    exp = (0.5 + o1) * math.sqrt(log_n * log_log_n)
    B = int(math.exp(exp))
    return B