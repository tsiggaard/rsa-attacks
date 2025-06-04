# Python Implementations of the Algorithms from the Bachelor Report "Attacks Against RSA in Theory and Practice"

## By Thomas Siggaard, Emil Hasseltoft Kjær Hansen, and Johan Enevoldsen

These are our implementations of the algorithms described in our bachelor report. The code is written in Python and is intended to be used for educational purposes.

## Requirements

- Python 3.x
- Required libraries: `gmpy2`, `sympy`, `numpy`, `matplotlib`, `rsa`

You can install the required libraries using pip:

```bash
pip install gmpy2 sympy numpy matplotlib rsa
```

## Running the Fermat Factorisation Method

To simulate the Fermat factorisation method iteration count, use the following command:

```bash
python3 fermat.py
```

## Running the Quadratic Sieve Factorisation Algorithm

To run the Quadratic Sieve factorisation algorithm, you can use the following command:

```bash
python3 qs_testing.py
```

## Running the Bleichenbacher Attack

To run the Bleichenbacher attack against a local padding oracle, with a random RSA key and message, use the following command:

```bash
python3 bleichenbacher.py
```
