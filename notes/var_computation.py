# Validate the analytic formula for the case of block size = 1 and permutation
import math
from functools import lru_cache
import itertools
import numpy as np
from scipy.special import binom

n = 9
m = 3
w = 6
k = 1
index_set = set(np.random.choice(range(1, n+1), size=w, replace=False))  # Arbitrary set of indices of x that are 1

# k = 1, sampling without replacement
counts = np.zeros(w + 1, dtype=int)
for c in itertools.combinations(range(1, n+1), m):
    num_1_in_first_w = sum(i <= w for i in c)
    counts[num_1_in_first_w] += 1

# Take order into account. The results is the same as above where we don't care about order.
counts_order = np.zeros(w + 1, dtype=int)
for c in itertools.permutations(range(1, n+1), m):
    num_1_in_first_w = sum(i in index_set for i in c)
    counts_order[num_1_in_first_w] += 1
print(counts / counts.sum() - counts_order / counts_order.sum())

# Use analytic formula instead of counting
w_prime = np.arange(w + 1)
counts_analytic = binom(w, w_prime) * binom(n - w, m - w_prime)
print(counts - counts_analytic)

# k > 1, sampling with replacement
def p_odd(n, w, k):
    t = np.arange(1, w + 1, 2)
    counts_odd = binom(w, t) * binom(n - w, k - t)
    return counts_odd.sum() / binom(n, k)

for k in range(n + 1):
    print(p_odd(n, w, k))


# k > 1, sampling without replacement
k = 2
counts_num_odd_blocks = np.zeros(w + 1, dtype=int)
counts_num_odd_blocks_per_wprime = np.zeros((w + 1, w + 1), dtype=int)
for c in itertools.permutations(range(1, n+1), m * k):
    wprime = sum(i in index_set for i in c)
    blocks = np.array(c).reshape(m, k)
    block_intersect_sizes = [sum(b in index_set for b in block) for block in blocks]
    num_odd_blocks = sum(size % 2 for size in block_intersect_sizes)
    counts_num_odd_blocks[num_odd_blocks] += 1
    counts_num_odd_blocks_per_wprime[wprime, num_odd_blocks] += 1

@lru_cache(maxsize=None)
def f(r, s, t, k):
    if r == 0:
        return int(t == 0)
    elif r * k < s:
        return 0
    else:
        total = 0
        for h in range(min(s, k) + 1):
            if t == 0 and (h % 2):
                continue
            coeff = binom(s, h) * binom(r * k - s, k - h)
            func = f(r - 1, s - h, t - (h % 2), k)
            total += coeff * func
        return total

counts_num_odd_blocks_per_wprime_analytic = np.zeros((w + 1, w + 1), dtype=int)
for wprime in range(min(w, m*k) + 1):
    for q in range(min(wprime, m) + 1):
        counts_num_odd_blocks_per_wprime_analytic[wprime, q] = binom(w, wprime) * binom(n - w, m*k - wprime) * f(m, wprime, q, k)
print(counts_num_odd_blocks_per_wprime_analytic / math.factorial(m*k) * math.factorial(k)**m / binom(n, m*k) -
      counts_num_odd_blocks_per_wprime / counts_num_odd_blocks_per_wprime.sum())
