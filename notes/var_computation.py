# Validate the analytic formula for the case of block size = 1 and permutation
from __future__ import division
import math
try:
    from functools import lru_cache
except ImportError:
    from lru_cache import lru_cache
import itertools
import numpy as np
from scipy.special import binom

@lru_cache(maxsize=None)
def g(r, s, t, k):
    """Recursive function defined in the notes.
    """
    if r == 0:
        return int(t == 0)
    elif r * k < s or t < 0:
        return 0
    else:
        h = np.arange(min(s, k) + 1)
        coeff = binom(s, h) * binom(r * k - s, k - h)
        return coeff.dot(np.array([g(r - 1, s - h_, t - (h_ % 2), k) for h_ in h]))

def get_Ax_zero_probs(n, m, w, k, f):
    """Probability of A(x - x') = 0 if x - x' has Hamming weight @w, A of size @m x @n
    is first chosen to have, on each row, @k entries that are one, then all entries are flipped
    with probability @f (i.e. random walk has length w).
    """
    assert n >= m * k
    counts = np.zeros((w + 1, w + 1))
    for wprime in range(min(w, m * k) + 1):
        for q in range(min(wprime, m) + 1):
            counts[wprime, q] = binom(w, wprime) * binom(n - w, m * k - wprime) * g(m, wprime, q, k)
    counts = counts.sum(axis=0)
    normalized_counts = counts / (math.factorial(m * k) * binom(n, m * k) / math.factorial(k)**m)
    assert np.allclose(normalized_counts.sum(), 1.0)
    q = np.arange(w + 1)
    random_walk_probs = (0.5 - 0.5 * (1 - 2 * f)**w)**q * (0.5 + 0.5 * (1 - 2 * f)**w)**(m - q)
    return normalized_counts.dot(random_walk_probs)

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

# for w in range(n + 1):
# for k in range(n + 1):
    # k = 1
    # print(p_odd(n, w, k) * (0.5 - 0.5 * (1 - 2 * f)**w) + (1 - p_odd(n, w, k)) * (0.5 + 0.5 * (1 - 2 * f)**w))


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



counts_num_odd_blocks_per_wprime_analytic = np.zeros((w + 1, w + 1))
for wprime in range(min(w, m*k) + 1):
    for q in range(min(wprime, m) + 1):
        counts_num_odd_blocks_per_wprime_analytic[wprime, q] = binom(w, wprime) * binom(n - w, m*k - wprime) * g(m, wprime, q, k)
print(counts_num_odd_blocks_per_wprime_analytic / math.factorial(m*k) * math.factorial(k)**m / binom(n, m*k) -
      counts_num_odd_blocks_per_wprime / counts_num_odd_blocks_per_wprime.sum())

import sys
sys.path.append('/Users/tridao/Documents/Research/SAT_18/low_density_parity_checks')
import random_walk

f = 0.1
vector_count, Ax_zero_probs = random_walk.get_Ax_zero_probs_incompleteCol(n, m, w, k, f)
a = np.array([(count / vector_count, prob) for prob, count in Ax_zero_probs.items()])
prob = a[:, 0].dot(a[:, 1])
prob_analytic = get_Ax_zero_probs(n, m, w, k, f)
print(prob_analytic - prob)
