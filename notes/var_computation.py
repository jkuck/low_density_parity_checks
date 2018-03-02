# Validate the analytic formula for the case of block size = 1 and permutation
import itertools
import numpy as np
from scipy.special import binom

n = 8
m = 3
w = 6
index_set = {1, 2, 5, 6, 7, 8}  # Arbitrary set of indices of x that are 1

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
counts_analytic = np.empty_like(counts)
for w_prime in range(min(w + 1, m + 1)):
    counts_analytic[w_prime] = binom(w, w_prime) * binom(n - w, m - w_prime)
print(counts - counts_analytic)
